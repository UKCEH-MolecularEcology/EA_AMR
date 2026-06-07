"""
Author: Susheel Bhanu BUSI
Affiliation: Molecular Ecology group, UKCEH
Date: [2024-10-16]
Run: snakemake -s workflow/rules/mge_metal.smk --use-conda --cores 64 -rp
Latest modification:
Purpose: Analyze ARGs alongside metal resistance genes (BacMet2) and mobile genetic
         elements (ISEScan) to investigate co-occurrence and potential co-selection.
         Addresses Reviewer 1 comment #6 — ARGs and metal resistance genes co-occur
         due to shared mode of action and may be co-selected and spread together via MGEs.
"""


############################################
rule mge_metal:
    input:
        os.path.join(RESULTS_DIR, "metal_resistance/combined_metal_resistance.tsv"),
        os.path.join(RESULTS_DIR, "mge_cooccurrence/combined_mge_arg_summary.tsv"),
        os.path.join(RESULTS_DIR, "genomad/combined_genomad_classification.tsv")
    output:
        touch("status/mge_metal.done")


############################################
localrules: download_bacmet_db


############################################
# Download and index the BacMet2 experimental database
rule download_bacmet_db:
    output:
        fasta=os.path.join(config.get("bacmet", {}).get("db_dir", "/hdd0/susbus/databases/bacmet"), "BacMet2_EXP_database.fasta"),
        dmnd=os.path.join(config.get("bacmet", {}).get("db_dir", "/hdd0/susbus/databases/bacmet"), "BacMet2_EXP_database.dmnd")
    log:
        os.path.join(RESULTS_DIR, "logs/download_bacmet_db.log")
    params:
        url=config.get("bacmet", {}).get("db_url", "http://bacmet.biomedicine.gu.se/download/BacMet2_EXP.fasta"),
        db_dir=config.get("bacmet", {}).get("db_dir", "/hdd0/susbus/databases/bacmet")
    conda:
        os.path.join(ENV_DIR, "bacmet.yaml")
    message:
        "Setup: downloading and indexing BacMet2 experimental database to {params.db_dir}"
    shell:
        "(date && "
        "mkdir -p {params.db_dir} && "
        "wget -O {output.fasta} {params.url} --no-check-certificate && "
        "diamond makedb --in {output.fasta} -d {output.dmnd} && "
        "date) &> >(tee {log})"


############################################
# Run DIAMOND blastp against BacMet2 per sample
rule metal_resistance_diamond:
    input:
        faa=os.path.join(RESULTS_DIR, "proteins/{sid}.faa"),
        db=os.path.join(config.get("bacmet", {}).get("db_dir", "/hdd0/susbus/databases/bacmet"), "BacMet2_EXP_database.dmnd")
    output:
        os.path.join(RESULTS_DIR, "metal_resistance/{sid}/{sid}_bacmet.tsv")
    log:
        os.path.join(RESULTS_DIR, "logs/metal_resistance.{sid}.log")
    threads:
        config["mmseqs2"]["threads"]
    params:
        pident=config.get("bacmet", {}).get("pident", 80),
        qcov=config.get("bacmet", {}).get("qcov", 80)
    conda:
        os.path.join(ENV_DIR, "bacmet.yaml")
    wildcard_constraints:
        sid="|".join(SAMPLES.index)
    message:
        "Metal resistance: DIAMOND BacMet2 search for {wildcards.sid}"
    shell:
        "(date && "
        "mkdir -p $(dirname {output}) && "
        "diamond blastp "
        "-q {input.faa} "
        "-d {input.db} "
        "-o {output} "
        "--threads {threads} "
        "-f 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore "
        "--id {params.pident} "
        "--query-cover {params.qcov} "
        "--more-sensitive && "
        "date) &> >(tee {log})"


############################################
# Summarise metal resistance hits across all samples
rule summarise_metal_resistance:
    input:
        expand(os.path.join(RESULTS_DIR, "metal_resistance/{sid}/{sid}_bacmet.tsv"), sid=SAMPLES.index)
    output:
        os.path.join(RESULTS_DIR, "metal_resistance/combined_metal_resistance.tsv")
    log:
        os.path.join(RESULTS_DIR, "logs/summarise_metal_resistance.log")
    message:
        "Metal resistance: summarising BacMet2 results across all samples"
    run:
        import os
        import pandas as pd

        COLS = ["qseqid", "sseqid", "pident", "len", "mismatch", "gap",
                "qs", "qe", "ss", "se", "evalue", "bitscore"]

        # Lookup dict: gene-name prefix -> metal/biocide class
        CLASS_LOOKUP = {
            "merA": "mercury", "merB": "mercury", "merC": "mercury",
            "merD": "mercury", "merE": "mercury", "merP": "mercury",
            "merT": "mercury",
            "copA": "copper",  "copB": "copper",  "copC": "copper",
            "cusA": "copper",  "cusB": "copper",
            "arsA": "arsenic", "arsB": "arsenic", "arsC": "arsenic",
            "arsD": "arsenic",
            "cadA": "cadmium", "cadC": "cadmium",
            "chrA": "chromate", "chrB": "chromate",
        }

        def assign_class(sseqid):
            for prefix, metal_class in CLASS_LOOKUP.items():
                if prefix.lower() in sseqid.lower():
                    return metal_class
            return "other"

        all_dfs = []
        for fp in input:
            sid = os.path.basename(fp).replace("_bacmet.tsv", "")
            try:
                df = pd.read_csv(fp, sep="\t", header=None, names=COLS)
                df["sid"] = sid
                df["metal_class"] = df["sseqid"].apply(assign_class)
                all_dfs.append(df)
            except (pd.errors.EmptyDataError, FileNotFoundError):
                pass

        if all_dfs:
            combined = pd.concat(all_dfs, ignore_index=True)
        else:
            combined = pd.DataFrame(columns=COLS + ["sid", "metal_class"])

        combined.to_csv(output[0], sep="\t", index=False)

        # Write per-sample x metal_class count pivot as a comment-free summary TSV
        summary_path = os.path.join(
            os.path.dirname(output[0]), "metal_resistance_summary.tsv"
        )
        if not combined.empty:
            pivot = (
                combined.groupby(["sid", "metal_class"])
                .size()
                .unstack(fill_value=0)
                .reset_index()
            )
            pivot.to_csv(summary_path, sep="\t", index=False)
        else:
            pd.DataFrame(columns=["sid"]).to_csv(summary_path, sep="\t", index=False)


############################################
# Run ISEScan per sample to detect insertion sequences
rule isescan_mge:
    input:
        fasta=os.path.join(RESULTS_DIR, "assembly_filtered/{sid}/{sid}_noOrganellar.fasta")
    output:
        outdir=directory(os.path.join(RESULTS_DIR, "isescan/{sid}")),
        # ISEScan mirrors the input directory structure inside outdir, so the
        # .sum file lands under assembly_filtered/ (the parent dir of the input fasta).
        summ=os.path.join(RESULTS_DIR, "isescan/{sid}/assembly_filtered/{sid}/{sid}_noOrganellar.fasta.sum")
    priority: -1
    log:
        os.path.join(RESULTS_DIR, "logs/isescan.{sid}.log")
    threads: 32
    conda:
        os.path.join(ENV_DIR, "isescan.yaml")
    wildcard_constraints:
        sid="|".join(SAMPLES.index)
    message:
        "MGE: ISEScan for {wildcards.sid}"
    shell:
        "(date && "
        "mkdir -p {output.outdir} && "
        "cd {output.outdir} && "
        "isescan.py --seqfile {input.fasta} --output {output.outdir} --nthread {threads} && "
        "date) &> >(tee {log})"


############################################
# Co-occurrence analysis: ISEs + geNomad plasmid/virus + ARGs on the same contig
# geNomad provides the broadest MGE classification (plasmids + phages);
# ISEScan adds insertion-sequence resolution within chromosomal contigs.
rule mge_arg_cooccurrence:
    input:
        is_sum=os.path.join(RESULTS_DIR, "isescan/{sid}/assembly_filtered/{sid}/{sid}_noOrganellar.fasta.sum"),
        rgi=os.path.join(RESULTS_DIR, "amr/{sid}/{sid}_rgi.txt"),
        genomad=os.path.join(
            RESULTS_DIR,
            "genomad/{sid}/{sid}_aggregated_classification/{sid}_aggregated_classification.tsv"
        )
    output:
        os.path.join(RESULTS_DIR, "mge_cooccurrence/{sid}/{sid}_mge_arg_contigs.tsv")
    log:
        os.path.join(RESULTS_DIR, "logs/mge_arg_cooccurrence.{sid}.log")
    wildcard_constraints:
        sid="|".join(SAMPLES.index)
    message:
        "MGE-ARG co-occurrence (ISEScan + geNomad): {wildcards.sid}"
    run:
        import os
        import pandas as pd

        # ── ISEScan insertion sequences ───────────────────────────────────────
        try:
            is_df = pd.read_csv(input.is_sum, sep=r"\s+", comment="#", engine="python")
            is_contigs = set(is_df["seqID"].astype(str)) if "seqID" in is_df.columns else set()
            is_families = (
                is_df.groupby("seqID")["isFamily"]
                .apply(lambda x: ",".join(sorted(set(x.astype(str)))))
                .to_dict()
            ) if "seqID" in is_df.columns and "isFamily" in is_df.columns else {}
        except (pd.errors.EmptyDataError, FileNotFoundError, KeyError):
            is_contigs, is_families = set(), {}

        # ── geNomad classification ────────────────────────────────────────────
        try:
            gd = pd.read_csv(input.genomad, sep="\t")
            seq_col = "seq_name" if "seq_name" in gd.columns else "sequence"
            plasmid_contigs = set(gd.loc[gd["classification"] == "Plasmid", seq_col].astype(str))
            virus_contigs   = set(gd.loc[gd["classification"] == "Virus",   seq_col].astype(str))
            gd_class  = gd.set_index(seq_col)["classification"].to_dict()
            gd_plas_s = gd.set_index(seq_col).get("plasmid_score", pd.Series(dtype=float)).to_dict()
            gd_vir_s  = gd.set_index(seq_col).get("virus_score",   pd.Series(dtype=float)).to_dict()
        except (pd.errors.EmptyDataError, FileNotFoundError, KeyError):
            plasmid_contigs = virus_contigs = set()
            gd_class = gd_plas_s = gd_vir_s = {}

        # ── RGI ARG annotations ───────────────────────────────────────────────
        try:
            rgi_df = pd.read_csv(input.rgi, sep="\t")
            if "Contig" in rgi_df.columns:
                rgi_df["contig_id"] = rgi_df["Contig"].astype(str).str.rsplit("_", n=1).str[0]
            else:
                rgi_df["contig_id"] = pd.Series(dtype=str)
            arg_contigs = set(rgi_df["contig_id"])
            arg_names = (
                rgi_df.groupby("contig_id")["Best_Hit_ARO"]
                .apply(lambda x: ",".join(sorted(set(x.astype(str)))))
                .to_dict()
            ) if "Best_Hit_ARO" in rgi_df.columns else {}
        except (pd.errors.EmptyDataError, FileNotFoundError, KeyError):
            arg_contigs, arg_names = set(), {}

        # A contig is "on_MGE" if ISEScan found an IS element, or geNomad
        # classified it as plasmid or virus.
        mge_contigs = is_contigs | plasmid_contigs | virus_contigs
        all_contigs = mge_contigs | arg_contigs

        rows = []
        for contig in sorted(all_contigs):
            rows.append({
                "contig_id":              contig,
                "has_IS":                 contig in is_contigs,
                "has_plasmid_genomad":    contig in plasmid_contigs,
                "has_virus_genomad":      contig in virus_contigs,
                "on_MGE":                 contig in mge_contigs,
                "has_ARG":                contig in arg_contigs,
                "IS_families":            is_families.get(contig, ""),
                "genomad_classification": gd_class.get(contig, "chromosome"),
                "genomad_plasmid_score":  gd_plas_s.get(contig, ""),
                "genomad_virus_score":    gd_vir_s.get(contig, ""),
                "ARG_names":              arg_names.get(contig, ""),
            })

        out_df = pd.DataFrame(rows, columns=[
            "contig_id", "has_IS", "has_plasmid_genomad", "has_virus_genomad",
            "on_MGE", "has_ARG", "IS_families",
            "genomad_classification", "genomad_plasmid_score", "genomad_virus_score",
            "ARG_names"
        ])
        os.makedirs(os.path.dirname(output[0]), exist_ok=True)
        out_df.to_csv(output[0], sep="\t", index=False)


############################################
# Summarise MGE-ARG co-occurrence across all samples
rule summarise_mge_arg:
    input:
        expand(os.path.join(RESULTS_DIR, "mge_cooccurrence/{sid}/{sid}_mge_arg_contigs.tsv"), sid=SAMPLES.index)
    output:
        os.path.join(RESULTS_DIR, "mge_cooccurrence/combined_mge_arg_summary.tsv")
    log:
        os.path.join(RESULTS_DIR, "logs/summarise_mge_arg.log")
    message:
        "MGE-ARG co-occurrence: summarising across all samples"
    run:
        import os
        import pandas as pd

        all_dfs = []
        for fp in input:
            sid = os.path.basename(fp).replace("_mge_arg_contigs.tsv", "")
            try:
                df = pd.read_csv(fp, sep="\t")
                df["sample"] = sid
                all_dfs.append(df)
            except (pd.errors.EmptyDataError, FileNotFoundError):
                pass

        if all_dfs:
            combined = pd.concat(all_dfs, ignore_index=True)
        else:
            combined = pd.DataFrame(columns=["contig_id", "has_IS", "has_ARG", "IS_families", "ARG_names", "sample"])

        # Per-sample summary statistics
        summary_rows = []
        for sid, grp in combined.groupby("sample"):
            n_is  = int(grp["has_IS"].sum())
            n_arg = int(grp["has_ARG"].sum())
            n_both = int((grp["has_IS"] & grp["has_ARG"]).sum())
            pct = (n_both / n_arg * 100) if n_arg > 0 else 0.0
            summary_rows.append({
                "sample":                     sid,
                "n_contigs_with_IS":          n_is,
                "n_contigs_with_ARG":         n_arg,
                "n_contigs_with_both":        n_both,
                "pct_ARG_contigs_with_IS":    round(pct, 2),
            })

        summary_df = pd.DataFrame(summary_rows, columns=[
            "sample", "n_contigs_with_IS", "n_contigs_with_ARG",
            "n_contigs_with_both", "pct_ARG_contigs_with_IS"
        ])

        # Write combined per-contig table
        combined.to_csv(output[0], sep="\t", index=False)

        # Write per-sample summary alongside
        summary_path = os.path.join(
            os.path.dirname(output[0]), "mge_arg_per_sample_summary.tsv"
        )
        summary_df.to_csv(summary_path, sep="\t", index=False)
