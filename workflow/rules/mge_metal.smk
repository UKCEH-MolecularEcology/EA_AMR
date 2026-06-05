"""
Author: Susheel Bhanu BUSI
Affiliation: Molecular Ecology group, UKCEH
Date: [2024-10-16]
Run: snakemake -s workflow/rules/mge_metal.smk --use-conda --cores 64 -rp
Latest modification:
Purpose: Analyse mobile genetic elements (MGEs), metal resistance genes, and their
         co-occurrence with ARGs — addressing Reviewer 1 comment #6:
         "it is necessary to analyze ARGs along with metal resistance genes and mobile
         genetic elements because ARGs and metal resistance genes co-occur due to shared
         mode of action, which could be co-selected and spread together."
         Tools:
           - ISEScan  : detection of insertion sequences (MGEs) in assemblies
           - BacMet2  : curated database of metal/biocide resistance genes (DIAMOND search)
           - Python   : contig-level co-occurrence analysis of IS elements and ARGs

Config additions required (config.yaml):
  bacmet:
    db_url: "http://bacmet.biomedicine.gu.se/download/BacMet2_EXP.fasta"
    pident: 80
    qcov: 80
    threads: 24
"""


############################################
# Aggregate / status rule — always listed first
rule mge_metal:
    input:
        expand(os.path.join(RESULTS_DIR, "metal_resistance/{sid}/{sid}_bacmet.tsv"), sid=SAMPLES.index),
        expand(os.path.join(RESULTS_DIR, "isescan/{sid}/{sid}.fasta.sum"), sid=SAMPLES.index),
        expand(os.path.join(RESULTS_DIR, "mge_cooccurrence/{sid}/{sid}_mge_arg_contigs.tsv"), sid=SAMPLES.index),
        os.path.join(RESULTS_DIR, "mge_cooccurrence/combined_mge_arg_summary.tsv"),
        os.path.join(RESULTS_DIR, "metal_resistance/combined_metal_resistance.tsv")
    output:
        touch("status/mge_metal.done")


############################################
localrules: download_bacmet_db


############################################
# Rule 1: Download BacMet2 experimentally verified resistance gene sequences and build DIAMOND db
rule download_bacmet_db:
    output:
        fasta=os.path.join(DB_DIR, "bacmet/BacMet2_EXP.fasta"),
        dmnd=os.path.join(DB_DIR, "bacmet/BacMet2_EXP.dmnd")
    log:
        os.path.join(RESULTS_DIR, "logs/download_bacmet_db.log")
    params:
        db_url=config["bacmet"]["db_url"],
        db_dir=os.path.join(DB_DIR, "bacmet")
    conda:
        os.path.join(ENV_DIR, "bacmet.yaml")
    message:
        "Setup: Download BacMet2 database and build DIAMOND index"
    shell:
        "(date && "
        "mkdir -p {params.db_dir} && "
        "wget -O {output.fasta} {params.db_url} --no-check-certificate && "
        "diamond makedb --in {output.fasta} -d {params.db_dir}/BacMet2_EXP && "
        "date) &> >(tee {log})"


############################################
# Rule 2: Search predicted proteins against BacMet2 using DIAMOND
rule metal_resistance_diamond:
    input:
        faa=os.path.join(RESULTS_DIR, "prodigal/{sid}/{sid}.faa"),
        db=os.path.join(DB_DIR, "bacmet/BacMet2_EXP.dmnd")
    output:
        os.path.join(RESULTS_DIR, "metal_resistance/{sid}/{sid}_bacmet.tsv")
    log:
        os.path.join(RESULTS_DIR, "logs/metal_resistance.{sid}.log")
    threads:
        config["mmseqs2"]["threads"]
    params:
        pident=config["bacmet"]["pident"],
        qcov=config["bacmet"]["qcov"]
    conda:
        os.path.join(ENV_DIR, "bacmet.yaml")
    wildcard_constraints:
        sid="|".join(SAMPLES.index)
    message:
        "Step: {wildcards.sid} — DIAMOND BacMet2 metal resistance search"
    shell:
        "(date && "
        "mkdir -p $(dirname {output}) && "
        "diamond blastp "
        "  -q {input.faa} "
        "  -d {input.db} "
        "  -o {output} "
        "  --threads {threads} "
        "  -f 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore "
        "  --id {params.pident} "
        "  --query-cover {params.qcov} "
        "  --more-sensitive && "
        "date) &> >(tee {log})"


############################################
# Rule 3: Detect insertion sequences (MGEs) with ISEScan
rule isescan_mge:
    input:
        os.path.join(RESULTS_DIR, "assembly/{sid}/{sid}.fasta")
    output:
        dir=directory(os.path.join(RESULTS_DIR, "isescan/{sid}/")),
        summary=os.path.join(RESULTS_DIR, "isescan/{sid}/{sid}.fasta.sum")
    log:
        os.path.join(RESULTS_DIR, "logs/isescan.{sid}.log")
    threads:
        8
    conda:
        os.path.join(ENV_DIR, "isescan.yaml")
    wildcard_constraints:
        sid="|".join(SAMPLES.index)
    message:
        "Step: {wildcards.sid} — ISEScan mobile genetic element detection"
    shell:
        "(date && "
        "mkdir -p {output.dir} && "
        "isescan.py "
        "  --seqfile {input} "
        "  --output {output.dir} "
        "  --nthread {threads} && "
        "date) &> >(tee {log})"


############################################
# Rule 4: Per-sample contig-level co-occurrence of IS elements and ARGs
rule mge_arg_cooccurrence:
    input:
        isescan_sum=os.path.join(RESULTS_DIR, "isescan/{sid}/{sid}.fasta.sum"),
        rgi_txt=os.path.join(RESULTS_DIR, "amr/{sid}/{sid}_rgi.txt")
    output:
        os.path.join(RESULTS_DIR, "mge_cooccurrence/{sid}/{sid}_mge_arg_contigs.tsv")
    log:
        os.path.join(RESULTS_DIR, "logs/mge_arg_cooccurrence.{sid}.log")
    wildcard_constraints:
        sid="|".join(SAMPLES.index)
    message:
        "Step: {wildcards.sid} — MGE-ARG contig co-occurrence analysis"
    run:
        import os
        import pandas as pd

        os.makedirs(os.path.dirname(output[0]), exist_ok=True)

        # --- Load ISEScan summary -----------------------------------------------
        # ISEScan .sum file is tab-separated with a header; column 'seqID' holds
        # the contig/sequence identifier.
        try:
            ise_df = pd.read_csv(input.isescan_sum, sep="\t", comment="#")
            if "seqID" not in ise_df.columns:
                # Some ISEScan versions use 'Sequence' as the first column name
                ise_df.rename(columns={ise_df.columns[0]: "seqID"}, inplace=True)
            # Collect IS family labels per contig (column 'isFamily' or 'family')
            family_col = "isFamily" if "isFamily" in ise_df.columns else (
                "family" if "family" in ise_df.columns else None
            )
            if family_col:
                ise_by_contig = (
                    ise_df.groupby("seqID")[family_col]
                    .apply(lambda x: ";".join(sorted(set(x.dropna().astype(str)))))
                    .reset_index()
                    .rename(columns={family_col: "IS_families"})
                )
            else:
                ise_by_contig = (
                    ise_df[["seqID"]]
                    .drop_duplicates()
                    .assign(IS_families="unknown")
                )
            is_contigs = set(ise_by_contig["seqID"].astype(str))
        except Exception as e:
            print(f"WARNING: could not parse ISEScan summary ({input.isescan_sum}): {e}")
            ise_by_contig = pd.DataFrame(columns=["seqID", "IS_families"])
            is_contigs = set()

        # --- Load RGI output ----------------------------------------------------
        # RGI txt output is tab-separated; column 'Contig' holds the contig ID
        # (which may include ORF suffix, e.g. contig_1_1 — strip the trailing ORF
        # number so it matches the assembly contig identifier).
        try:
            rgi_df = pd.read_csv(input.rgi_txt, sep="\t")
            if "Contig" not in rgi_df.columns:
                # Older RGI versions use 'ORF_ID'; derive Contig from it
                rgi_df["Contig"] = rgi_df["ORF_ID"].str.rsplit("_", n=1).str[0]
            else:
                # Strip trailing ORF index appended by RGI (e.g. "_1")
                rgi_df["Contig"] = rgi_df["Contig"].str.rsplit("_", n=1).str[0]

            arg_col = "Best_Hit_ARO" if "Best_Hit_ARO" in rgi_df.columns else rgi_df.columns[1]
            arg_by_contig = (
                rgi_df.groupby("Contig")[arg_col]
                .apply(lambda x: ";".join(sorted(set(x.dropna().astype(str)))))
                .reset_index()
                .rename(columns={arg_col: "ARG_names"})
            )
            arg_contigs = set(arg_by_contig["Contig"].astype(str))
        except Exception as e:
            print(f"WARNING: could not parse RGI output ({input.rgi_txt}): {e}")
            arg_by_contig = pd.DataFrame(columns=["Contig", "ARG_names"])
            arg_contigs = set()

        # --- Build union of all relevant contigs --------------------------------
        all_contigs = is_contigs | arg_contigs
        result_rows = []
        for contig in sorted(all_contigs):
            has_is = contig in is_contigs
            has_arg = contig in arg_contigs
            is_fam_str = (
                ise_by_contig.loc[ise_by_contig["seqID"] == contig, "IS_families"].values[0]
                if has_is else ""
            )
            arg_name_str = (
                arg_by_contig.loc[arg_by_contig["Contig"] == contig, "ARG_names"].values[0]
                if has_arg else ""
            )
            result_rows.append({
                "contig_id": contig,
                "has_IS_element": has_is,
                "has_ARG": has_arg,
                "IS_families": is_fam_str,
                "ARG_names": arg_name_str
            })

        result_df = pd.DataFrame(result_rows, columns=[
            "contig_id", "has_IS_element", "has_ARG", "IS_families", "ARG_names"
        ])
        result_df.to_csv(output[0], sep="\t", index=False)
        print(
            f"Sample co-occurrence written: {len(result_rows)} contigs "
            f"({len(is_contigs & arg_contigs)} with both IS element and ARG)"
        )


############################################
# Rule 5: Aggregate per-sample co-occurrence tables into a combined summary
rule summarise_mge_arg:
    input:
        expand(
            os.path.join(RESULTS_DIR, "mge_cooccurrence/{sid}/{sid}_mge_arg_contigs.tsv"),
            sid=SAMPLES.index
        )
    output:
        os.path.join(RESULTS_DIR, "mge_cooccurrence/combined_mge_arg_summary.tsv")
    log:
        os.path.join(RESULTS_DIR, "logs/summarise_mge_arg.log")
    message:
        "Summary: Aggregate MGE-ARG co-occurrence across all samples"
    run:
        import os
        import pandas as pd

        frames = []
        sample_ids = list(SAMPLES.index)
        for sid, fpath in zip(sample_ids, input):
            try:
                df = pd.read_csv(fpath, sep="\t")
                df.insert(0, "sample", sid)
                frames.append(df)
            except Exception as e:
                print(f"WARNING: could not read {fpath}: {e}")

        if not frames:
            combined = pd.DataFrame(columns=[
                "sample", "contig_id", "has_IS_element", "has_ARG",
                "IS_families", "ARG_names",
                "pct_ARG_contigs_with_IS"
            ])
            combined.to_csv(output[0], sep="\t", index=False)
        else:
            combined = pd.concat(frames, ignore_index=True)

            # Per-sample statistic: % of ARG-containing contigs that also carry an IS element
            stats = []
            for sid, grp in combined.groupby("sample"):
                arg_contigs = grp[grp["has_ARG"] == True]
                n_arg = len(arg_contigs)
                n_arg_with_is = int((arg_contigs["has_IS_element"] == True).sum())
                pct = round(100.0 * n_arg_with_is / n_arg, 2) if n_arg > 0 else 0.0
                stats.append({"sample": sid,
                               "n_arg_contigs": n_arg,
                               "n_arg_contigs_with_IS": n_arg_with_is,
                               "pct_ARG_contigs_with_IS": pct})
            stats_df = pd.DataFrame(stats)

            # Merge stat columns back onto the per-contig table
            combined = combined.merge(
                stats_df[["sample", "pct_ARG_contigs_with_IS"]],
                on="sample", how="left"
            )
            combined.to_csv(output[0], sep="\t", index=False)
            print(f"Combined MGE-ARG summary written to {output[0]}")
            print(stats_df.to_string(index=False))


############################################
# Rule 6: Aggregate per-sample BacMet2 hits into a combined metal resistance table
rule summarise_metal_resistance:
    input:
        expand(
            os.path.join(RESULTS_DIR, "metal_resistance/{sid}/{sid}_bacmet.tsv"),
            sid=SAMPLES.index
        )
    output:
        os.path.join(RESULTS_DIR, "metal_resistance/combined_metal_resistance.tsv")
    log:
        os.path.join(RESULTS_DIR, "logs/summarise_metal_resistance.log")
    message:
        "Summary: Aggregate BacMet2 metal resistance hits across all samples"
    run:
        import os
        import re
        import pandas as pd

        # Heuristic mapping from gene-name prefix to metal/compound class.
        # Expand as needed; keys are lower-case prefixes of the subject sequence ID.
        METAL_MAP = {
            "mera": "mercury",  "merb": "mercury",  "merc": "mercury",
            "merd": "mercury",  "mere": "mercury",  "merp": "mercury",
            "copa": "copper",   "copb": "copper",   "copc": "copper",
            "cusr": "copper",   "cusa": "copper",   "cusb": "copper",
            "cusc": "copper",   "cusd": "copper",   "cuse": "copper",
            "arca": "arsenic",  "arsb": "arsenic",  "arsc": "arsenic",
            "arsr": "arsenic",  "acr1": "arsenic",  "acr2": "arsenic",
            "acr3": "arsenic",  "arsa": "arsenic",
            "cad":  "cadmium",  "cada": "cadmium",  "cadc": "cadmium",
            "znt":  "zinc",     "znta": "zinc",     "czca": "zinc/cobalt",
            "czcb": "zinc/cobalt/cadmium", "czcc": "zinc/cobalt/cadmium",
            "ncc":  "nickel/cobalt/cadmium",
            "chr":  "chromate", "chra": "chromate", "chrb": "chromate",
            "pbr":  "lead",     "pbra": "lead",     "pbrb": "lead",
            "silp": "silver",   "silb": "silver",   "sile": "silver",
            "tela": "tellurite","telb": "tellurite",
            "arse": "arsenite", "arsH": "arsenic",
            "feo":  "iron",     "feoa": "iron",     "feob": "iron",
        }

        def infer_metal_class(gene_name: str) -> str:
            """Return a metal/compound class from the gene name prefix."""
            name_lower = str(gene_name).lower()
            for prefix, metal in METAL_MAP.items():
                if name_lower.startswith(prefix):
                    return metal
            return "unknown"

        blast_cols = [
            "qseqid", "sseqid", "pident", "length", "mismatch",
            "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore"
        ]

        sample_ids = list(SAMPLES.index)
        frames = []
        for sid, fpath in zip(sample_ids, input):
            try:
                df = pd.read_csv(fpath, sep="\t", header=None, names=blast_cols)
                df.insert(0, "sample", sid)
                # Parse gene name from subject sequence ID (first field before '|' or space)
                df["gene_name"] = df["sseqid"].str.split(r"[|\s]").str[0]
                df["metal_class"] = df["gene_name"].apply(infer_metal_class)
                frames.append(df)
            except Exception as e:
                print(f"WARNING: could not read {fpath}: {e}")

        if not frames:
            combined = pd.DataFrame(columns=["sample"] + blast_cols + ["gene_name", "metal_class"])
        else:
            combined = pd.concat(frames, ignore_index=True)

        combined.to_csv(output[0], sep="\t", index=False)
        print(f"Combined metal resistance table written to {output[0]}")
        if frames:
            summary = (
                combined.groupby(["sample", "metal_class"])
                .size()
                .reset_index(name="n_hits")
            )
            print(summary.to_string(index=False))
