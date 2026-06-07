"""
Author: Susheel Bhanu BUSI
Affiliation: Molecular Ecology group, UKCEH
Date: [2026-06-04]
Run: snakemake -s workflow/rules/binning.smk --use-conda --cores 64 -rp
Latest modification:
Purpose: MAG generation from metagenomic assemblies to complement contig-based ARG
         analysis (addresses Reviewer 1 #3&#4 and Reviewer 2 #3).
         Pipeline: MetaBAT2 → CheckM2 quality filtering → dRep dereplication → GTDB-Tk taxonomy.
"""


############################################
rule binning:
    input:
        os.path.join(RESULTS_DIR, "gtdbtk/classify/gtdbtk.bac120.summary.tsv"),
        expand(os.path.join(RESULTS_DIR, "checkm/{sid}/quality_report.tsv"), sid=SAMPLES.index)
    output:
        touch("status/binning.done")


############################################
# Per-sample rules
rule metabat2_binning:
    input:
        fasta = os.path.join(RESULTS_DIR, "assembly_filtered/{sid}/{sid}_noOrganellar.fasta"),
        depth = os.path.join(RESULTS_DIR, "coverage/{sid}/{sid}_depth.txt")
    output:
        flag  = touch(os.path.join(RESULTS_DIR, "bins/{sid}/metabat2.done"))
    params:
        outdir      = os.path.join(RESULTS_DIR, "bins/{sid}/{sid}"),
        min_contig  = config.get("metabat2", {}).get("min_contig", 2000)
    threads:
        config["metabat2"]["threads"]
    conda:
        os.path.join(ENV_DIR, "metabat2.yaml")
    wildcard_constraints:
        sid = "|".join(SAMPLES.index)
    log:
        os.path.join(RESULTS_DIR, "logs/metabat2.{sid}.log")
    message:
        "MetaBAT2 binning: {wildcards.sid}"
    shell:
        "(date && "
        "mkdir -p $(dirname {params.outdir}) && "
        "metabat2 -i {input.fasta} -a {input.depth} -o {params.outdir} "
        "--minContig {params.min_contig} -t {threads} --unbinned && "
        "date) &> >(tee {log})"


############################################
rule checkm2_quality:
    input:
        flag = os.path.join(RESULTS_DIR, "bins/{sid}/metabat2.done")
    output:
        report = os.path.join(RESULTS_DIR, "checkm/{sid}/quality_report.tsv")
    params:
        bin_dir = os.path.join(RESULTS_DIR, "bins/{sid}"),
        out_dir = os.path.join(RESULTS_DIR, "checkm/{sid}")
    threads:
        config["checkm"]["threads"]
    conda:
        "/prj/DECODE/conda_envs/checkm2"
    wildcard_constraints:
        sid = "|".join(SAMPLES.index)
    log:
        os.path.join(RESULTS_DIR, "logs/checkm2.{sid}.log")
    message:
        "CheckM2 quality assessment: {wildcards.sid}"
    shell:
        "(date && "
        "checkm2 predict --threads {threads} --input {params.bin_dir} "
        "--output-directory {params.out_dir} -x fa --force && "
        "date) &> >(tee {log})"


############################################
localrules: collect_hq_bins

rule collect_hq_bins:
    input:
        reports = expand(os.path.join(RESULTS_DIR, "checkm/{sid}/quality_report.tsv"), sid=SAMPLES.index)
    output:
        hq_manifest = os.path.join(RESULTS_DIR, "bins_hq/hq_bins_manifest.tsv")
    log:
        os.path.join(RESULTS_DIR, "logs/collect_hq_bins.log")
    message:
        "Collecting high-quality bins (completeness >= 50, contamination <= 10)"
    run:
        import os
        import pandas as pd

        hq_dir = os.path.join(RESULTS_DIR, "bins_hq")
        os.makedirs(hq_dir, exist_ok=True)

        records = []
        for report_path in input.reports:
            # Derive sample ID from path: .../checkm/{sid}/quality_report.tsv
            sid = os.path.basename(os.path.dirname(report_path))
            df = pd.read_csv(report_path, sep="\t")
            # Normalise column names (CheckM2 uses 'Completeness'/'Contamination')
            df.columns = [c.strip() for c in df.columns]
            hq = df[
                (df["Completeness"] >= 50) &
                (df["Contamination"] <= 10)
            ]
            for _, row in hq.iterrows():
                bin_name = row["Name"]  # e.g. sample1.1, sample1.unbinned
                # Derive bin number from the bin name (last field after the last dot)
                bin_number = bin_name.rsplit(".", 1)[-1]
                bin_path = os.path.join(RESULTS_DIR, "bins", sid, "{}.{}.fa".format(sid, bin_number))
                link_target = os.path.join(hq_dir, "{}.{}.fa".format(sid, bin_number))
                if os.path.exists(bin_path) and not os.path.exists(link_target):
                    os.symlink(bin_path, link_target)
                records.append({
                    "sample":        sid,
                    "bin_id":        bin_name,
                    "bin_path":      bin_path,
                    "completeness":  row["Completeness"],
                    "contamination": row["Contamination"]
                })

        manifest_df = pd.DataFrame(records, columns=["sample", "bin_id", "bin_path", "completeness", "contamination"])
        manifest_df.to_csv(output.hq_manifest, sep="\t", index=False)


############################################
rule dereplicate_bins:
    input:
        hq_manifest = os.path.join(RESULTS_DIR, "bins_hq/hq_bins_manifest.tsv")
    output:
        flag    = touch(os.path.join(RESULTS_DIR, "drep/drep.done")),
        drep_dir = directory(os.path.join(RESULTS_DIR, "drep"))
    params:
        bins_dir = os.path.join(RESULTS_DIR, "bins_hq"),
        comp     = config["drep"]["comp"],
        cont     = config["drep"]["cont"]
    threads:
        config["drep"]["threads"]
    conda:
        "/prj/DECODE/conda_envs/viwrap/ViWrap-dRep"
    log:
        os.path.join(RESULTS_DIR, "logs/drep.log")
    message:
        "dRep dereplication of high-quality bins"
    shell:
        "(date && "
        "dRep dereplicate {output.drep_dir} -g {params.bins_dir}/*.fa "
        "-comp {params.comp} -con {params.cont} -p {threads} --S_algorithm fastANI && "
        "date) &> >(tee {log})"


############################################
rule gtdbtk_classify:
    input:
        drep_flag = os.path.join(RESULTS_DIR, "drep/drep.done")
    output:
        summary = os.path.join(RESULTS_DIR, "gtdbtk/classify/gtdbtk.bac120.summary.tsv")
    params:
        genome_dir = os.path.join(RESULTS_DIR, "drep/dereplicated_genomes"),
        out_dir    = os.path.join(RESULTS_DIR, "gtdbtk"),
        db_path    = config["gtdbtk"]["path"]
    threads:
        config["gtdbtk"]["threads"]
    conda:
        "/prj/DECODE/conda_envs/viwrap/ViWrap-GTDBTk"
    log:
        os.path.join(RESULTS_DIR, "logs/gtdbtk.log")
    message:
        "GTDB-Tk taxonomic classification of dereplicated MAGs"
    shell:
        "(date && "
        "export GTDBTK_DATA_PATH={params.db_path} && "
        "gtdbtk classify_wf --genome_dir {params.genome_dir} "
        "--out_dir {params.out_dir} --cpus {threads} --skip_ani_screen -x fa && "
        "date) &> >(tee {log})"
