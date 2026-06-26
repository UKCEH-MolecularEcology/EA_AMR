"""
Author: Susheel Bhanu BUSI
Affiliation: Molecular Ecology group, UKCEH
Date: [2024-06-04]
Run: snakemake -s workflow/rules/binning.smk --use-conda --cores 4 -rp
Latest modification:
Purpose: MAG-based analysis of ARG hosts to complement the contig-based approach.
         Addresses Reviewer 1 comments #3 & #4 and Reviewer 2 comment #3.
         Covers: MetaBAT2 binning, CheckM2 quality, HQ bin filtering, dRep
         dereplication, and GTDB-Tk taxonomy assignment.
         ARG annotation on MAGs is handled in a separate bin_amr.smk.

Config additions required:
  metabat2:
    threads: 56   # already present
    min_contig: 2000
  checkm:
    threads: 64   # already present
    extension: fa
  drep:
    threads: 64   # already present
  gtdbtk:
    threads: 64   # already present
    path: "/hdd0/susbus/databases/gtdbtk/release214"   # already present
"""


############################################
# Aggregate rule — must come first
rule binning:
    input:
        expand(os.path.join(RESULTS_DIR, "checkm/{sid}/quality_report.tsv"), sid=SAMPLES.index),
        os.path.join(RESULTS_DIR, "gtdbtk/classify/gtdbtk.bac120.summary.tsv")
    output:
        touch("status/binning.done")


############################################
wildcard_constraints:
    sid="|".join(SAMPLES.index)


############################################
localrules: filter_hq_bins


############################################
# 1. Per-sample binning with MetaBAT2
rule metabat2_binning:
    input:
        fasta=os.path.join(RESULTS_DIR, "assembly/{sid}/{sid}.fasta"),
        depth=os.path.join(RESULTS_DIR, "coverage/{sid}/{sid}_depth.txt")
    output:
        bin_dir=directory(os.path.join(RESULTS_DIR, "bins/{sid}/")),
        done=touch(os.path.join(RESULTS_DIR, "bins/{sid}/.done"))
    conda:
        os.path.join(ENV_DIR, "metabat2.yaml")
    threads:
        config["metabat2"]["threads"]
    log:
        os.path.join(RESULTS_DIR, "logs/metabat2.{sid}.log")
    message:
        "Step: {wildcards.sid} — MetaBAT2 binning"
    benchmark:
        os.path.join(RESULTS_DIR, "benchmarks/metabat2.{sid}.txt")
    shell:
        "(date && "
        "mkdir -p {output.bin_dir} && "
        "metabat2 -i {input.fasta} -a {input.depth} "
        "-o {output.bin_dir}/{wildcards.sid} "
        "--minContig 2000 -t {threads} --unbinned && "
        "date) &> >(tee {log})"


############################################
# 2. Per-sample bin quality assessment with CheckM2
rule checkm2_quality:
    input:
        bin_dir=os.path.join(RESULTS_DIR, "bins/{sid}/")
    output:
        report=os.path.join(RESULTS_DIR, "checkm/{sid}/quality_report.tsv"),
        dir=directory(os.path.join(RESULTS_DIR, "checkm/{sid}/"))
    conda:
        os.path.join(ENV_DIR, "checkm.yaml")
    threads:
        config["checkm"]["threads"]
    log:
        os.path.join(RESULTS_DIR, "logs/checkm2.{sid}.log")
    message:
        "Step: {wildcards.sid} — CheckM2 bin quality assessment"
    benchmark:
        os.path.join(RESULTS_DIR, "benchmarks/checkm2.{sid}.txt")
    shell:
        "(date && "
        "checkm2 predict --threads {threads} "
        "--input {input.bin_dir}/*.fa "
        "--output-directory {output.dir} "
        "-x fa && "
        "date) &> >(tee {log})"


############################################
# 3. Filter HQ bins (completeness >= 50, contamination <= 10)
rule filter_hq_bins:
    input:
        report=os.path.join(RESULTS_DIR, "checkm/{sid}/quality_report.tsv"),
        bin_dir=os.path.join(RESULTS_DIR, "bins/{sid}/")
    output:
        hq_list=os.path.join(RESULTS_DIR, "checkm/{sid}/hq_bins.txt"),
        hq_dir=directory(os.path.join(RESULTS_DIR, "bins_hq/{sid}/"))
    log:
        os.path.join(RESULTS_DIR, "logs/filter_hq_bins.{sid}.log")
    message:
        "Step: {wildcards.sid} — Filtering HQ bins (completeness >= 50, contamination <= 10)"
    run:
        import shutil
        import csv

        os.makedirs(output.hq_dir, exist_ok=True)

        hq_paths = []
        with open(input.report, newline="") as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            for row in reader:
                try:
                    completeness = float(row["Completeness"])
                    contamination = float(row["Contamination"])
                except (KeyError, ValueError):
                    continue
                if completeness >= 50 and contamination <= 10:
                    bin_name = row["Name"]
                    # CheckM2 reports bins without extension; try .fa
                    src = os.path.join(input.bin_dir, bin_name + ".fa")
                    if not os.path.exists(src):
                        # fallback: bin name may already carry extension
                        src = os.path.join(input.bin_dir, bin_name)
                    if os.path.exists(src):
                        dest = os.path.join(output.hq_dir, os.path.basename(src))
                        shutil.copy2(src, dest)
                        hq_paths.append(dest)

        with open(output.hq_list, "w") as out_fh:
            for p in hq_paths:
                out_fh.write(p + "\n")

        with open(log[0], "w") as log_fh:
            log_fh.write(
                "Filtered {:d} HQ bins for sample {}\n".format(
                    len(hq_paths), wildcards.sid
                )
            )


############################################
# 4. Global dereplication with dRep across all samples
rule dereplicate_bins:
    input:
        bins=expand(
            os.path.join(RESULTS_DIR, "bins_hq/{sid}/"),
            sid=SAMPLES.index
        )
    output:
        dir=directory(os.path.join(RESULTS_DIR, "drep/dereplicated_genomes/"))
    conda:
        os.path.join(ENV_DIR, "drep.yaml")
    threads:
        config["drep"]["threads"]
    log:
        os.path.join(RESULTS_DIR, "logs/drep.log")
    message:
        "Step: global — dRep dereplication of HQ bins across all samples"
    benchmark:
        os.path.join(RESULTS_DIR, "benchmarks/drep.txt")
    shell:
        "(date && "
        "dRep dereplicate {output.dir} "
        "-g {input.bins}/*.fa "
        "-comp 50 -con 10 "
        "-p {threads} "
        "--S_algorithm fastANI && "
        "date) &> >(tee {log})"


############################################
# 5. Taxonomy assignment with GTDB-Tk
rule gtdbtk_classify:
    input:
        dir=os.path.join(RESULTS_DIR, "drep/dereplicated_genomes/")
    output:
        summary=os.path.join(RESULTS_DIR, "gtdbtk/classify/gtdbtk.bac120.summary.tsv"),
        dir=directory(os.path.join(RESULTS_DIR, "gtdbtk/"))
    conda:
        os.path.join(ENV_DIR, "gtdbtk.yaml")
    threads:
        config["gtdbtk"]["threads"]
    log:
        os.path.join(RESULTS_DIR, "logs/gtdbtk.log")
    params:
        db=config["gtdbtk"]["path"]
    message:
        "Step: global — GTDB-Tk taxonomy classification of dereplicated MAGs"
    benchmark:
        os.path.join(RESULTS_DIR, "benchmarks/gtdbtk.txt")
    shell:
        "(date && "
        "export GTDBTK_DATA_PATH={params.db} && "
        "gtdbtk classify_wf "
        "--genome_dir {input.dir} "
        "--out_dir {output.dir} "
        "--cpus {threads} "
        "--skip_ani_screen "
        "-x fa && "
        "date) &> >(tee {log})"
