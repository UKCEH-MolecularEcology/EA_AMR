"""
Author: Susheel Bhanu BUSI
Affiliation: Molecular Ecology group, UKCEH
Date: [2024-06-04]
Run: snakemake -s workflow_revised/Snakefile --configfile config/config.yaml --use-conda -rp
Latest modification:
Purpose: Run geNomad on assembled contigs to classify each contig as chromosome,
         plasmid, or virus (phage). Output feeds into contig_summary.smk for the
         integrated contig-level master table, and into mge_metal.smk for MGE-ARG
         and ARG-metal co-occurrence analysis.

# ---- config.yaml additions required ----------------------------------------
# genomad:
#   threads: 16
#   min_score: 0.7        # minimum classification score (0–1)
#   min_plasmid_marker_enrichment: 0.0   # default, keep as-is
#   db_url: "https://portal.nersc.gov/genomad/genomad_db_v1.7.tar.gz"
# -----------------------------------------------------------------------------
"""


############################################
rule genomad:
    input:
        expand(
            os.path.join(
                RESULTS_DIR,
                "genomad/{sid}/{sid}_noOrganellar_aggregated_classification/{sid}_noOrganellar_aggregated_classification.tsv"
            ),
            sid=SAMPLES.index
        )
    output:
        touch("status/genomad.done")


############################################
localrules: download_genomad_db


############################################
# Download geNomad database — only needed if db_path is not set in config.
# If config["genomad"]["db_path"] points to an existing database, this rule
# is never triggered (run_genomad uses the pre-existing path directly).
rule download_genomad_db:
    output:
        flag=os.path.join(DB_DIR, "genomad/genomad_db.done")
    params:
        db_dir=os.path.join(DB_DIR, "genomad")
    log:
        os.path.join(RESULTS_DIR, "logs/genomad/download_genomad_db.log")
    conda:
        os.path.join(ENV_DIR, "genomad.yaml")
    message:
        "Setup: downloading geNomad database"
    shell:
        "(date && "
        "mkdir -p {params.db_dir} && "
        "genomad download-database {params.db_dir} && "
        "touch {output.flag} && "
        "date) &> >(tee {log})"


############################################
# Run geNomad end-to-end per sample.
# db_path is resolved from config["genomad"]["db_path"] if set (use pre-existing
# database), otherwise falls back to the downloaded copy in DB_DIR/genomad/.
rule run_genomad:
    input:
        fasta=os.path.join(RESULTS_DIR, "assembly_filtered/{sid}/{sid}_noOrganellar.fasta")
    priority: -1
    output:
        classification=os.path.join(
            RESULTS_DIR,
            "genomad/{sid}/{sid}_noOrganellar_aggregated_classification/{sid}_noOrganellar_aggregated_classification.tsv"
        ),
        plasmid_summary=os.path.join(
            RESULTS_DIR,
            "genomad/{sid}/{sid}_noOrganellar_summary/{sid}_noOrganellar_plasmid_summary.tsv"
        ),
        virus_summary=os.path.join(
            RESULTS_DIR,
            "genomad/{sid}/{sid}_noOrganellar_summary/{sid}_noOrganellar_virus_summary.tsv"
        )
    params:
        db_dir=config.get("genomad", {}).get(
            "db_path", os.path.join(DB_DIR, "genomad/genomad_db")
        ),
        out_dir=os.path.join(RESULTS_DIR, "genomad/{sid}"),
        min_score=config.get("genomad", {}).get("min_score", 0.7),
        bin=config.get("genomad", {}).get("bin", "genomad")
    threads:
        config.get("genomad", {}).get("threads", 16)
    conda:
        os.path.join(ENV_DIR, "genomad.yaml")
    log:
        os.path.join(RESULTS_DIR, "logs/genomad/run_genomad.{sid}.log")
    wildcard_constraints:
        sid="|".join(SAMPLES.index)
    message:
        "geNomad: classifying contigs as chromosome/plasmid/virus for {wildcards.sid}"
    shell:
        "(date && "
        "{params.bin} end-to-end "
        "--min-score {params.min_score} "
        "--threads {threads} "
        "--cleanup "
        "{input.fasta} "
        "{params.out_dir} "
        "{params.db_dir} && "
        "date) &> >(tee {log})"


############################################
# Combine per-sample geNomad classifications into a single table
rule combine_genomad:
    input:
        expand(
            os.path.join(
                RESULTS_DIR,
                "genomad/{sid}/{sid}_noOrganellar_aggregated_classification/{sid}_noOrganellar_aggregated_classification.tsv"
            ),
            sid=SAMPLES.index
        )
    output:
        os.path.join(RESULTS_DIR, "genomad/combined_genomad_classification.tsv")
    log:
        os.path.join(RESULTS_DIR, "logs/genomad/combine_genomad.log")
    message:
        "geNomad: combining per-sample classifications"
    run:
        import os
        import pandas as pd

        dfs = []
        for fp in input:
            sid = fp.split("/genomad/")[1].split("/")[0]
            try:
                df = pd.read_csv(fp, sep="\t")
                df.insert(0, "sample", sid)
                dfs.append(df)
            except (pd.errors.EmptyDataError, FileNotFoundError):
                pass

        combined = pd.concat(dfs, ignore_index=True) if dfs else pd.DataFrame()
        os.makedirs(os.path.dirname(str(output[0])), exist_ok=True)
        combined.to_csv(str(output[0]), sep="\t", index=False)
