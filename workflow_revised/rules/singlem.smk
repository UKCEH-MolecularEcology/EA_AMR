"""
Author: Susheel Bhanu BUSI
Affiliation: Molecular Ecology group, UKCEH
Date: [2024-02-12]
Run: snakemake -s workflow/rules/singlem.smk --use-conda --cores 64 -rp
Latest modification:
Purpose: To run singlem on raw reads (singlem step) and on organellar-filtered
         assembled contigs (singlem_contigs step).
         Contig mode uses -f/--genome-fasta-files with genome-tuned defaults:
           --min-orf-length 300 (vs 72 for reads)
           --min-taxon-coverage 0.1 (vs 0.35 for reads)
         Output feeds into taxonomy_validation.smk as SingleM_contigs.
"""


############################################
rule singlem:
    input:
        os.path.join(RESULTS_DIR, "singlem/combined_singlem_otu.csv"),
        os.path.join(RESULTS_DIR, "singlem/combined_singlem_relab.csv")
    output:
        touch("status/singlem.done")


############################################
localrules: setup_singlem_db


############################################
# Download singleM database
rule setup_singlem_db:
    output:
        db=directory(os.path.join(DB_DIR, "singlem")), 
        dummy=os.path.join(DB_DIR, "singlem/db.done")
    log:
        os.path.join(RESULTS_DIR, "logs/setup.singlem.db.log")
    conda:
        "singlem"
    message:
        "Setup: download singleM database"
    shell:
        "(date && mkdir -p {output.db} && "
        "singlem data --output-directory {output.db} && "
        "touch {output.dummy} && date) &> >(tee {log})"

rule run_singlem:
    input:
        in1=lambda wildcards: SAMPLES.loc[wildcards.sid, "sR1"], 
        in2=lambda wildcards: SAMPLES.loc[wildcards.sid, "sR2"],
        dummy=os.path.join(DB_DIR, "singlem/db.done")
    output:
        profile=os.path.join(RESULTS_DIR, "singlem/{sid}_singlem_profile.tsv"),
        table=os.path.join(RESULTS_DIR, "singlem/{sid}_singlem_otu.csv")
    log:
        os.path.join(RESULTS_DIR, "logs/singlem/{sid}.log")
    conda:
        "singlem"
    threads:
        config["singlem"]["threads"]
    params:
        db=os.path.join(DB_DIR, "singlem")
    message:
        "Running singlem on: {wildcards.sid}"
    shell:
        "(date && export SINGLEM_METAPACKAGE_PATH={params.db}/{config[singlem][db]} && "
        "singlem pipe -1 {input[0]} -2 {input[1]} -p {output.profile} --otu-table {output.table} --threads {threads} && "
        "date) &> >(tee {log})"

rule summarise_singlem:
    input:
        table=expand(os.path.join(RESULTS_DIR, "singlem/{sid}_singlem_otu.csv"), sid=SAMPLES.index),
        profile=expand(os.path.join(RESULTS_DIR, "singlem/{sid}_singlem_profile.tsv"), sid=SAMPLES.index)
    output:
        df_otu=os.path.join(RESULTS_DIR, "singlem/combined_singlem_otu.csv"),
        df_relab=os.path.join(RESULTS_DIR, "singlem/combined_singlem_relab.csv")
    log:
        os.path.join(RESULTS_DIR, "logs/single/combine.log")
    conda:
        "singlem"
    message:
        "Combined all singlem outputs"
    shell:
        "(date && singlem summarise --input-otu-tables {input.table} --output-otu-table {output.df_otu} && "
        "singlem summarise --input-otu-tables {input.table} --input-taxonomic-profiles {input.profile} --output-species-by-site-relative-abundance {output.df_relab} && "
        "date) &> >(tee {log})"


############################################
# ── CONTIG-BASED SINGLEM ─────────────────────────────────────────────────────
# Uses -f/--genome-fasta-files (genome mode) on organellar-filtered assemblies.
# genome mode defaults: --min-orf-length 300, --min-taxon-coverage 0.1
############################################

rule singlem_contigs:
    input:
        os.path.join(RESULTS_DIR, "singlem_contigs/combined_singlem_contigs_relab.csv")
    output:
        touch("status/singlem_contigs.done")


rule run_singlem_contigs:
    input:
        fasta=os.path.join(RESULTS_DIR, "assembly_filtered/{sid}/{sid}_noOrganellar.fasta"),
        dummy=os.path.join(DB_DIR, "singlem/db.done")
    output:
        profile=os.path.join(RESULTS_DIR, "singlem_contigs/{sid}_singlem_contigs_profile.tsv"),
        table=os.path.join(RESULTS_DIR, "singlem_contigs/{sid}_singlem_contigs_otu.csv")
    log:
        os.path.join(RESULTS_DIR, "logs/singlem_contigs/{sid}.log")
    conda:
        "singlem"
    threads:
        config["singlem"]["threads"]
    priority: 1
    params:
        db=os.path.join(DB_DIR, "singlem")
    wildcard_constraints:
        sid="|".join(SAMPLES.index)
    message:
        "SingleM (contigs): {wildcards.sid}"
    shell:
        "(date && export SINGLEM_METAPACKAGE_PATH={params.db}/{config[singlem][db]} && "
        "singlem pipe -f {input.fasta} -p {output.profile} --otu-table {output.table} "
        "--output-extras --threads {threads} && "
        "date) &> >(tee {log})"


rule summarise_singlem_contigs:
    input:
        table=expand(os.path.join(RESULTS_DIR, "singlem_contigs/{sid}_singlem_contigs_otu.csv"), sid=SAMPLES.index),
        profile=expand(os.path.join(RESULTS_DIR, "singlem_contigs/{sid}_singlem_contigs_profile.tsv"), sid=SAMPLES.index)
    output:
        df_otu=os.path.join(RESULTS_DIR, "singlem_contigs/combined_singlem_contigs_otu.csv"),
        df_relab=os.path.join(RESULTS_DIR, "singlem_contigs/combined_singlem_contigs_relab.csv")
    log:
        os.path.join(RESULTS_DIR, "logs/singlem_contigs/combine.log")
    conda:
        "singlem"
    message:
        "Combining all SingleM contig outputs"
    shell:
        "(date && singlem summarise --input-otu-tables {input.table} --output-otu-table {output.df_otu} && "
        "singlem summarise --input-otu-tables {input.table} --input-taxonomic-profiles {input.profile} --output-species-by-site-relative-abundance {output.df_relab} && "
        "date) &> >(tee {log})"
