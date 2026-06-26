"""
Author: Susheel Bhanu BUSI
Affiliation: Molecular Ecology group, UKCEH
Date: [2024-06-04]
Run: snakemake -s workflow/rules/organellar_filter.smk --use-conda --cores 4 -rp
Latest modification:
Purpose: Filter assembled contigs against chloroplast and mitochondrial reference
         sequences to remove organellar contamination.
         Addresses Reviewer 2 comment #5: the pipeline previously filtered only
         human reads (GRCh38) but not chloroplast (Cyanobacteriota) or mitochondrial
         (Pseudomonadota) sequences, which could inflate contig-level phylum counts
         for those two phyla.

# ---- config.yaml additions required ----------------------------------------
# organellar:
#   chloroplast_url: "https://ftp.ncbi.nlm.nih.gov/refseq/release/plastid/plastid.1.1.genomic.fna.gz"
#   mitochondria_url: "https://ftp.ncbi.nlm.nih.gov/refseq/release/mitochondrion/mitochondrion.1.1.genomic.fna.gz"
#   threads: 16
#   min_k: 31
# -----------------------------------------------------------------------------
"""


############################################
# Aggregate rule — always defined first
rule organellar_filter:
    input:
        expand(
            os.path.join(RESULTS_DIR, "assembly_filtered/{sid}/{sid}_noOrganellar.fasta"),
            sid=SAMPLES.index
        ),
        expand(
            os.path.join(RESULTS_DIR, "assembly_filtered/{sid}/{sid}_organellar_stats.txt"),
            sid=SAMPLES.index
        )
    output:
        touch("status/organellar_filter.done")


############################################
# Download chloroplast and mitochondrial reference sequences from NCBI RefSeq
rule download_organellar_db:
    output:
        chloroplast=os.path.join(DB_DIR, "organellar/chloroplast.fna"),
        mitochondria=os.path.join(DB_DIR, "organellar/mitochondria.fna"),
        combined=os.path.join(DB_DIR, "organellar/organellar_refs.fna")
    params:
        chloroplast_url=config["organellar"]["chloroplast_url"],
        mitochondria_url=config["organellar"]["mitochondria_url"],
        outdir=os.path.join(DB_DIR, "organellar")
    log:
        os.path.join(RESULTS_DIR, "logs/organellar/download_organellar_db.log")
    message:
        "Downloading organellar reference sequences (chloroplast + mitochondria) from NCBI RefSeq"
    shell:
        "(date && "
        "mkdir -p {params.outdir} && "
        "wget -q -O {params.outdir}/chloroplast.fna.gz '{params.chloroplast_url}' && "
        "gunzip -f {params.outdir}/chloroplast.fna.gz && "
        "wget -q -O {params.outdir}/mitochondria.fna.gz '{params.mitochondria_url}' && "
        "gunzip -f {params.outdir}/mitochondria.fna.gz && "
        "cat {output.chloroplast} {output.mitochondria} > {output.combined} && "
        "date) &> >(tee {log})"


############################################
# Per-sample: filter assembled contigs against the combined organellar reference
rule filter_organellar_contigs:
    input:
        fasta=os.path.join(RESULTS_DIR, "assembly/{sid}/{sid}.fasta"),
        db=os.path.join(DB_DIR, "organellar/organellar_refs.fna")
    output:
        clean=os.path.join(RESULTS_DIR, "assembly_filtered/{sid}/{sid}_noOrganellar.fasta"),
        organellar=os.path.join(RESULTS_DIR, "assembly_filtered/{sid}/{sid}_organellar_hits.fasta"),
        stats=os.path.join(RESULTS_DIR, "assembly_filtered/{sid}/{sid}_organellar_stats.txt")
    conda:
        os.path.join(ENV_DIR, "bbmap.yaml")
    threads:
        config["organellar"]["threads"]
    params:
        min_k=config["organellar"]["min_k"],
        outdir=lambda wildcards, output: os.path.dirname(output.clean)
    log:
        os.path.join(RESULTS_DIR, "logs/organellar/filter_organellar_contigs.{sid}.log")
    wildcard_constraints:
        sid="|".join(SAMPLES.index)
    message:
        "Step: {wildcards.sid}"
    benchmark:
        os.path.join(RESULTS_DIR, "benchmarks/organellar_filter.{sid}.txt")
    shell:
        "(date && "
        "mkdir -p {params.outdir} && "
        "bbduk.sh "
        "in={input.fasta} "
        "out={output.clean} "
        "outm={output.organellar} "
        "ref={input.db} "
        "k={params.min_k} "
        "hdist=1 "
        "stats={output.stats} "
        "threads={threads} && "
        "date) &> >(tee {log})"
