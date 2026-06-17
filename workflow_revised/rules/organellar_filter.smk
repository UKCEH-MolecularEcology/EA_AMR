"""
Author: Susheel Bhanu BUSI
Affiliation: Molecular Ecology group, UKCEH
Date: [2024-06-05]
Run: snakemake -s workflow_revised/Snakefile --configfile config/config.yaml --use-conda -rp
Latest modification:
Purpose: Filter assembled contigs against plastid and mitochondrion reference
         sequences to remove organellar contamination.
         Addresses Reviewer 2 comment #5: the pipeline previously filtered only
         human reads (GRCh38) but not plastid (Cyanobacteriota chloroplasts) or
         mitochondrial (Pseudomonadota) sequences, which could inflate contig-level
         phylum counts for those two most-discussed phyla.

         Method: BWA mem (same approach as human read filtering in preprocessing.smk).
         Contigs that map to the combined plastid+mitochondrion reference are removed;
         unmapped contigs are retained. Faster than BBDuk for assembled contigs.

# ---- config.yaml additions required ----------------------------------------
# organellar:
#   plastid_url: "https://ftp.ncbi.nlm.nih.gov/refseq/release/plastid/plastid.1.1.genomic.fna.gz"
#   mitochondrion_url: "https://ftp.ncbi.nlm.nih.gov/refseq/release/mitochondrion/mitochondrion.1.1.genomic.fna.gz"
#   threads: 16
# -----------------------------------------------------------------------------
"""


############################################
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
localrules: download_organellar_db


############################################
# Download plastid and mitochondrion reference sequences from NCBI RefSeq
rule download_organellar_db:
    output:
        plastid=os.path.join(DB_DIR, "organellar/plastid.fna"),
        mitochondrion=os.path.join(DB_DIR, "organellar/mitochondrion.fna"),
        combined=os.path.join(DB_DIR, "organellar/organellar_refs.fna")
    params:
        plastid_url=config["organellar"]["plastid_url"],
        mitochondrion_url=config["organellar"]["mitochondrion_url"],
        outdir=os.path.join(DB_DIR, "organellar")
    log:
        os.path.join(RESULTS_DIR, "logs/organellar/download_organellar_db.log")
    message:
        "Downloading organellar reference sequences (plastid + mitochondrion) from NCBI RefSeq"
    shell:
        "(date && "
        "mkdir -p {params.outdir} && "
        "wget -q -O {params.outdir}/plastid.fna.gz '{params.plastid_url}' && "
        "gunzip -f {params.outdir}/plastid.fna.gz && "
        "wget -q -O {params.outdir}/mitochondrion.fna.gz '{params.mitochondrion_url}' && "
        "gunzip -f {params.outdir}/mitochondrion.fna.gz && "
        "cat {output.plastid} {output.mitochondrion} > {output.combined} && "
        "date) &> >(tee {log})"


############################################
# Index the combined organellar reference for BWA
rule bwa_index_organellar:
    input:
        os.path.join(DB_DIR, "organellar/organellar_refs.fna")
    output:
        expand(
            os.path.join(DB_DIR, "organellar/organellar_refs.fna.{ext}"),
            ext=BWA_IDX_EXT
        )
    conda:
        os.path.join(ENV_DIR, "bwa.yaml")
    threads:
        config["bwa"]["threads"]
    log:
        os.path.join(RESULTS_DIR, "logs/organellar/bwa_index_organellar.log")
    message:
        "Indexing organellar reference with BWA"
    shell:
        "(date && bwa index {input} && date) &> >(tee {log})"


############################################
# Per-sample: map contigs to organellar reference, keep unmapped
rule filter_organellar_contigs:
    input:
        fasta=os.path.join(DATA_DIR, "assembly/{sid}.fasta"),
        idx=expand(
            os.path.join(DB_DIR, "organellar/organellar_refs.fna.{ext}"),
            ext=BWA_IDX_EXT
        )
    output:
        clean=os.path.join(RESULTS_DIR, "assembly_filtered/{sid}/{sid}_noOrganellar.fasta"),
        organellar=os.path.join(RESULTS_DIR, "assembly_filtered/{sid}/{sid}_organellar_hits.fasta"),
        stats=os.path.join(RESULTS_DIR, "assembly_filtered/{sid}/{sid}_organellar_stats.txt")
    conda:
        os.path.join(ENV_DIR, "bwa.yaml")
    threads:
        config["organellar"]["threads"]
    params:
        idx_prefix=lambda wildcards, input: os.path.splitext(input.idx[0])[0]
    log:
        os.path.join(RESULTS_DIR, "logs/organellar/filter_organellar_contigs.{sid}.log")
    wildcard_constraints:
        sid="|".join(SAMPLES.index)
    message:
        "Organellar filter (BWA): {wildcards.sid}"
    shell:
        "(date && "
        "mkdir -p $(dirname {output.clean}) && "
        # Map contigs to combined organellar reference; store temp sorted BAM
        "tmpbam=$(dirname {output.clean})/{wildcards.sid}_organellar_tmp.bam && "
        "bwa mem -t {threads} {params.idx_prefix} {input.fasta} 2>>{log} | "
        "samtools sort -@ {threads} -o $tmpbam && "
        # Unmapped contigs (flag 4 set) → clean FASTA
        "samtools view -b -f 4 -@ {threads} $tmpbam | samtools fasta - > {output.clean} && "
        # Mapped contigs (flag 4 not set) → organellar FASTA for QC
        "samtools view -b -F 4 -@ {threads} $tmpbam | samtools fasta - > {output.organellar} && "
        # Flagstat for stats file
        "samtools flagstat $tmpbam > {output.stats} && "
        "rm -f $tmpbam && "
        "date) &> >(tee {log})"
