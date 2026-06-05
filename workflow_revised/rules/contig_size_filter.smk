"""
Author: Susheel Bhanu BUSI
Affiliation: Molecular Ecology group, UKCEH
Date: [2024-06-04]
Run: snakemake -s workflow/rules/contig_size_filter.smk --use-conda --cores 4 -rp
Latest modification:
Purpose: Filter assembled contigs by minimum length before ARG detection and
         contig-level taxonomic assignment (addresses Reviewer 2 comment #4).
         Two thresholds are applied:
           - >= 2 kbp : standard cutoff prior to ARG/AMR analysis and contig taxonomy
           - >= 10 kbp: more stringent cutoff for sensitivity analysis of contig taxonomy
         Per-sample stats are produced alongside each filtered FASTA, and a
         combined summary TSV is generated across all samples and thresholds.

         Add the following block to config.yaml:
           seqkit:
             threads: 4
             min_length_amr: 2000
             min_length_taxonomy: 10000
"""


############################################
rule contig_size_filter:
    input:
        expand(os.path.join(RESULTS_DIR, "assembly_filtered/{sid}/{sid}_min2kb.fasta"),  sid=SAMPLES.index),
        expand(os.path.join(RESULTS_DIR, "assembly_filtered/{sid}/{sid}_min10kb.fasta"), sid=SAMPLES.index),
        os.path.join(RESULTS_DIR, "assembly_filtered/contig_size_stats.tsv")
    output:
        touch("status/contig_size_filter.done")


############################################
# Filter contigs >= 2 kbp (for ARG / AMR analysis)
rule filter_contigs_2kb:
    input:
        os.path.join(DATA_DIR, "assembly/{sid}.fasta")
    output:
        fasta=os.path.join(RESULTS_DIR, "assembly_filtered/{sid}/{sid}_min2kb.fasta"),
        stats=os.path.join(RESULTS_DIR, "assembly_filtered/{sid}/{sid}_min2kb.stats")
    conda:
        os.path.join(ENV_DIR, "seqkit.yaml")
    threads:
        config['seqkit']['threads']
    log:
        os.path.join(RESULTS_DIR, "logs/contig_size_filter/{sid}_min2kb.log")
    message:
        "Step: {wildcards.sid} -- filtering contigs >= 2 kbp (AMR threshold)"
    wildcard_constraints:
        sid="|".join(SAMPLES.index)
    benchmark:
        os.path.join(RESULTS_DIR, "benchmarks/filter_contigs_1kb.{sid}.txt")
    shell:
        "(date && "
        "mkdir -p $(dirname {output.fasta}) && "
        "seqkit seq -m {config[seqkit][min_length_amr]} --threads {threads} {input} > {output.fasta} && "
        "seqkit stats -a --threads {threads} {output.fasta} > {output.stats} && "
        "date) &> >(tee {log})"


############################################
# Filter contigs >= 10 kbp (for contig-level taxonomic assignment)
rule filter_contigs_10kb:
    input:
        os.path.join(DATA_DIR, "assembly/{sid}.fasta")
    output:
        fasta=os.path.join(RESULTS_DIR, "assembly_filtered/{sid}/{sid}_min10kb.fasta"),
        stats=os.path.join(RESULTS_DIR, "assembly_filtered/{sid}/{sid}_min10kb.stats")
    conda:
        os.path.join(ENV_DIR, "seqkit.yaml")
    threads:
        config['seqkit']['threads']
    log:
        os.path.join(RESULTS_DIR, "logs/contig_size_filter/{sid}_min10kb.log")
    message:
        "Step: {wildcards.sid} -- filtering contigs >= 10 kbp (taxonomy threshold)"
    wildcard_constraints:
        sid="|".join(SAMPLES.index)
    benchmark:
        os.path.join(RESULTS_DIR, "benchmarks/filter_contigs_10kb.{sid}.txt")
    shell:
        "(date && "
        "mkdir -p $(dirname {output.fasta}) && "
        "seqkit seq -m {config[seqkit][min_length_taxonomy]} --threads {threads} {input} > {output.fasta} && "
        "seqkit stats -a --threads {threads} {output.fasta} > {output.stats} && "
        "date) &> >(tee {log})"


############################################
# Combine per-sample stats into a single summary TSV
rule contig_filter_stats:
    input:
        fasta_1kb=expand(os.path.join(RESULTS_DIR, "assembly_filtered/{sid}/{sid}_min2kb.fasta"),  sid=SAMPLES.index),
        fasta_10kb=expand(os.path.join(RESULTS_DIR, "assembly_filtered/{sid}/{sid}_min10kb.fasta"), sid=SAMPLES.index)
    output:
        os.path.join(RESULTS_DIR, "assembly_filtered/contig_size_stats.tsv")
    conda:
        os.path.join(ENV_DIR, "seqkit.yaml")
    threads:
        config['seqkit']['threads']
    log:
        os.path.join(RESULTS_DIR, "logs/contig_size_filter/contig_size_stats.log")
    message:
        "Combining contig size filter stats across all samples and thresholds"
    shell:
        "(date && "
        "seqkit stats -a --threads {threads} {input.fasta_1kb} {input.fasta_10kb} > {output} && "
        "date) &> >(tee {log})"
