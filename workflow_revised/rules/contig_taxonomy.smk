"""
Author: Susheel Bhanu BUSI
Affiliation: Molecular Ecology group, UKCEH
Date: [2023-09-25]
Run: snakemake -s workflow/rules/contig_taxonomy.smk --use-conda --cores 4 -rp
Latest modification:
Purpose: To run Kraken2+BRACKEN on metagenome assemblies, i.e. contigs
"""


############################################
rule contig_taxonomy:
    input:
        # Unfiltered contigs
        os.path.join(RESULTS_DIR, "mpa_report/contig/combined_output.tsv"),
        os.path.join(RESULTS_DIR, "bracken/contig/combined_bracken.txt"),
        # >= 2 kb contigs
        os.path.join(RESULTS_DIR, "mpa_report/contig_2kb/combined_output.tsv"),
        os.path.join(RESULTS_DIR, "bracken/contig_2kb/combined_bracken.txt"),
        # >= 10 kb contigs (sensitivity analysis)
        os.path.join(RESULTS_DIR, "mpa_report/contig_10kb/combined_output.tsv"),
        os.path.join(RESULTS_DIR, "bracken/contig_10kb/combined_bracken.txt")
    output:
        touch("status/contig_taxonomy.done")


############################################
localrules: contig_combine_bracken, contig_combine_mpa, contig_combine_bracken_2kb, contig_combine_mpa_2kb, contig_combine_bracken_10kb, contig_combine_mpa_10kb


############################################
# Taxonomic classification using KRAKEN2
rule contig_kraken2:
    input:
        os.path.join(RESULTS_DIR, "assembly_filtered/{sid}/{sid}_noOrganellar.fasta"),
    output:
        report=os.path.join(RESULTS_DIR, "kraken2/contig/{sid}_kraken.report"),
        summary=os.path.join(RESULTS_DIR, "kraken2/contig/{sid}_kraken.out")
    conda:
        os.path.join(ENV_DIR, "kraken2.yaml")
    threads:
        config['kraken2']['threads']
    params:
        db=config['kraken2']['db'],
        confidence=config['kraken2']['contig_confidence']
    log:
        os.path.join(RESULTS_DIR, "logs/contig/kraken2.{sid}.log")
    wildcard_constraints:
        sid="|".join(SAMPLES.index)
    message:
        "Running kraken2 on {wildcards.sid}"
    shell:
        "(date && kraken2 --threads {threads} --db {params.db} --confidence {params.confidence} --output {output.summary} --report {output.report} {input} && date) &> >(tee {log})"

# Running Struo2 database
use rule contig_kraken2 as contig_struo2_kraken2 with:
    output:
        report=os.path.join(RESULTS_DIR, "kraken2/contig/struo2_{sid}_kraken.report"),
        summary=os.path.join(RESULTS_DIR, "kraken2/contig/struo2_{sid}_kraken.out")
    threads:
        config['struo2_kraken2']['threads']
    params:
        db=config['struo2_kraken2']['db'],
        confidence=config['kraken2']['confidence']
    log:
        os.path.join(RESULTS_DIR, "logs/contig/struo2_kraken2.{sid}.log")
    wildcard_constraints:
        sid="|".join(SAMPLES.index)
    message:
        "Running {wildcards.sid} with the Struo2_Kraken2 db"

# Running KRAKEN2+BRACKEN as suggested by KRAKEN2 website
rule contig_bracken:
    input:
        report=rules.contig_kraken2.output.report
    output:
        bracken=os.path.join(RESULTS_DIR, "bracken/contig/{sid}.bracken"),
        report=os.path.join(RESULTS_DIR, "bracken/contig/{sid}_bracken.report")
    threads:
        config['kraken2']['threads']
    conda:
        os.path.join(ENV_DIR, "bracken.yaml")
    params:
        db=config['kraken2']['db'],
        read=config['kraken2']['read'],
        level=config['kraken2']['contig_level'],
        bracken=config['bracken']['bin'],
        header="name\ttaxonomy_id\ttaxonomy_lvl\tkraken_assigned_reads\tadded_reads\tnew_est_reads\tfraction_total_reads"
    log:
        os.path.join(RESULTS_DIR, "logs/contig/bracken.{sid}.log")
    wildcard_constraints:
        sid="|".join(SAMPLES.index)
    message:
        "Running kraken & bracken for {wildcards.sid}"
    shell:
        "(date && "
        "{params.bracken} -d {params.db} -i {input.report} -o {output.bracken} -w {output.report} -r {params.read} -l {params.level} || "
        "(echo 'WARNING: Bracken found no reads at {params.level} level for {wildcards.sid} — writing empty output' && "
        " printf '{params.header}\\n' > {output.bracken} && "
        " cp {input.report} {output.report}) && "
        "date) &> >(tee {log})"

rule contig_remove_uncultured:
    input:
        bracken=os.path.join(RESULTS_DIR, "bracken/contig/{sid}.bracken")
    output:
        edited=os.path.join(RESULTS_DIR, "bracken/contig/{sid}_edited.bracken")
    log:
        os.path.join(RESULTS_DIR, "logs/contig/edited_bracken_{sid}")
    wildcard_constraints:
        sid="|".join(SAMPLES.index)
    message:
        "Removing 'uncultured' taxa from bracken output from {wildcards.sid} due to combining issues"
    shell:
        "(date && grep -v 'uncultured' {input.bracken} | grep -v 'endosymbionts' | grep -v 'Incertae Sedis' > {output.edited} && date) &> >(tee {log})"

rule contig_combine_bracken:
    input:
        bracken=expand(os.path.join(RESULTS_DIR, "bracken/contig/{sid}_edited.bracken"), sid=SAMPLES.index)
    output:
        out=os.path.join(RESULTS_DIR, "bracken/contig/combined_bracken.txt")
    conda:
        os.path.join(ENV_DIR, "python2.yaml")
    params:
        combine=config['bracken']['combine']
    log:
        os.path.join(RESULTS_DIR, "logs/contig/bracken_combine.log")
    message:
        "Combining all the output from BRACKEN"
    shell:
        "(date && python {params.combine} --files {input.bracken} -o {output.out} && date)  &> >(tee {log})"


#########################
### MPA-style report ###
rule contig_mpa_report:
    input:
        report=os.path.join(RESULTS_DIR, "bracken/contig/{sid}_bracken.report")
    output:
        mpa=os.path.join(RESULTS_DIR, "mpa_report/contig/{sid}_mpa.tsv")
    conda:
        os.path.join(ENV_DIR, "bracken_new.yaml")
    log:
        os.path.join(RESULTS_DIR, "logs/contig/mpa_{sid}.log")
    wildcard_constraints:
        sid="|".join(SAMPLES.index)
    message:
        "Creating mpa-style report for {wildcards.sid}"
    shell:
        "(date && kreport2mpa.py -r {input.report} -o {output.mpa} && date)  &> >(tee {log})"

rule contig_combine_mpa:
    input:
        mpa=expand(os.path.join(RESULTS_DIR, "mpa_report/contig/{sid}_mpa.tsv"), sid=SAMPLES.index)
    output:
        combined=os.path.join(RESULTS_DIR, "mpa_report/contig/combined_output.tsv")
    conda:
        os.path.join(ENV_DIR, "krakentools.yaml")
    params:
        combine=os.path.join(SRC_DIR, "combine_mpa_modified.py")
    log:
        os.path.join(RESULTS_DIR, "logs/contig/mpa_combine.log")
    message:
        "Creating a combined mpa-style report"
    shell:
        "(date && {params.combine} -i {input.mpa} -d $(dirname {output.combined}) && date)  &> >(tee {log})"


#########################
### Rules for size-filtered contigs (>=2 kb) ###

rule contig_kraken2_2kb:
    input:
        os.path.join(RESULTS_DIR, "assembly_filtered/{sid}/{sid}_min2kb.fasta"),
    output:
        report=os.path.join(RESULTS_DIR, "kraken2/contig_2kb/{sid}_kraken.report"),
        summary=os.path.join(RESULTS_DIR, "kraken2/contig_2kb/{sid}_kraken.out")
    conda:
        os.path.join(ENV_DIR, "kraken2.yaml")
    threads:
        config['kraken2']['threads']
    params:
        db=config['kraken2']['db'],
        confidence=config['kraken2']['contig_confidence']
    log:
        os.path.join(RESULTS_DIR, "logs/contig_2kb/kraken2.{sid}.log")
    wildcard_constraints:
        sid="|".join(SAMPLES.index)
    message:
        "Running kraken2 on size-filtered (>=2 kb) contigs for {wildcards.sid}"
    shell:
        "(date && kraken2 --threads {threads} --db {params.db} --confidence {params.confidence} --output {output.summary} --report {output.report} {input} && date) &> >(tee {log})"

rule contig_bracken_2kb:
    input:
        report=rules.contig_kraken2_2kb.output.report
    output:
        bracken=os.path.join(RESULTS_DIR, "bracken/contig_2kb/{sid}.bracken"),
        report=os.path.join(RESULTS_DIR, "bracken/contig_2kb/{sid}_bracken.report")
    threads:
        config['kraken2']['threads']
    conda:
        os.path.join(ENV_DIR, "bracken.yaml")
    params:
        db=config['kraken2']['db'],
        read=config['kraken2']['read'],
        level=config['kraken2']['contig_level'],
        bracken=config['bracken']['bin'],
        header="name\ttaxonomy_id\ttaxonomy_lvl\tkraken_assigned_reads\tadded_reads\tnew_est_reads\tfraction_total_reads"
    log:
        os.path.join(RESULTS_DIR, "logs/contig_2kb/bracken.{sid}.log")
    wildcard_constraints:
        sid="|".join(SAMPLES.index)
    message:
        "Running kraken & bracken on size-filtered (>=2 kb) contigs for {wildcards.sid}"
    shell:
        "(date && "
        "{params.bracken} -d {params.db} -i {input.report} -o {output.bracken} -w {output.report} -r {params.read} -l {params.level} || "
        "(echo 'WARNING: Bracken found no reads at {params.level} level for {wildcards.sid} (2kb) — writing empty output' && "
        " printf '{params.header}\\n' > {output.bracken} && "
        " cp {input.report} {output.report}) && "
        "date) &> >(tee {log})"

rule contig_remove_uncultured_2kb:
    input:
        bracken=os.path.join(RESULTS_DIR, "bracken/contig_2kb/{sid}.bracken")
    output:
        edited=os.path.join(RESULTS_DIR, "bracken/contig_2kb/{sid}_edited.bracken")
    log:
        os.path.join(RESULTS_DIR, "logs/contig_2kb/edited_bracken_{sid}")
    wildcard_constraints:
        sid="|".join(SAMPLES.index)
    message:
        "Removing 'uncultured' taxa from bracken output from {wildcards.sid} (2kb contigs)"
    shell:
        "(date && grep -v 'uncultured' {input.bracken} | grep -v 'endosymbionts' | grep -v 'Incertae Sedis' > {output.edited} && date) &> >(tee {log})"

rule contig_combine_bracken_2kb:
    input:
        bracken=expand(os.path.join(RESULTS_DIR, "bracken/contig_2kb/{sid}_edited.bracken"), sid=SAMPLES.index)
    output:
        out=os.path.join(RESULTS_DIR, "bracken/contig_2kb/combined_bracken.txt")
    conda:
        os.path.join(ENV_DIR, "python2.yaml")
    params:
        combine=config['bracken']['combine']
    log:
        os.path.join(RESULTS_DIR, "logs/contig_2kb/bracken_combine.log")
    message:
        "Combining all BRACKEN output for size-filtered (>=2 kb) contigs"
    shell:
        "(date && python {params.combine} --files {input.bracken} -o {output.out} && date)  &> >(tee {log})"

rule contig_mpa_report_2kb:
    input:
        report=os.path.join(RESULTS_DIR, "bracken/contig_2kb/{sid}_bracken.report")
    output:
        mpa=os.path.join(RESULTS_DIR, "mpa_report/contig_2kb/{sid}_mpa.tsv")
    conda:
        os.path.join(ENV_DIR, "bracken_new.yaml")
    log:
        os.path.join(RESULTS_DIR, "logs/contig_2kb/mpa_{sid}.log")
    wildcard_constraints:
        sid="|".join(SAMPLES.index)
    message:
        "Creating mpa-style report for size-filtered (>=2 kb) contigs for {wildcards.sid}"
    shell:
        "(date && kreport2mpa.py -r {input.report} -o {output.mpa} && date)  &> >(tee {log})"

rule contig_combine_mpa_2kb:
    input:
        mpa=expand(os.path.join(RESULTS_DIR, "mpa_report/contig_2kb/{sid}_mpa.tsv"), sid=SAMPLES.index)
    output:
        combined=os.path.join(RESULTS_DIR, "mpa_report/contig_2kb/combined_output.tsv")
    conda:
        os.path.join(ENV_DIR, "krakentools.yaml")
    params:
        combine=os.path.join(SRC_DIR, "combine_mpa_modified.py")
    log:
        os.path.join(RESULTS_DIR, "logs/contig_2kb/mpa_combine.log")
    message:
        "Creating a combined mpa-style report for size-filtered (>=2 kb) contigs"
    shell:
        "(date && {params.combine} -i {input.mpa} -d $(dirname {output.combined}) && date)  &> >(tee {log})"


#########################
### Rules for size-filtered contigs (>=10 kb) ###

# Taxonomic classification using KRAKEN2 on size-filtered (>=10 kb) contigs
rule contig_kraken2_10kb:
    input:
        os.path.join(RESULTS_DIR, "assembly_filtered/{sid}/{sid}_min10kb.fasta"),
    output:
        report=os.path.join(RESULTS_DIR, "kraken2/contig_10kb/{sid}_kraken.report"),
        summary=os.path.join(RESULTS_DIR, "kraken2/contig_10kb/{sid}_kraken.out")
    conda:
        os.path.join(ENV_DIR, "kraken2.yaml")
    threads:
        config['kraken2']['threads']
    params:
        db=config['kraken2']['db'],
        confidence=config['kraken2']['contig_confidence']
    log:
        os.path.join(RESULTS_DIR, "logs/contig_10kb/kraken2.{sid}.log")
    wildcard_constraints:
        sid="|".join(SAMPLES.index)
    message:
        "Running kraken2 on size-filtered (>=10 kb) contigs for {wildcards.sid}"
    shell:
        "(date && kraken2 --threads {threads} --db {params.db} --confidence {params.confidence} --output {output.summary} --report {output.report} {input} && date) &> >(tee {log})"

rule contig_bracken_10kb:
    input:
        report=rules.contig_kraken2_10kb.output.report
    output:
        bracken=os.path.join(RESULTS_DIR, "bracken/contig_10kb/{sid}.bracken"),
        report=os.path.join(RESULTS_DIR, "bracken/contig_10kb/{sid}_bracken.report")
    threads:
        config['kraken2']['threads']
    conda:
        os.path.join(ENV_DIR, "bracken.yaml")
    params:
        db=config['kraken2']['db'],
        read=config['kraken2']['read'],
        level=config['kraken2']['contig_level'],
        bracken=config['bracken']['bin'],
        header="name\ttaxonomy_id\ttaxonomy_lvl\tkraken_assigned_reads\tadded_reads\tnew_est_reads\tfraction_total_reads"
    log:
        os.path.join(RESULTS_DIR, "logs/contig_10kb/bracken.{sid}.log")
    wildcard_constraints:
        sid="|".join(SAMPLES.index)
    message:
        "Running kraken & bracken on size-filtered (>=10 kb) contigs for {wildcards.sid}"
    shell:
        "(date && "
        "{params.bracken} -d {params.db} -i {input.report} -o {output.bracken} -w {output.report} -r {params.read} -l {params.level} || "
        "(echo 'WARNING: Bracken found no reads at {params.level} level for {wildcards.sid} (10kb) — writing empty output' && "
        " printf '{params.header}\\n' > {output.bracken} && "
        " cp {input.report} {output.report}) && "
        "date) &> >(tee {log})"

rule contig_remove_uncultured_10kb:
    input:
        bracken=os.path.join(RESULTS_DIR, "bracken/contig_10kb/{sid}.bracken")
    output:
        edited=os.path.join(RESULTS_DIR, "bracken/contig_10kb/{sid}_edited.bracken")
    log:
        os.path.join(RESULTS_DIR, "logs/contig_10kb/edited_bracken_{sid}")
    wildcard_constraints:
        sid="|".join(SAMPLES.index)
    message:
        "Removing 'uncultured' taxa from bracken output from {wildcards.sid} (10kb contigs) due to combining issues"
    shell:
        "(date && grep -v 'uncultured' {input.bracken} | grep -v 'endosymbionts' | grep -v 'Incertae Sedis' > {output.edited} && date) &> >(tee {log})"

rule contig_combine_bracken_10kb:
    input:
        bracken=expand(os.path.join(RESULTS_DIR, "bracken/contig_10kb/{sid}_edited.bracken"), sid=SAMPLES.index)
    output:
        out=os.path.join(RESULTS_DIR, "bracken/contig_10kb/combined_bracken.txt")
    conda:
        os.path.join(ENV_DIR, "python2.yaml")
    params:
        combine=config['bracken']['combine']
    log:
        os.path.join(RESULTS_DIR, "logs/contig_10kb/bracken_combine.log")
    message:
        "Combining all BRACKEN output for size-filtered (>=10 kb) contigs"
    shell:
        "(date && python {params.combine} --files {input.bracken} -o {output.out} && date)  &> >(tee {log})"

rule contig_mpa_report_10kb:
    input:
        report=os.path.join(RESULTS_DIR, "bracken/contig_10kb/{sid}_bracken.report")
    output:
        mpa=os.path.join(RESULTS_DIR, "mpa_report/contig_10kb/{sid}_mpa.tsv")
    conda:
        os.path.join(ENV_DIR, "bracken_new.yaml")
    log:
        os.path.join(RESULTS_DIR, "logs/contig_10kb/mpa_{sid}.log")
    wildcard_constraints:
        sid="|".join(SAMPLES.index)
    message:
        "Creating mpa-style report for size-filtered (>=10 kb) contigs for {wildcards.sid}"
    shell:
        "(date && kreport2mpa.py -r {input.report} -o {output.mpa} && date)  &> >(tee {log})"

rule contig_combine_mpa_10kb:
    input:
        mpa=expand(os.path.join(RESULTS_DIR, "mpa_report/contig_10kb/{sid}_mpa.tsv"), sid=SAMPLES.index)
    output:
        combined=os.path.join(RESULTS_DIR, "mpa_report/contig_10kb/combined_output.tsv")
    conda:
        os.path.join(ENV_DIR, "krakentools.yaml")
    params:
        combine=os.path.join(SRC_DIR, "combine_mpa_modified.py")
    log:
        os.path.join(RESULTS_DIR, "logs/contig_10kb/mpa_combine.log")
    message:
        "Creating a combined mpa-style report for size-filtered (>=10 kb) contigs"
    shell:
        "(date && {params.combine} -i {input.mpa} -d $(dirname {output.combined}) && date)  &> >(tee {log})"
