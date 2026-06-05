"""
Author: Susheel Bhanu BUSI
Affiliation: Molecular Ecology group, UKCEH
Date: [2024-06-04]
Run: snakemake -s workflow_revised/Snakefile --configfile config/config.yaml --use-conda -rp
Latest modification:
Purpose: (1) Run Bracken on pre-existing read-level Kraken2 reports (DATA_DIR/kraken2/)
             and combine outputs for use in taxonomy_validation.smk.
         (2) Re-run Kraken2 at higher confidence (default 0.5 vs original 0.2) and
             run Bracken on the new reports as a sensitivity analysis.
             Addresses Reviewer 2 comment #4 on Kraken2 accuracy: results consistent
             across confidence thresholds strengthen the taxonomic assignments.

         Standard chain  (existing reports, confidence=0.2):
           DATA_DIR/kraken2/{sid}_kraken.report → Bracken → combined_bracken.txt
         High-confidence chain (rerun from reads, confidence=0.5):
           DATA_DIR/reads/{sid}_filtered.R{1,2}.fq → Kraken2(0.5) → Bracken
             → combined_bracken_highconf.txt

         Add to config.yaml:
           kraken2:
             highconf_confidence: 0.5   # sensitivity analysis threshold
"""


############################################
# Confidence levels to run (stored as strings to use safely in paths).
# Map them to float values for the Kraken2 --confidence flag.
KRAKEN2_CONF_LEVELS = config.get("kraken2", {}).get("highconf_levels", ["05", "07"])
KRAKEN2_CONF_FLOAT  = {"05": 0.5, "07": 0.7, "02": 0.2}

# Build a mapping {sid: (R1_path, R2_path)} from the actual .fq.gz files spread
# across ea_biofilm_results, ea_biofilm_results_AT, ea_biofilm_results_AT_additional.
# Symlinks in DATA_DIR/reads/ point to .fq paths that do not exist — only .fq.gz does.
import glob as _glob

_READ_MAP = {}
for _fq in _glob.glob("/prj/DECODE/ea_biofilm_results*/preprocessed/reads/*/*_filtered.R1.fq.gz"):
    _sid = os.path.basename(_fq).replace("_filtered.R1.fq.gz", "")
    _r2  = _fq.replace("_filtered.R1.fq.gz", "_filtered.R2.fq.gz")
    if os.path.exists(_r2):
        _READ_MAP[_sid] = (_fq, _r2)

def _get_read(wildcards, read):
    if wildcards.sid not in _READ_MAP:
        raise ValueError(f"No reads found for sample {wildcards.sid} in ea_biofilm_results*/")
    return _READ_MAP[wildcards.sid][0 if read == "R1" else 1]


############################################
rule bracken_reads:
    input:
        os.path.join(RESULTS_DIR, "bracken/combined_bracken.txt")
        # High-confidence reruns (c=0.5, c=0.7) are defined below but require
        # reads at DATA_DIR/reads/ — enable by adding "bracken_highconf" to steps.
    output:
        touch("status/bracken_reads.done")


############################################
localrules: combine_bracken_reads, remove_uncultured_reads


############################################
# Run Bracken on the existing per-sample Kraken2 read-level reports
rule bracken_from_existing:
    input:
        report=os.path.join(DATA_DIR, "kraken2/{sid}_kraken.report")
    output:
        bracken=os.path.join(RESULTS_DIR, "bracken/{sid}.bracken"),
        report=os.path.join(RESULTS_DIR, "bracken/{sid}_bracken.report")
    conda:
        os.path.join(ENV_DIR, "bracken.yaml")
    threads:
        config['kraken2']['threads']
    params:
        db=config['kraken2']['db'],
        read=config['kraken2']['read'],
        level=config['kraken2']['level'],
        bracken=config['bracken']['bin']
    log:
        os.path.join(RESULTS_DIR, "logs/bracken_reads/bracken.{sid}.log")
    wildcard_constraints:
        sid="|".join(SAMPLES.index)
    message:
        "Bracken (from existing Kraken2 reports): {wildcards.sid}"
    shell:
        "(date && "
        "{params.bracken} "
        "-d {params.db} "
        "-i {input.report} "
        "-o {output.bracken} "
        "-w {output.report} "
        "-r {params.read} "
        "-l {params.level} && "
        "date) &> >(tee {log})"


############################################
# Strip uncultured / endosymbiont / Incertae Sedis entries before combining
rule remove_uncultured_reads:
    input:
        bracken=os.path.join(RESULTS_DIR, "bracken/{sid}.bracken")
    output:
        edited=os.path.join(RESULTS_DIR, "bracken/{sid}_edited.bracken")
    log:
        os.path.join(RESULTS_DIR, "logs/bracken_reads/edited_bracken_{sid}.log")
    wildcard_constraints:
        sid="|".join(SAMPLES.index)
    message:
        "Removing uncultured/endosymbiont entries from Bracken output: {wildcards.sid}"
    shell:
        "(date && "
        "grep -v 'uncultured' {input.bracken} "
        "| grep -v 'endosymbionts' "
        "| grep -v 'Incertae Sedis' "
        "> {output.edited} && "
        "date) &> >(tee {log})"


############################################
# Combine all per-sample Bracken outputs into one table
rule combine_bracken_reads:
    input:
        bracken=expand(
            os.path.join(RESULTS_DIR, "bracken/{sid}_edited.bracken"),
            sid=SAMPLES.index
        )
    output:
        out=os.path.join(RESULTS_DIR, "bracken/combined_bracken.txt")
    conda:
        os.path.join(ENV_DIR, "python2.yaml")
    params:
        combine=config['bracken']['combine']
    log:
        os.path.join(RESULTS_DIR, "logs/bracken_reads/bracken_combine.log")
    message:
        "Combining all Bracken read-level outputs"
    shell:
        "(date && "
        "python {params.combine} --files {input.bracken} -o {output.out} && "
        "date) &> >(tee {log})"


############################################
# ── HIGH-CONFIDENCE SENSITIVITY ANALYSIS ─────────────────────────────────────
# Runs Kraken2 at multiple confidence thresholds (default: 0.5 and 0.7) to
# validate that taxonomic profiles are robust to more stringent classification.
# Note: --confidence is a Kraken2 k-mer level flag; it cannot be applied
# retroactively in Bracken. Bracken's -t flag is a minimum read-count filter,
# not a classification confidence.
#
# Add to config.yaml to customise:
#   kraken2:
#     highconf_levels: ["05", "07"]   # string keys, no decimal point
############################################

rule kraken2_reads_highconf:
    input:
        r1=lambda wildcards: SAMPLES.loc[wildcards.sid, "sR1"],
        r2=lambda wildcards: SAMPLES.loc[wildcards.sid, "sR2"]
    output:
        report=os.path.join(RESULTS_DIR, "kraken2_c{conf}/{sid}_kraken.report"),
        summary=os.path.join(RESULTS_DIR, "kraken2_c{conf}/{sid}_kraken.out")
    conda:
        os.path.join(ENV_DIR, "kraken2.yaml")
    threads:
        config['kraken2']['threads']
    params:
        db=config['kraken2']['db'],
        confidence=lambda wildcards: KRAKEN2_CONF_FLOAT.get(wildcards.conf, 0.5)
    log:
        os.path.join(RESULTS_DIR, "logs/bracken_reads/kraken2_c{conf}.{sid}.log")
    wildcard_constraints:
        sid="|".join(SAMPLES.index),
        conf="|".join(KRAKEN2_CONF_LEVELS)
    message:
        "Kraken2 confidence={params.confidence}: {wildcards.sid}"
    shell:
        "(date && "
        "kraken2 "
        "--threads {threads} "
        "--db {params.db} "
        "--confidence {params.confidence} "
        "--paired "
        "--output {output.summary} "
        "--report {output.report} "
        "{input.r1} {input.r2} && "
        "date) &> >(tee {log})"


rule bracken_highconf:
    input:
        report=os.path.join(RESULTS_DIR, "kraken2_c{conf}/{sid}_kraken.report")
    output:
        bracken=os.path.join(RESULTS_DIR, "bracken_c{conf}/{sid}.bracken"),
        report=os.path.join(RESULTS_DIR, "bracken_c{conf}/{sid}_bracken.report")
    conda:
        os.path.join(ENV_DIR, "bracken.yaml")
    threads:
        config['kraken2']['threads']
    params:
        db=config['kraken2']['db'],
        read=config['kraken2']['read'],
        level=config['kraken2']['level'],
        bracken=config['bracken']['bin']
    log:
        os.path.join(RESULTS_DIR, "logs/bracken_reads/bracken_c{conf}.{sid}.log")
    wildcard_constraints:
        sid="|".join(SAMPLES.index),
        conf="|".join(KRAKEN2_CONF_LEVELS)
    message:
        "Bracken (confidence={wildcards.conf}): {wildcards.sid}"
    shell:
        "(date && "
        "{params.bracken} "
        "-d {params.db} "
        "-i {input.report} "
        "-o {output.bracken} "
        "-w {output.report} "
        "-r {params.read} "
        "-l {params.level} && "
        "date) &> >(tee {log})"


rule remove_uncultured_highconf:
    input:
        bracken=os.path.join(RESULTS_DIR, "bracken_c{conf}/{sid}.bracken")
    output:
        edited=os.path.join(RESULTS_DIR, "bracken_c{conf}/{sid}_edited.bracken")
    log:
        os.path.join(RESULTS_DIR, "logs/bracken_reads/edited_bracken_c{conf}_{sid}.log")
    wildcard_constraints:
        sid="|".join(SAMPLES.index),
        conf="|".join(KRAKEN2_CONF_LEVELS)
    message:
        "Removing uncultured entries (confidence={wildcards.conf}): {wildcards.sid}"
    shell:
        "(date && "
        "grep -v 'uncultured' {input.bracken} "
        "| grep -v 'endosymbionts' "
        "| grep -v 'Incertae Sedis' "
        "> {output.edited} && "
        "date) &> >(tee {log})"


rule combine_bracken_highconf:
    input:
        bracken=lambda wildcards: expand(
            os.path.join(RESULTS_DIR, "bracken_c{conf}/{sid}_edited.bracken"),
            sid=SAMPLES.index,
            conf=wildcards.conf
        )
    output:
        out=os.path.join(RESULTS_DIR, "bracken_c{conf}/combined_bracken_c{conf}.txt")
    conda:
        os.path.join(ENV_DIR, "python2.yaml")
    params:
        combine=config['bracken']['combine']
    log:
        os.path.join(RESULTS_DIR, "logs/bracken_reads/bracken_combine_c{conf}.log")
    wildcard_constraints:
        conf="|".join(KRAKEN2_CONF_LEVELS)
    message:
        "Combining Bracken outputs (confidence={wildcards.conf})"
    shell:
        "(date && "
        "python {params.combine} --files {input.bracken} -o {output.out} && "
        "date) &> >(tee {log})"
