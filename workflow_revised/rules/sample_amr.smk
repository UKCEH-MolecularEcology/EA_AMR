"""
Author: Susheel Bhanu BUSI
Affiliation: Molecular Ecology group, UKCEH
Date: [2024-10-16]
Run: snakemake -s workflow/rules/amr.smk --use-conda --cores 64 -rp
Latest modification:
Purpose: To run AMR on individual assembly
"""


############################################
rule amr:
    input:
        expand(os.path.join(RESULTS_DIR, "amr/{sid}/{sid}_rgi.txt"), sid=SAMPLES.index),
        expand(os.path.join(RESULTS_DIR, "amr/coverage/{sid}/{sid}_rgi_coverage.txt"), sid=SAMPLES.index),
        expand(os.path.join(RESULTS_DIR, "fetchMG/{sid}/{sid}.faa.fetchMGs.scores"), sid=SAMPLES.index),
        os.path.join(RESULTS_DIR, "amr/normalized/arg_normalized_abundance.tsv"),
        os.path.join(RESULTS_DIR, "amr/normalized/genome_equivalents.tsv")
    output:
        touch("status/sample_amr.done")


############################################
localrules: download_rgi_db, setup_rgi_db


############################################
# Download RGI data
rule download_rgi_db:
    output:
        archive=temp(os.path.join(DB_DIR, "rgi/card-data.tar.bz2")),
        json=os.path.join(DB_DIR, "rgi/card.json")
    log:
        os.path.join(RESULTS_DIR, "logs/setup.rgi.db.log")
    params:
        db_url=config["rgi"]["db_url"]
    message:
        "Setup: download RGI data"
    shell:
        "(date && "
        "wget -O {output.archive} {params.db_url} --no-check-certificate && "
        "tar -C $(dirname {output.archive}) -xvf {output.archive} && "
        "date) &> >(tee {log})"

# Setup RGI: load required DB
# NOTE: to make sure that the same DB is used for all targets
rule setup_rgi_db:
    input:
        os.path.join(DB_DIR, "rgi/card.json")
    output:
        "status/rgi_setup.done"
    log:
        os.path.join(RESULTS_DIR, "logs/rgi.setup.log")
    conda:
        os.path.join(ENV_DIR, "rgi.yaml")
    message:
        "Setup: load RGI DB"
    shell:
        "(date && "
        "rgi clean --local && "
        "rgi load --card_json {input} --local && "
        "rgi database --version --local && "
        "touch {output} && date) &> >(tee {log})"

# Run RGI: Assembly (DNA)
rule sample_annotation_rgi:
    input:
        fna=os.path.join(DATA_DIR, "assembly/{sid}.fasta"),
        db=os.path.join(DB_DIR, "rgi/card.json"),
        setup="status/rgi_setup.done" # NOTE: to make sure that the same DB is used for all targets
    output:
        txt=os.path.join(RESULTS_DIR, "amr/{sid}/{sid}_rgi.txt")
    log:
        os.path.join(RESULTS_DIR, "logs/rgi.{sid}.log")
    threads:
        config["prodigal"]["threads"]
    params:
        alignment_tool="DIAMOND"
    conda:
        os.path.join(ENV_DIR, "rgi.yaml")
    message:
        "Annotation: RGI: {wildcards.sid}"
    shell:
        "(date && "
        "rgi database --version --local && "
        # NOTE: https://github.com/arpcard/rgi/issues/93: KeyError: 'snp' --> re-run
        "rgi main --input_sequence {input.fna} --output_file $(dirname {output})/$(basename -s '.txt' {output}) --local -a {params.alignment_tool} --clean --low_quality -n {threads} || "
        "rgi main --input_sequence {input.fna} --output_file $(dirname {output})/$(basename -s '.txt' {output}) --local -a {params.alignment_tool} --clean --low_quality -n {threads} && "
        "date) &> >(tee {log})"


# Merging RGI output with `gene_coverage.txt` 
#rule merge_rgi_coverage:
#    input:
#        rgi_file=os.path.join(RESULTS_DIR, "amr/{sid}/{sid}_rgi.txt"),
#        coverage_file=os.path.join(DATA_DIR, "coverages/{sid}_gene_coverage.txt")
#    output:
#        os.path.join(RESULTS_DIR, "amr/coverage/{sid}/{sid}_rgi_coverage.txt")
#    log:
#        os.path.join(RESULTS_DIR, "logs/rgi.coverage.{sid}.log")
#    params:
#        sample="{sid}"
#    run:
#        import os
#        import pandas as pd
#
#        # Read the RGI file
#        rgi_file = input.rgi_file
#        print(f"Reading RGI file: {rgi_file}")
#        rgi_df = pd.read_csv(rgi_file, sep='\t')
#        print("RGI DataFrame head:")
#        print(rgi_df.head())
#
#        # Read the coverage file (ignoring the first column)
#        coverage_file = input.coverage_file
#        print(f"Reading coverage file: {coverage_file}")
#        coverage_df = pd.read_csv(coverage_file, sep=' ', header=None, usecols=[1, 2], names=['Contig', 'Coverage'])
#        print("Coverage DataFrame head:")
#        print(coverage_df.head())
#
#        # Merge the two dataframes based on Contig
#        print("Merging RGI and coverage DataFrames based on Contig")
#        merged_df = pd.merge(rgi_df, coverage_df, on='Contig', how='left')
#        print("Merged DataFrame head:")
#        print(merged_df.head())
#
#        # Add the sample column
#        merged_df['Sample'] = params.sample
#        print(f"Added sample column: {params.sample}")
#
#        # Write the merged dataframe to the output file
#        output_file = output[0]
#        print(f"Writing merged DataFrame to: {output_file}")
#        merged_df.to_csv(output_file, sep='\t', index=False)
#        print("Merge complete.")

rule merge_rgi_coverage:
    input:
        rgi_file=os.path.join(RESULTS_DIR, "amr/{sid}/{sid}_rgi.txt"),
        coverage_file=os.path.join(DATA_DIR, "coverages/{sid}_gene_coverage.txt")
    output:
        os.path.join(RESULTS_DIR, "amr/coverage/{sid}/{sid}_rgi_coverage.txt")
    log:
        os.path.join(RESULTS_DIR, "logs/rgi.coverage.{sid}.log")
    params:
        sample="{sid}"
    run:
        import os
        import pandas as pd

        # Read the RGI file
        rgi_file = input.rgi_file
        print(f"Reading RGI file: {rgi_file}")
        rgi_df = pd.read_csv(rgi_file, sep='\t')
        print("RGI DataFrame head:")
        print(rgi_df.head())

        # Read the coverage file (ignoring the first column)
        coverage_file = input.coverage_file
        print(f"Reading coverage file: {coverage_file}")
        coverage_df = pd.read_csv(coverage_file, sep=' ', header=None, usecols=[1, 2], names=['Contig', 'Coverage'])
        print("Coverage DataFrame head:")
        print(coverage_df.head())

        # Merge the two dataframes based on Contig
        print("Merging RGI and coverage DataFrames based on Contig")
        merged_df = pd.merge(rgi_df, coverage_df, on='Contig', how='left')
        print("Merged DataFrame head:")
        print(merged_df.head())

        # Calculate total coverage for the entire sample
        total_coverage = coverage_df['Coverage'].sum()
        print(f"Total coverage for the entire sample: {total_coverage}")

        # Add the total coverage as a new column to the merged dataframe
        merged_df['Total_Coverage'] = total_coverage
        print("Added Total_Coverage column to the DataFrame.")

        # Add the sample column
        merged_df['Sample'] = params.sample
        print(f"Added sample column: {params.sample}")

        # Write the merged dataframe to the output file
        output_file = output[0]
        print(f"Writing merged DataFrame to: {output_file}")
        merged_df.to_csv(output_file, sep='\t', index=False)
        print("Merge complete.")

############################################
# Normalise ARG abundance by single-copy marker genes (fetchMG genome equivalents)
# Method: count proteins assigned to each of the 40 universal COG marker genes
# per sample; genome_equivalents = median across all marker COGs; then
# normalised_coverage = ARG_coverage / genome_equivalents.
rule normalize_arg_abundance:
    input:
        rgi_coverage=expand(os.path.join(RESULTS_DIR, "amr/coverage/{sid}/{sid}_rgi_coverage.txt"), sid=SAMPLES.index),
        fetchmg=expand(os.path.join(RESULTS_DIR, "fetchMG/{sid}/{sid}.faa.fetchMGs.scores"), sid=SAMPLES.index)
    output:
        normalized=os.path.join(RESULTS_DIR, "amr/normalized/arg_normalized_abundance.tsv"),
        genome_equiv=os.path.join(RESULTS_DIR, "amr/normalized/genome_equivalents.tsv")
    log:
        os.path.join(RESULTS_DIR, "logs/normalize_arg_abundance.log")
    message:
        "Normalising ARG abundance by fetchMG single-copy marker gene genome equivalents"
    run:
        import os
        import numpy as np
        import pandas as pd

        os.makedirs(os.path.join(RESULTS_DIR, "amr/normalized"), exist_ok=True)

        # ── Compute genome equivalents from fetchMG scores files ──────────────
        # fetchMGs scores file: tab-separated, columns protein_id / COG_id / score
        # Each row = one protein assigned to a universal marker COG.
        # genome_equivalents = median per-COG detection count across all 40 COGs.
        ge_rows = []
        sid_to_fetchmg = {
            os.path.basename(fp).split(".faa")[0]: fp
            for fp in input.fetchmg
        }
        for sid, fp in sid_to_fetchmg.items():
            try:
                scores = pd.read_csv(
                    fp, sep="\t", header=None,
                    names=["protein_id", "COG_id", "score"],
                    usecols=[0, 1, 2]
                )
                cog_counts = scores.groupby("COG_id")["protein_id"].count()
                ge = float(np.median(cog_counts.values)) if len(cog_counts) > 0 else np.nan
            except Exception:
                ge = np.nan
            ge_rows.append({"sample": sid, "genome_equivalents": ge})

        ge_df = pd.DataFrame(ge_rows).set_index("sample")
        ge_df.to_csv(str(output.genome_equiv), sep="\t")

        # ── Load and normalise all per-sample RGI coverage files ─────────────
        all_rgi = []
        for fp in input.rgi_coverage:
            try:
                df = pd.read_csv(fp, sep="\t")
                all_rgi.append(df)
            except Exception:
                pass

        if not all_rgi:
            pd.DataFrame().to_csv(str(output.normalized), sep="\t", index=False)
        else:
            combined = pd.concat(all_rgi, ignore_index=True)
            combined["genome_equivalents"] = combined["Sample"].map(
                ge_df["genome_equivalents"]
            )
            # ARGs per genome equivalent — the normalised abundance metric
            combined["arg_per_genome_equivalent"] = (
                combined["Coverage"] / combined["genome_equivalents"]
            )
            combined.to_csv(str(output.normalized), sep="\t", index=False)


############################################
# Getting single copy genes using fetchMGs: marker gene extraction from protein FASTA
rule fetchMG:
    input:
        faa=os.path.join(RESULTS_DIR, "proteins/{sid}.faa")
    output:
        outdir=os.path.join(RESULTS_DIR, "fetchMG/{sid}/{sid}.faa.fetchMGs.scores")
    log:
        os.path.join(RESULTS_DIR, "logs/fetchMG.{sid}.log")
    threads:
        config["fetchMG"]["threads"]
    params:
        mode="extraction",
        target="gene"
    conda:
        os.path.join(ENV_DIR, "fetchMG.yaml")
    message:
        "Marker genes: fetchMGs {wildcards.sid}"
    shell:
        # fetchMGs syntax: fetchMGs extraction <input.faa> gene <output_dir> -t <threads>
        "(date && "
        "mkdir -p $(dirname {output.outdir}) && "
        "fetchMGs {params.mode} {input.faa} {params.target} $(dirname {output.outdir}) -t {threads} && "
        "date) &> >(tee {log})"

