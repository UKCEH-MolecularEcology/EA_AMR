"""
Author: Susheel Bhanu BUSI
Affiliation: Molecular Ecology group, UKCEH
Date: [2024-10-16]
Run: snakemake -s workflow/Snakefile --use-conda --cores 64 -rp
Latest modification:
Purpose: Run CARD-RGI on dereplicated MAGs (produced by binning.smk) to provide a
         MAG-anchored complement to the contig-based resistome analysis.
         Addresses Reviewer 1 comments #3 & #4 and Reviewer 2 comment #3.
"""


############################################
# Dynamic wildcard resolution:
# MAG filenames are only known after binning+dereplication has run.
# We use glob_wildcards at parse time on the dereplicated_genomes directory.
# If dRep has not been run yet the list will be empty; Snakemake will
# correctly resolve to zero targets (no error).  Run the 'binning' step
# first, then invoke 'bin_amr' in a second snakemake call, or rely on
# the checkpoint mechanism if binning.smk exposes one.
DREP_DIR = os.path.join(RESULTS_DIR, "drep/dereplicated_genomes")
MAGS, = glob_wildcards(os.path.join(DREP_DIR, "{mag}.fa"))


############################################
# Aggregate target rule — always defined first
rule bin_amr:
    input:
        os.path.join(RESULTS_DIR, "bin_amr/bin_rgi_with_taxonomy.tsv"),
        os.path.join(RESULTS_DIR, "bin_amr/arg_per_mag_phylum.tsv")
    output:
        touch("status/bin_amr.done")


############################################
# Per-MAG RGI annotation (DNA, DIAMOND aligner)
rule bin_rgi:
    input:
        fna=os.path.join(DREP_DIR, "{mag}.fa"),
        setup="status/rgi_setup.done"
    output:
        txt=os.path.join(RESULTS_DIR, "bin_amr/{mag}/{mag}_rgi.txt")
    log:
        os.path.join(RESULTS_DIR, "logs/bin_rgi.{mag}.log")
    threads:
        config["rgi"]["threads"]
    conda:
        os.path.join(ENV_DIR, "rgi.yaml")
    message:
        "Step: bin_rgi {wildcards.mag}"
    shell:
        "(date && "
        "mkdir -p $(dirname {output.txt}) && "
        "rgi database --version --local && "
        # NOTE: https://github.com/arpcard/rgi/issues/93: KeyError: 'snp' --> re-run on failure
        "rgi main "
        "--input_sequence {input.fna} "
        "--output_file $(dirname {output.txt})/$(basename -s '.txt' {output.txt}) "
        "--local -a DIAMOND --clean --low_quality -n {threads} || "
        "rgi main "
        "--input_sequence {input.fna} "
        "--output_file $(dirname {output.txt})/$(basename -s '.txt' {output.txt}) "
        "--local -a DIAMOND --clean --low_quality -n {threads} && "
        "date) &> >(tee {log})"


############################################
# Aggregate all per-MAG RGI results into a single TSV
rule combine_bin_rgi:
    input:
        expand(os.path.join(RESULTS_DIR, "bin_amr/{mag}/{mag}_rgi.txt"), mag=MAGS)
    output:
        os.path.join(RESULTS_DIR, "bin_amr/combined_bin_rgi.tsv")
    log:
        os.path.join(RESULTS_DIR, "logs/combine_bin_rgi.log")
    message:
        "Step: combine_bin_rgi"
    run:
        import os
        import glob
        import pandas as pd

        frames = []
        # Collect every *_rgi.txt produced under bin_amr/
        rgi_pattern = os.path.join(RESULTS_DIR, "bin_amr", "*", "*_rgi.txt")
        for rgi_file in sorted(glob.glob(rgi_pattern)):
            # Derive MAG name from the filename stem (remove _rgi suffix)
            basename = os.path.basename(rgi_file)           # e.g. MAG_001_rgi.txt
            mag_name = basename.replace("_rgi.txt", "")     # e.g. MAG_001
            try:
                df = pd.read_csv(rgi_file, sep="\t")
                if df.empty:
                    continue
                df.insert(0, "MAG", mag_name)
                frames.append(df)
            except Exception as e:
                print(f"WARNING: could not parse {rgi_file}: {e}")

        if frames:
            combined = pd.concat(frames, ignore_index=True)
        else:
            combined = pd.DataFrame(columns=["MAG"])

        combined.to_csv(output[0], sep="\t", index=False)
        with open(log[0], "w") as fh:
            fh.write(f"Combined {len(frames)} RGI files into {output[0]}\n")


############################################
# Merge combined RGI table with GTDB-Tk taxonomy
rule merge_bin_rgi_taxonomy:
    input:
        rgi=os.path.join(RESULTS_DIR, "bin_amr/combined_bin_rgi.tsv"),
        gtdbtk=os.path.join(RESULTS_DIR, "gtdbtk/classify/gtdbtk.bac120.summary.tsv")
    output:
        os.path.join(RESULTS_DIR, "bin_amr/bin_rgi_with_taxonomy.tsv")
    log:
        os.path.join(RESULTS_DIR, "logs/merge_bin_rgi_taxonomy.log")
    message:
        "Step: merge_bin_rgi_taxonomy"
    run:
        import os
        import pandas as pd

        rgi_df = pd.read_csv(input.rgi, sep="\t")

        # GTDB-Tk summary: user_genome column holds the MAG stem name
        gtdb_df = pd.read_csv(input.gtdbtk, sep="\t")

        # Rename for clarity
        gtdb_df = gtdb_df.rename(columns={
            "user_genome": "MAG",
            "classification": "gtdb_classification"
        })

        # Select columns that are reliably present in every GTDB-Tk output
        gtdb_cols = ["MAG", "gtdb_classification"]

        # Optionally carry completeness from CheckM2 columns if present
        for col in ["checkm2_completeness", "checkm_completeness", "completeness"]:
            if col in gtdb_df.columns:
                gtdb_cols.append(col)
                gtdb_df = gtdb_df.rename(columns={col: "completeness"})
                break

        gtdb_sub = gtdb_df[[c for c in gtdb_cols if c in gtdb_df.columns]].copy()

        # Left-join: keep all RGI rows, add taxonomy where available
        merged = rgi_df.merge(gtdb_sub, on="MAG", how="left")

        merged.to_csv(output[0], sep="\t", index=False)
        with open(log[0], "w") as fh:
            fh.write(f"Merged {len(rgi_df)} RGI rows with GTDB-Tk taxonomy -> {len(merged)} rows\n")


############################################
# Summarise ARG occurrences per GTDB phylum
rule bin_amr_summary:
    input:
        os.path.join(RESULTS_DIR, "bin_amr/bin_rgi_with_taxonomy.tsv")
    output:
        os.path.join(RESULTS_DIR, "bin_amr/arg_per_mag_phylum.tsv")
    log:
        os.path.join(RESULTS_DIR, "logs/bin_amr_summary.log")
    message:
        "Step: bin_amr_summary"
    run:
        import os
        import re
        import pandas as pd

        df = pd.read_csv(input[0], sep="\t")

        def parse_phylum(classification):
            """
            Extract phylum from GTDB classification string of the form
            d__Bacteria;p__Pseudomonadota;c__...
            Returns the phylum token or 'Unclassified' if not parseable.
            """
            if pd.isna(classification) or str(classification).strip() == "":
                return "Unclassified"
            match = re.search(r'p__([^;]+)', str(classification))
            if match:
                phylum = match.group(1).strip()
                return phylum if phylum else "Unclassified"
            return "Unclassified"

        df["gtdb_phylum"] = df.get("gtdb_classification", pd.Series(dtype=str)).apply(parse_phylum)

        # ARG family column names vary between RGI versions; try in order of preference
        arg_family_col = None
        for candidate in ["AMR Gene Family", "Drug Class", "Resistance Mechanism", "Best_Hit_ARO"]:
            if candidate in df.columns:
                arg_family_col = candidate
                break

        if arg_family_col is None:
            # Fall back to first non-MAG, non-taxonomy column
            exclude = {"MAG", "gtdb_classification", "gtdb_phylum", "completeness"}
            fallback_cols = [c for c in df.columns if c not in exclude]
            arg_family_col = fallback_cols[0] if fallback_cols else None

        if arg_family_col:
            summary = (
                df.groupby(["gtdb_phylum", arg_family_col])
                .size()
                .reset_index(name="ARG_count")
                .sort_values(["gtdb_phylum", "ARG_count"], ascending=[True, False])
            )
            summary.columns = ["gtdb_phylum", "ARG_family", "ARG_count"]
        else:
            # No ARG family column found — emit per-phylum MAG counts only
            summary = (
                df.groupby("gtdb_phylum")["MAG"]
                .nunique()
                .reset_index(name="MAG_count")
            )

        summary.to_csv(output[0], sep="\t", index=False)
        with open(log[0], "w") as fh:
            fh.write(f"Summary written to {output[0]} using ARG family column: {arg_family_col}\n")
            fh.write(f"Phyla detected: {sorted(df['gtdb_phylum'].unique())}\n")
