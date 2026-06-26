"""
Author: Susheel Bhanu BUSI
Affiliation: Molecular Ecology group, UKCEH
Date: [2024-10-16]
Run: snakemake -s workflow/Snakefile --use-conda --cores 64 -rp
Latest modification: 2025-06-04
Purpose: Run CARD-RGI on dereplicated MAGs (produced by binning.smk) and merge
         results with GTDB-Tk taxonomy to deliver host-resolved ARG analysis.
         Addresses Reviewer 1 comments #3 & #4 and Reviewer 2 comment #3.

Dependencies (from binning.smk):
  - Dereplicated MAGs: RESULTS_DIR/drep/dereplicated_genomes/*.fa
  - GTDB-Tk summary:   RESULTS_DIR/gtdbtk/classify/gtdbtk.bac120.summary.tsv
  - RGI DB loaded:     status/rgi_setup.done + DB_DIR/rgi/card.json (sample_amr.smk)

NOTE: MAG filenames are only known after dRep has completed.  We therefore use a
checkpoint + glob_wildcards pattern so that Snakemake re-evaluates the DAG once
dereplicated_genomes/ is populated.
"""

import glob as _glob


############################################
# Helper: dereplicated genomes directory
DREP_DIR = os.path.join(RESULTS_DIR, "drep/dereplicated_genomes")


############################################
# Helper function: collect per-MAG RGI outputs after checkpoint resolves
def get_mag_rgi_outputs(wildcards):
    """
    Wait for checkpoint_drep_complete, then discover all *.fa MAGs via
    glob_wildcards and return the expected per-MAG RGI output paths.
    """
    # Trigger checkpoint — raises IncompleteCheckpointException if not done yet
    checkpoints.checkpoint_drep_complete.get()

    # Discover MAG stems from the dereplicated_genomes directory
    mag_files = _glob.glob(os.path.join(DREP_DIR, "*.fa"))
    mag_stems = [os.path.splitext(os.path.basename(f))[0] for f in mag_files]

    return expand(
        os.path.join(RESULTS_DIR, "bin_amr/{mag}/{mag}_rgi.txt"),
        mag=mag_stems
    )


############################################
# 0. Aggregate target rule — always defined first
rule bin_amr:
    input:
        get_mag_rgi_outputs,
        os.path.join(RESULTS_DIR, "bin_amr/combined_bin_rgi.tsv"),
        os.path.join(RESULTS_DIR, "bin_amr/bin_rgi_with_taxonomy.tsv"),
        os.path.join(RESULTS_DIR, "bin_amr/arg_per_phylum_summary.tsv")
    output:
        touch("status/bin_amr.done")


############################################
# 1. Checkpoint: signal that dRep has finished and MAG files are in place
checkpoint checkpoint_drep_complete:
    input:
        os.path.join(RESULTS_DIR, "drep/drep.done")
    output:
        flag=touch(os.path.join(RESULTS_DIR, "drep/checkpoint.done"))


############################################
# 2. Per-MAG RGI annotation (DNA, DIAMOND aligner)
rule bin_rgi:
    input:
        fna=os.path.join(DREP_DIR, "{mag}.fa"),
        db=os.path.join(DB_DIR, "rgi/card.json"),
        setup="status/rgi_setup.done"
    output:
        txt=os.path.join(RESULTS_DIR, "bin_amr/{mag}/{mag}_rgi.txt")
    log:
        os.path.join(RESULTS_DIR, "logs/bin_rgi.{mag}.log")
    threads:
        config["rgi"]["threads"]
    params:
        alignment_tool="DIAMOND"
    conda:
        os.path.join(ENV_DIR, "rgi.yaml")
    message:
        "Step: bin_rgi — RGI on MAG {wildcards.mag}"
    shell:
        "(date && "
        "mkdir -p $(dirname {output.txt}) && "
        "rgi database --version --local && "
        # NOTE: https://github.com/arpcard/rgi/issues/93 KeyError: 'snp' -> retry once
        "rgi main "
        "--input_sequence {input.fna} "
        "--output_file $(dirname {output.txt})/$(basename -s '.txt' {output.txt}) "
        "--local -a {params.alignment_tool} --clean --low_quality -n {threads} || "
        "rgi main "
        "--input_sequence {input.fna} "
        "--output_file $(dirname {output.txt})/$(basename -s '.txt' {output.txt}) "
        "--local -a {params.alignment_tool} --clean --low_quality -n {threads} && "
        "date) &> >(tee {log})"


############################################
# 3. Combine all per-MAG RGI tables into one TSV
rule combine_bin_rgi:
    input:
        get_mag_rgi_outputs
    output:
        os.path.join(RESULTS_DIR, "bin_amr/combined_bin_rgi.tsv")
    log:
        os.path.join(RESULTS_DIR, "logs/combine_bin_rgi.log")
    message:
        "Step: combine_bin_rgi — merging all per-MAG RGI results"
    run:
        import os
        import pandas as pd

        frames = []
        for rgi_file in sorted(input):
            # Derive MAG stem: strip _rgi.txt suffix from filename
            basename = os.path.basename(rgi_file)       # e.g. MAG_001_rgi.txt
            mag_name = basename.replace("_rgi.txt", "") # e.g. MAG_001
            try:
                df = pd.read_csv(rgi_file, sep="\t")
                if df.empty:
                    continue
                df.insert(0, "MAG", mag_name)
                frames.append(df)
            except Exception as exc:
                print(f"WARNING: could not parse {rgi_file}: {exc}")

        if frames:
            combined = pd.concat(frames, ignore_index=True)
        else:
            combined = pd.DataFrame(columns=["MAG"])

        combined.to_csv(output[0], sep="\t", index=False)

        with open(log[0], "w") as fh:
            fh.write(f"Combined {len(frames)} RGI files into {output[0]}\n")


############################################
# 4. Merge combined RGI table with GTDB-Tk taxonomy
rule merge_bin_rgi_taxonomy:
    input:
        rgi=os.path.join(RESULTS_DIR, "bin_amr/combined_bin_rgi.tsv"),
        tax=os.path.join(RESULTS_DIR, "gtdbtk/classify/gtdbtk.bac120.summary.tsv")
    output:
        os.path.join(RESULTS_DIR, "bin_amr/bin_rgi_with_taxonomy.tsv")
    log:
        os.path.join(RESULTS_DIR, "logs/merge_bin_rgi_taxonomy.log")
    message:
        "Step: merge_bin_rgi_taxonomy — joining RGI results with GTDB-Tk classification"
    run:
        import os
        import re
        import pandas as pd

        rgi_df = pd.read_csv(input.rgi, sep="\t")

        # GTDB-Tk summary: 'user_genome' holds the MAG stem name
        gtdb_df = pd.read_csv(input.tax, sep="\t")

        # Keep only the columns needed for host-resolved analysis
        keep_cols = ["user_genome", "classification"]
        if "fastani_reference_radius" in gtdb_df.columns:
            keep_cols.append("fastani_reference_radius")

        gtdb_sub = gtdb_df[keep_cols].copy()
        gtdb_sub = gtdb_sub.rename(columns={"user_genome": "MAG"})

        # Parse phylum from GTDB classification string
        # e.g. "d__Bacteria;p__Pseudomonadota;c__..." -> "Pseudomonadota"
        def _parse_phylum(classification):
            if pd.isna(classification) or str(classification).strip() == "":
                return "Unclassified"
            parts = str(classification).split(";")
            for part in parts:
                part = part.strip()
                if part.startswith("p__"):
                    phylum = part[3:].strip()
                    return phylum if phylum else "Unclassified"
            return "Unclassified"

        gtdb_sub["gtdb_phylum"] = gtdb_sub["classification"].apply(_parse_phylum)

        # Left-join: retain all RGI rows, annotate with taxonomy where available
        merged = rgi_df.merge(gtdb_sub, on="MAG", how="left")

        merged.to_csv(output[0], sep="\t", index=False)

        with open(log[0], "w") as fh:
            fh.write(
                f"Merged {len(rgi_df)} RGI rows with GTDB-Tk taxonomy "
                f"-> {len(merged)} rows written to {output[0]}\n"
            )


############################################
# 5. Summarise ARG families per GTDB phylum (pivot: rows=phyla, cols=ARG families)
rule bin_amr_phylum_summary:
    input:
        os.path.join(RESULTS_DIR, "bin_amr/bin_rgi_with_taxonomy.tsv")
    output:
        os.path.join(RESULTS_DIR, "bin_amr/arg_per_phylum_summary.tsv")
    log:
        os.path.join(RESULTS_DIR, "logs/bin_amr_phylum_summary.log")
    message:
        "Step: bin_amr_phylum_summary — ARG family counts per GTDB phylum"
    run:
        import os
        import pandas as pd

        df = pd.read_csv(input[0], sep="\t")

        # Ensure phylum column is present (may have been added in merge step)
        if "gtdb_phylum" not in df.columns:
            if "classification" in df.columns:
                def _parse_phylum(classification):
                    if pd.isna(classification) or str(classification).strip() == "":
                        return "Unclassified"
                    parts = str(classification).split(";")
                    for part in parts:
                        part = part.strip()
                        if part.startswith("p__"):
                            phylum = part[3:].strip()
                            return phylum if phylum else "Unclassified"
                    return "Unclassified"
                df["gtdb_phylum"] = df["classification"].apply(_parse_phylum)
            else:
                df["gtdb_phylum"] = "Unclassified"

        # Resolve ARG Gene Family column (RGI column names vary across versions)
        arg_family_col = None
        for candidate in ["AMR Gene Family", "Drug Class", "Resistance Mechanism", "Best_Hit_ARO"]:
            if candidate in df.columns:
                arg_family_col = candidate
                break

        if arg_family_col is None:
            exclude = {"MAG", "classification", "gtdb_phylum", "gtdb_classification",
                       "fastani_reference_radius", "completeness"}
            fallback = [c for c in df.columns if c not in exclude]
            arg_family_col = fallback[0] if fallback else None

        if arg_family_col:
            # Count unique MAGs carrying each ARG family within each phylum
            grouped = (
                df.groupby(["gtdb_phylum", arg_family_col])["MAG"]
                .nunique()
                .reset_index(name="MAG_count")
            )

            # Pivot: rows = phyla, columns = ARG families, values = unique MAG count
            pivot = grouped.pivot_table(
                index="gtdb_phylum",
                columns=arg_family_col,
                values="MAG_count",
                aggfunc="sum",
                fill_value=0
            )

            # Sort rows by total ARG count (descending)
            pivot["__total__"] = pivot.sum(axis=1)
            pivot = pivot.sort_values("__total__", ascending=False).drop(columns="__total__")
            pivot.columns.name = None
            pivot = pivot.reset_index()
            pivot = pivot.rename(columns={"gtdb_phylum": "Phylum"})
        else:
            # No ARG family column found — fall back to per-phylum MAG counts
            pivot = (
                df.groupby("gtdb_phylum")["MAG"]
                .nunique()
                .reset_index(name="MAG_count")
                .rename(columns={"gtdb_phylum": "Phylum"})
                .sort_values("MAG_count", ascending=False)
            )

        pivot.to_csv(output[0], sep="\t", index=False)

        with open(log[0], "w") as fh:
            fh.write(
                f"Summary written to {output[0]}\n"
                f"ARG family column used: {arg_family_col}\n"
                f"Phyla detected: {sorted(df['gtdb_phylum'].unique())}\n"
            )
