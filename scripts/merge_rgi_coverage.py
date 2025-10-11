# Script for mering RGI output file with `{sid}_gene_coverage.txt` file

import pandas as pd
import os

# Accessing Snakemake variables
rgi_file = snakemake.input.rgi_file
coverage_file = snakemake.input.coverage_file
output_file = snakemake.output.merged_file
sample_name = snakemake.params.sample

def merge_files(rgi_file, coverage_file, output_file, sample_name):
    # Check if the rgi file is empty
    if os.stat(rgi_file).st_size == 0:
        print(f"{rgi_file} is empty. Creating a dummy file.")
        with open(output_file, 'w') as f:
            f.write("Sample,ORF_ID,Coverage\n")
        return

    # Load the RGI and coverage files
    rgi_df = pd.read_csv(rgi_file, sep="\t")
    coverage_df = pd.read_csv(coverage_file, sep="\s+", header=None, names=["Ignore", "ORF_ID", "Coverage"])
    
    # Drop the first column from coverage file
    coverage_df = coverage_df.drop(columns=["Ignore"])
    
    # Merge on the ORF_ID column
    merged_df = pd.merge(rgi_df, coverage_df, on="ORF_ID", how="left")
    
    # Add the sample name as a new column
    merged_df["Sample"] = sample_name

    # Reorder columns so that Sample appears first
    cols = ["Sample"] + merged_df.columns.tolist()
    cols.remove("Sample")
    merged_df = merged_df[["Sample"] + cols]

    # Write to output file
    merged_df.to_csv(output_file, sep="\t", index=False)

# Call the merge function with the variables provided by Snakemake
merge_files(rgi_file, coverage_file, output_file, sample_name)

