# EA RSN: AMR and Coverage Snakemake Workflows

**Author:** Susheel Bhanu BUSI  
**Affiliation:** Molecular Ecology Group, UK Centre for Ecology & Hydrology (UKCEH)  
**Date:** 2024-10-16  
**Last Modified:** 2024-10-16  

---

## Overview
This repository contains two Snakemake workflows designed to process metagenomic assemblies for antimicrobial resistance (AMR) and gene coverage estimation.  
The workflows automate the steps required to annotate assemblies with the CARD-RGI tool, calculate contig and gene coverage, and merge annotation results with coverage statistics.

---

## Workflow 1: AMR (`workflow/rules/sample_amr.smk`)

**Purpose:**  
To identify and annotate antimicrobial resistance genes (ARGs) using CARD-RGI, estimate marker gene presence using fetchMGs, and merge coverage results.

**Run command:**
```bash
snakemake -s workflow/rules/amr.smk --use-conda --cores 64 -rp
```

### Key Rules
1. **download_rgi_db** – Downloads the CARD RGI database.  
2. **setup_rgi_db** – Loads and verifies the local RGI database.  
3. **sample_annotation_rgi** – Runs RGI on individual assemblies to annotate ARGs.  
4. **merge_rgi_coverage** – Merges RGI output with gene coverage information.  
5. **arg_prodigal** – Predicts open reading frames using Prodigal.  
6. **fetchMG** – Extracts marker genes using fetchMGs.

---

## Workflow 2: Coverage (`workflow/rules/coverage.smk`)

**Purpose:**  
To calculate contig- and gene-level coverage from metagenomic assemblies.

**Run command:**
```bash
# dry-run
snakemake -s workflow/Snakefile --cores 72 --jobs 6 --use-conda --conda-frontend conda -rpn

# full run
snakemake -s workflow/Snakefile --cores 72 --jobs 6 --use-conda --conda-frontend conda -rp 
```

### Key Rules
1. **mapping_index** – Indexes assemblies for read mapping using BWA.  
2. **mapping** – Maps preprocessed reads to assemblies and generates sorted BAM files.  
3. **summarise_depth** – Estimates coverage depth using `jgi_summarize_bam_contig_depths`.  
4. **contig_length** – Calculates contig lengths using `fastaNamesSizes.pl`.  
5. **gene_depth** – Calculates gene-level depth and length.  
6. **contig_gene_link** – Links genes to contigs and coverage values.

---

## Environment Setup
Both workflows require Conda environments defined in the `ENV_DIR` path.  
Each rule specifies its own environment YAML file.  
Ensure that all dependencies (BWA, samtools, RGI, Prodigal, fetchMGs, etc.) are properly installed or accessible via conda.

---

## Input and Output Structure

**Input directories:**
- `DATA_DIR`: Contains assemblies and preprocessed reads.  
- `DB_DIR`: Contains RGI database files.  

**Output directories:**
- `RESULTS_DIR/amr/`: RGI annotations and merged results.  
- `RESULTS_DIR/coverage/`: Coverage summaries and statistics.  
- `RESULTS_DIR/fetchMG/`: Marker gene extraction results.  

---

## Status Files
Each workflow generates status files to ensure reproducibility:

- `status/rgi_setup.done` – Indicates successful RGI setup.  
- `status/sample_amr.done` – Marks AMR workflow completion.  
- `status/coverage.done` – Marks coverage workflow completion.

---

## Logging
All steps produce detailed log files stored under:  
`RESULTS_DIR/logs/`  
These logs contain timestamps and command outputs for reproducibility and troubleshooting.

---

## Citation and Acknowledgement
If using this workflow, please cite:

- CARD RGI: Alcock et al., *Nucleic Acids Research* (2020).  
- fetchMGs: Sunagawa et al., *Nature* (2013).  
- Prodigal: Hyatt et al., *BMC Bioinformatics* (2010).  

Developed by **Susheel Bhanu BUSI**, Molecular Ecology Group, UKCEH.
