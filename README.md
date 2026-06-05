# EA RSN: Freshwater Biofilm Resistome — Snakemake Pipeline

**Authors:** Susheel Bhanu Busi and Amy Thorpe  
**Affiliation:** Molecular Ecology Group, UK Centre for Ecology & Hydrology (UKCEH)  
**Last modified:** 2026-06-04

---

## Overview

This repository contains Snakemake workflows for the metagenomic analysis of antimicrobial
resistance genes (ARGs) in freshwater biofilm samples collected across England's river
network (~450 samples), as described in:

> Spurr et al. *"Environmental drivers shape the freshwater biofilm resistome across a
> national river network"* — submitted to *Nature Communications*.

The **revised workflow** (`workflow_revised/`) incorporates all bioinformatics additions
requested during peer review. The original workflow (`workflow/`) is retained for reference.

---

## Repository structure

```
EA_AMR/
├── workflow_revised/          # Revised pipeline (use this)
│   ├── Snakefile
│   ├── rules/                 # All Snakemake rule files
│   ├── envs/                  # Conda environment YAMLs
│   └── scripts/               # Helper scripts
├── workflow/                  # Original pipeline (reference only)
├── config/
│   ├── config.yaml            # Main configuration file — edit steps here
│   ├── 450_samples.tsv        # Sample manifest with absolute read paths
│   └── slurm.yaml             # (unused — no SLURM on this machine)
├── schemas/                   # Config/sample validation schemas
├── RGI_results/               # Pipeline outputs
└── paper_analysis/            # Manuscript and reviewer notes
```

---

## Quick start

### 1. Configure steps

Edit `config/config.yaml` — set the `steps` list to the analyses you want to run:

```yaml
steps: [
  "sample_amr",
  "organellar_filter", "contig_size_filter", "contig_taxonomy",
  "bracken_reads", "bracken_highconf", "taxonomy_validation",
  "genomad", "mge_metal", "contig_summary"
]
```

### 2. Dry run

```bash
snakemake -s workflow_revised/Snakefile \
  --configfile config/config.yaml \
  --conda-frontend conda --use-conda \
  --cores 64 \
  --conda-prefix /prj/DECODE/conda_envs \
  --jobs 8 \
  -rpkn
```

### 3. Full run

```bash
snakemake -s workflow_revised/Snakefile \
  --configfile config/config.yaml \
  --conda-frontend conda --use-conda \
  --cores 64 \
  --conda-prefix /prj/DECODE/conda_envs \
  --jobs 8 \
  -rpk
```

---

## Detailed documentation

See [`workflow_revised/README.md`](workflow_revised/README.md) for full documentation
including all pipeline steps, database setup, output descriptions, and citations.
