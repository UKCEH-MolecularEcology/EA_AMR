# EA RSN Resistome Pipeline — Revised Workflow

**Authors:** Susheel Bhanu Busi and Amy Thorpe  
**Affiliation:** Molecular Ecology Group, UK Centre for Ecology & Hydrology (UKCEH)  
**Last modified:** 2026-06-04

---

## Contents

1. [Overview](#overview)
2. [Prerequisites](#prerequisites)
3. [Configuration](#configuration)
4. [Database setup](#database-setup)
5. [Pipeline steps](#pipeline-steps)
6. [Running the pipeline](#running-the-pipeline)
7. [Key outputs](#key-outputs)
8. [Reviewer revision additions](#reviewer-revision-additions)
9. [Citations](#citations)

---

## Overview

This Snakemake workflow processes 450 freshwater biofilm metagenomes from England's river
network to characterise the antimicrobial resistance gene (ARG) resistome and its
environmental drivers. It was revised in response to peer review at *Nature Communications*
to address concerns about taxonomic accuracy, contig-level analysis, mobile genetic
elements, and metal resistance co-selection.

The pipeline produces:
- ARG annotations (CARD-RGI) with genome-equivalent-normalised abundances (fetchMG)
- Contig-level taxonomy at three length cutoffs (unfiltered, ≥2 kb, ≥10 kb)
- Kraken2 read-level taxonomy at three confidence thresholds (0.2, 0.5, 0.7)
- Taxonomy validation: Spearman correlations across Kraken2, SingleM, and Sylph
- Contig mobility classification (geNomad: chromosome / plasmid / virus)
- Metal resistance gene annotation (BacMet2) and ARG–metal co-selection analysis
- MAG generation, quality control, and ARG profiling (MetaBAT2 → CheckM2 → GTDB-Tk)
- An integrated per-contig master table joining all of the above

---

## Prerequisites

- **Snakemake:** `/hdd0/susbus/tools/conda_envs/snakemake/bin/snakemake`
- **Conda/Mamba:** `/home/susbus/miniforge3`
- **Conda environments** are built automatically on first run using the YAMLs in `envs/`

---

## Configuration

All settings live in `config/config.yaml`. Key fields:

| Field | Description |
|---|---|
| `steps` | List of pipeline steps to run (see [Pipeline steps](#pipeline-steps)) |
| `work_dir` | Root of this repository |
| `data_dir` | Pre-existing assemblies, reads, Kraken2 reports, coverages |
| `results_dir` | All pipeline outputs are written here |
| `env_dir` | Path to `workflow_revised/envs/` |
| `samples` | Path to `config/450_samples.tsv` |

**Sample manifest** (`config/450_samples.tsv`): tab-separated with columns
`Sample_ID`, `sR1`, `sR2`. The `sR1`/`sR2` columns contain absolute paths to the
human-filtered paired-end reads (`.fq.gz`) spread across:
- `/prj/DECODE/ea_biofilm_results/preprocessed/reads/{sid}/`
- `/prj/DECODE/ea_biofilm_results_AT/preprocessed/reads/{sid}/`
- `/prj/DECODE/ea_biofilm_results_AT_additional/preprocessed/reads/{sid}/`

---

## Database setup

Run these once before the main pipeline. They are `localrules` and execute on the
head/login node (no job submission needed).

```bash
BASE="snakemake -s workflow_revised/Snakefile --configfile config/config.yaml \
  --conda-frontend conda --use-conda \
  --conda-prefix /prj/DECODE/conda_envs --cores 1"

# Organellar reference (plastid + mitochondrion from NCBI RefSeq)
$BASE download_organellar_db

# BacMet2 metal resistance database → /hdd0/susbus/databases/bacmet/
$BASE download_bacmet_db

# geNomad database — pre-existing at:
# /hdd0/susbus/databases/ViWrap_db/genomad_db  (already configured in config.yaml)
```

The CARD-RGI database setup (`download_rgi_db`, `setup_rgi_db`) runs automatically as
part of the `sample_amr` step.

---

## Pipeline steps

Enable steps by adding their names to `steps` in `config/config.yaml`.

### Original steps

| Step | Rule file | Description |
|---|---|---|
| `preprocessing` | `preprocessing.smk` | Trim Galore QC, human read filtering (GRCh38), FastQC/MultiQC |
| `taxonomy` | `taxonomy.smk` + `singlem.smk` | Kraken2+Bracken read-level taxonomy; SingleM marker gene profiling |
| `assembly` | `assembly.smk` | MEGAHIT metagenome assembly |
| `annotation` | `annotation.smk` | Prodigal ORF prediction; EggNOG functional annotation |
| `coverage` | `coverage.smk` | BWA read mapping; contig and gene-level coverage (jgi_summarize) |
| `sample_amr` | `sample_amr.smk` | CARD-RGI ARG annotation (DIAMOND); fetchMG marker genes; coverage-normalised ARG abundances |
| `sylph` | `sylph.smk` | Sylph read-based profiling against GTDB r220 (prokaryotes, viruses, fungi) |
| `contig_taxonomy` | `contig_taxonomy.smk` | Kraken2+Bracken on assembled contigs; unfiltered, ≥2 kb, and ≥10 kb chains |

### Reviewer revision additions

| Step | Rule file | Addresses | Description |
|---|---|---|---|
| `organellar_filter` | `organellar_filter.smk` | R2 #5 | BBDuk filter of plastid and mitochondrion contigs before taxonomy |
| `contig_size_filter` | `contig_size_filter.smk` | R2 #4 | seqkit contig filtering at ≥2 kb (AMR) and ≥10 kb (taxonomy sensitivity) |
| `bracken_reads` | `bracken_reads.smk` | R2 #4 | Bracken on existing Kraken2 read reports (confidence=0.2) for taxonomy validation |
| `bracken_highconf` | `bracken_reads.smk` | R2 #4 | Kraken2 reruns at confidence=0.5 and 0.7 from reads; Bracken on new reports |
| `taxonomy_validation` | `taxonomy_validation.smk` | R2 #4 | Spearman correlations at phylum/class/genus: Kraken2 (3 confidence levels, 3 contig cutoffs) vs SingleM vs Sylph |
| `genomad` | `genomad.smk` | R1 #2, R1 #3 | geNomad contig classification: chromosome / plasmid / virus |
| `mge_metal` | `mge_metal.smk` | R1 #6 | ISEScan insertion sequences; BacMet2 metal resistance genes (DIAMOND); ARG+MGE co-occurrence using geNomad + ISEScan |
| `binning` | `binning.smk` | R1 #3, R2 #3 | MetaBAT2 binning; CheckM2 QC (≥50% completeness, ≤10% contamination); dRep dereplication; GTDB-Tk taxonomy |
| `bin_amr` | `bin_amr.smk` | R1 #3, R2 #3 | CARD-RGI on dereplicated MAGs; merged with GTDB-Tk taxonomy; ARG-per-phylum summary |
| `contig_summary` | `contig_summary.smk` | All | Integrated per-contig master table + ARG–metal co-selection analysis |

### Step dependency order

```
assembly ──────────────────┬──► organellar_filter
                           ├──► contig_size_filter ──► contig_taxonomy ──► taxonomy_validation
                           └──► genomad ──────────────────────────────────────────────┐
                                                                                       │
annotation + coverage ─────────────────────────────────────────────────────────────────┤
                                                                                       │
sample_amr (RGI + fetchMG) ────────────────────────────────────────────────────────────┤
                                                                                       ▼
bracken_reads ─────────────────────────────────────► taxonomy_validation     mge_metal ──► contig_summary
bracken_highconf ──────────────────────────────────► taxonomy_validation
binning ──► bin_amr
```

Snakemake resolves all dependencies automatically — list all desired steps in `config.yaml`
and run once.

---

## Running the pipeline

### Dry run (always run first)

```bash
snakemake -s workflow_revised/Snakefile \
  --configfile config/config.yaml \
  --conda-frontend conda --use-conda \
  --cores 64 \
  --conda-prefix /prj/DECODE/conda_envs \
  --jobs 8 \
  -rpkn
```

### Full run

```bash
snakemake -s workflow_revised/Snakefile \
  --configfile config/config.yaml \
  --conda-frontend conda --use-conda \
  --cores 64 \
  --conda-prefix /prj/DECODE/conda_envs \
  --jobs 8 \
  -rpk
```

**Flags:**
- `--conda-frontend conda` — use conda (not mamba) to build environments
- `--conda-prefix /prj/DECODE/conda_envs` — shared environment cache; avoids rebuilding envs between runs
- `--cores 64` — total CPU pool across all running jobs
- `--jobs 8` — maximum parallel rule executions
- `-r` — print the reason each job is run
- `-p` — print shell commands
- `-k` — keep going on independent job failures
- `-n` — dry-run only (omit for real run)

### Useful one-liners

```bash
BASE="snakemake -s workflow_revised/Snakefile --configfile config/config.yaml \
  --conda-frontend conda --use-conda \
  --conda-prefix /prj/DECODE/conda_envs"

# Generate the job DAG as a PDF
$BASE --dag -n | dot -Tpdf > dag.pdf

# Force re-run a specific rule for one sample
$BASE --cores 64 --jobs 8 -rpk --forcerun run_genomad \
  RGI_results/genomad/A21_2022/A21_2022_aggregated_classification/A21_2022_aggregated_classification.tsv

# Only build conda environments without running any jobs
$BASE --cores 8 --conda-create-envs-only

# Run only database downloads
$BASE --cores 1 download_organellar_db
$BASE --cores 1 download_bacmet_db
```

---

## Key outputs

All outputs are written to `results_dir` (default: `/prj/DECODE/EA_AMR/RGI_results/`).

### ARG analysis

| Path | Description |
|---|---|
| `amr/{sid}/{sid}_rgi.txt` | Per-sample RGI annotations (DIAMOND, CARD database) |
| `amr/coverage/{sid}/{sid}_rgi_coverage.txt` | RGI merged with gene-level coverage |
| `amr/normalized/arg_normalized_abundance.tsv` | ARGs per genome equivalent (fetchMG normalisation) |
| `amr/normalized/genome_equivalents.tsv` | Per-sample genome equivalent estimates |
| `fetchMG/{sid}/{sid}.faa.fetchMGs.scores` | Single-copy marker gene scores |

### Contig taxonomy

| Path | Description |
|---|---|
| `kraken2/contig/{sid}_kraken.{report,out}` | Kraken2 on unfiltered contigs |
| `kraken2/contig_2kb/{sid}_kraken.{report,out}` | Kraken2 on ≥2 kb contigs |
| `kraken2/contig_10kb/{sid}_kraken.{report,out}` | Kraken2 on ≥10 kb contigs (sensitivity analysis) |
| `bracken/contig*/combined_bracken.txt` | Combined Bracken outputs per cutoff |

### Read-level taxonomy

| Path | Description |
|---|---|
| `bracken/combined_bracken.txt` | Bracken on existing Kraken2 reports (c=0.2) |
| `bracken_c05/combined_bracken_c05.txt` | Bracken on Kraken2 reruns at confidence=0.5 |
| `bracken_c07/combined_bracken_c07.txt` | Bracken on Kraken2 reruns at confidence=0.7 |

### Taxonomy validation

| Path | Description |
|---|---|
| `taxonomy_validation/{rank}_comparison_wide.tsv` | Relative abundances per method × sample (rank = phylum/class/genus) |
| `taxonomy_validation/{rank}_correlations.tsv` | Pairwise Spearman r between all methods per taxon |
| `taxonomy_validation/validation_summary.tsv` | Summary: Kraken2 contigs vs SingleM and Sylph; high-confidence flag |

### geNomad

| Path | Description |
|---|---|
| `genomad/{sid}/{sid}_aggregated_classification/{sid}_aggregated_classification.tsv` | Per-contig chromosome/plasmid/virus classification |
| `genomad/combined_genomad_classification.tsv` | All samples combined |

### Metal resistance and MGEs

| Path | Description |
|---|---|
| `metal_resistance/{sid}/{sid}_bacmet.tsv` | BacMet2 DIAMOND hits per sample |
| `metal_resistance/combined_metal_resistance.tsv` | All samples combined with metal class annotations |
| `isescan/{sid}/{sid}.fasta.sum` | ISEScan insertion sequence summary |
| `mge_cooccurrence/{sid}/{sid}_mge_arg_contigs.tsv` | Per-contig MGE + ARG co-occurrence (geNomad + ISEScan) |
| `mge_cooccurrence/combined_mge_arg_summary.tsv` | Per-sample MGE-ARG co-occurrence statistics |

### Integrated contig summary

| Path | Description |
|---|---|
| `contig_summary/master_table.tsv` | One row per (sample, contig, ARG). Columns: ARG name/family/drug class/mechanism/identity/coverage/normalised abundance, metal resistance gene/class/identity, geNomad classification/scores, ISEScan IS families, Kraken2 taxonomy at 3 cutoffs, co-selection flag, ARG and metal gene counts per contig |
| `contig_summary/coselection/coselected_pairs.tsv` | All (ARG, metal gene) pairs co-occurring on the same contig |
| `contig_summary/coselection/coselection_by_classification.tsv` | Fisher's exact test: co-selection enrichment on plasmids/viruses vs chromosomes |
| `contig_summary/coselection/coselection_statistics.tsv` | ARG family × metal class pair frequencies, mean abundances, mobility fractions, top taxa |

### MAGs

| Path | Description |
|---|---|
| `bins/{sid}/` | MetaBAT2 bins per sample |
| `checkm/{sid}/quality_report.tsv` | CheckM2 completeness/contamination per bin |
| `bins_hq/` | High-quality bins (completeness ≥50%, contamination ≤10%) |
| `drep/dereplicated_genomes/` | Dereplicated MAGs (dRep) |
| `gtdbtk/classify/gtdbtk.bac120.summary.tsv` | GTDB-Tk taxonomy for all MAGs |
| `bin_amr/combined_bin_rgi.tsv` | RGI on all MAGs combined |
| `bin_amr/bin_rgi_with_taxonomy.tsv` | MAG ARGs joined with GTDB taxonomy |
| `bin_amr/arg_per_phylum_summary.tsv` | ARG family counts per GTDB phylum |

### Organellar filtering

| Path | Description |
|---|---|
| `assembly_filtered/{sid}/{sid}_noOrganellar.fasta` | Contigs with plastid/mitochondrion sequences removed |
| `assembly_filtered/{sid}/{sid}_organellar_stats.txt` | BBDuk filtering statistics |

---

## Reviewer revision additions

This workflow was extended in response to reviewers at *Nature Communications*:

| Reviewer comment | Implementation |
|---|---|
| R2 #4: Kraken2 accuracy / contig size cutoff | `contig_size_filter.smk` (≥2 kb, ≥10 kb); `bracken_reads.smk` (Kraken2 reruns at c=0.5, c=0.7); `taxonomy_validation.smk` (multi-method Spearman comparison at phylum/class/genus) |
| R2 #5: Chloroplast/mitochondria contamination | `organellar_filter.smk` (BBDuk against NCBI RefSeq plastid + mitochondrion) |
| R1 #2, R1 #3 + R2 #3: MAG-based ARG analysis | `binning.smk` (MetaBAT2 → CheckM2 → dRep → GTDB-Tk); `bin_amr.smk` (RGI on MAGs) |
| R1 #3 + R2 #3: Contig mobility | `genomad.smk` (geNomad chromosome/plasmid/virus classification) |
| R1 #6: Metal resistance co-occurrence | `mge_metal.smk` (BacMet2 + ISEScan + geNomad MGE co-occurrence) |
| All: Integrated summary | `contig_summary.smk` (master table + co-selection statistics) |

---

## Citations

If using this workflow, please cite the following tools:

**ARG annotation**
- CARD/RGI: Alcock et al., *Nucleic Acids Research* (2023). https://doi.org/10.1093/nar/gkac920

**Taxonomy**
- Kraken2: Wood et al., *Genome Biology* (2019). https://doi.org/10.1186/s13059-019-1891-0
- Bracken: Lu et al., *PeerJ Computer Science* (2017). https://doi.org/10.7717/peerj-cs.104
- SingleM / Sandpiper: Woodcroft et al., *Nature Biotechnology* (2024). https://doi.org/10.1038/s41587-024-02311-0
- Sylph: Shaw & Yu, *Nature Biotechnology* (2024). https://doi.org/10.1038/s41587-024-02412-y

**Assembly and annotation**
- MEGAHIT: Li et al., *Bioinformatics* (2015). https://doi.org/10.1093/bioinformatics/btv033
- Prodigal: Hyatt et al., *BMC Bioinformatics* (2010). https://doi.org/10.1186/1471-2105-11-119
- fetchMG: Sunagawa et al., *Nature* (2013). https://doi.org/10.1038/nature12229

**MAG analysis**
- MetaBAT2: Kang et al., *PeerJ* (2019). https://doi.org/10.7717/peerj.7359
- CheckM2: Chklovski et al., *Nature Methods* (2023). https://doi.org/10.1038/s41592-023-01940-w
- dRep: Olm et al., *ISME Journal* (2017). https://doi.org/10.1038/ismej.2017.126
- GTDB-Tk: Chaumeil et al., *Bioinformatics* (2022). https://doi.org/10.1093/bioinformatics/btac672

**Mobile elements and metal resistance**
- geNomad: Camargo et al., *Nature Biotechnology* (2023). https://doi.org/10.1038/s41587-023-01953-y
- ISEScan: Xie & Tang, *Bioinformatics* (2017). https://doi.org/10.1093/bioinformatics/btx433
- BacMet2: Pal et al., *Nucleic Acids Research* (2014). https://doi.org/10.1093/nar/gkt1252

**Workflow management**
- Snakemake: Mölder et al., *F1000Research* (2021). https://doi.org/10.12688/f1000research.29032.2

---

Developed by **Susheel Bhanu Busi** and **Amy Thorpe**, Molecular Ecology Group, UKCEH.
