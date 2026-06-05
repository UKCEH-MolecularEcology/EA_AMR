"""
Author: Susheel Bhanu BUSI
Affiliation: Molecular Ecology group, UKCEH
Date: [2024-06-04]
Run: snakemake -s workflow/rules/taxonomy_validation.smk --use-conda --cores 4 -rp
Latest modification:
Purpose: Validate Kraken2/Bracken contig taxonomy against SingleM and Sylph
         (Reviewer 2 comment #4): compute Spearman rank correlations at phylum,
         class, and genus levels between Kraken2_reads, Kraken2_contigs, SingleM,
         and Sylph across all samples.
"""

# Taxonomic ranks to compare: display name -> lineage prefix
TAXON_RANKS = {"phylum": "p__", "class": "c__", "genus": "g__"}


############################################
rule taxonomy_validation:
    input:
        expand(os.path.join(RESULTS_DIR, "taxonomy_validation/{rank}_comparison_wide.tsv"), rank=TAXON_RANKS.keys()),
        expand(os.path.join(RESULTS_DIR, "taxonomy_validation/{rank}_correlations.tsv"), rank=TAXON_RANKS.keys()),
        os.path.join(RESULTS_DIR, "taxonomy_validation/validation_summary.tsv")
    output:
        touch("status/taxonomy_validation.done")


############################################
localrules: compare_taxonomy_methods, taxonomy_validation_report


############################################
# Rule 1: aggregate taxon abundances and compute pairwise Spearman correlations
#         at phylum, class, and genus levels
rule compare_taxonomy_methods:
    input:
        bracken_reads=os.path.join(RESULTS_DIR, "bracken/combined_bracken.txt"),
        bracken_contigs=os.path.join(RESULTS_DIR, "bracken/contig/combined_bracken.txt"),
        singlem=os.path.join(RESULTS_DIR, "singlem/combined_singlem_relab.csv"),
        sylph=os.path.join(RESULTS_DIR, "sylph/profiling_prokaryotes.tsv")
    output:
        wide=expand(os.path.join(RESULTS_DIR, "taxonomy_validation/{rank}_comparison_wide.tsv"), rank=TAXON_RANKS.keys()),
        corr=expand(os.path.join(RESULTS_DIR, "taxonomy_validation/{rank}_correlations.tsv"), rank=TAXON_RANKS.keys())
    message:
        "Step: compare taxonomy methods at phylum, class, and genus levels"
    run:
        import os
        import re
        import numpy as np
        import pandas as pd
        from scipy.stats import spearmanr

        os.makedirs(os.path.join(RESULTS_DIR, "taxonomy_validation"), exist_ok=True)

        # ------------------------------------------------------------------ #
        # Helper: extract a taxonomic name at the requested rank from a lineage
        # prefix examples: "p__" (phylum), "c__" (class), "g__" (genus)
        # Handles GTDB-style "d__Bacteria;p__X;c__Y;g__Z" and pipe-separated
        # ------------------------------------------------------------------ #
        def extract_rank(lineage, prefix):
            if not isinstance(lineage, str):
                return None
            for token in re.split(r"[;|,]", lineage):
                token = token.strip()
                if token.startswith(prefix):
                    name = token[len(prefix):].strip()
                    return name if name else None
            m = re.search(re.escape(prefix) + r"([^|;,]+)", lineage)
            if m:
                return m.group(1).strip() or None
            return None

        # ------------------------------------------------------------------ #
        # Parse combined Bracken file — returns index=sample, cols=taxon names
        # ------------------------------------------------------------------ #
        def parse_combined_bracken(path, prefix):
            try:
                df = pd.read_csv(path, sep="\t")
            except Exception:
                return pd.DataFrame()
            frac_cols = [c for c in df.columns if c.endswith("_frac")]
            tax_src = "name" if "name" in df.columns else ("taxonomy" if "taxonomy" in df.columns else None)
            if not frac_cols or tax_src is None:
                return pd.DataFrame()
            samples_bc = [c.replace("_frac", "") for c in frac_cols]
            df["_taxon"] = df[tax_src].apply(lambda x: extract_rank(x, prefix))
            df = df.dropna(subset=["_taxon"])
            result = {}
            for sid, fc in zip(samples_bc, frac_cols):
                result[sid] = df.groupby("_taxon")[fc].sum()
            out = pd.DataFrame(result).T
            out.index.name = "sample"
            return out

        # ------------------------------------------------------------------ #
        # Parse SingleM combined relative abundance CSV
        # ------------------------------------------------------------------ #
        def parse_singlem(path, prefix):
            try:
                df = pd.read_csv(path, sep=None, engine="python")
            except Exception:
                return pd.DataFrame()
            tax_col = df.columns[0]
            df["_taxon"] = df[tax_col].apply(lambda x: extract_rank(x, prefix))
            df = df.dropna(subset=["_taxon"])
            sample_cols = [c for c in df.columns if c not in [tax_col, "_taxon"]]
            result = {}
            for sid in sample_cols:
                result[sid] = df.groupby("_taxon")[sid].sum()
            out = pd.DataFrame(result).T
            out.index.name = "sample"
            row_sums = out.sum(axis=1).replace(0, np.nan)
            out = out.div(row_sums, axis=0).fillna(0)
            return out

        # ------------------------------------------------------------------ #
        # Parse Sylph profiling TSV (GTDB r220)
        # ------------------------------------------------------------------ #
        def parse_sylph(path, prefix):
            try:
                df = pd.read_csv(path, sep="\t")
            except Exception:
                return pd.DataFrame()
            tax_col = next(
                (c for c in ["clade_name", "lineage", "Lineage", "taxonomy",
                              "Taxonomic_lineage"] if c in df.columns),
                df.columns[-1]
            )
            abund_col = next(
                (c for c in ["taxonomic_abundance", "Taxonomic_abundance",
                              "relative_abundance", "sequence_abundance"] if c in df.columns),
                df.columns[2]
            )
            sample_col = next(
                (c for c in ["sample_file", "Sample_file", "sample", "Sample"] if c in df.columns),
                None
            )
            df["_taxon"] = df[tax_col].apply(lambda x: extract_rank(x, prefix))
            df = df.dropna(subset=["_taxon"])
            df[abund_col] = pd.to_numeric(df[abund_col], errors="coerce").fillna(0)
            if sample_col is not None:
                df["_sid"] = df[sample_col].apply(
                    lambda x: os.path.splitext(os.path.basename(str(x)))[0]
                )
                result = {}
                for sid, grp in df.groupby("_sid"):
                    taxon_abund = grp.groupby("_taxon")[abund_col].sum()
                    total = taxon_abund.sum()
                    result[sid] = taxon_abund / total if total > 0 else taxon_abund
                out = pd.DataFrame(result).T
            else:
                taxon_abund = df.groupby("_taxon")[abund_col].sum()
                total = taxon_abund.sum()
                out = (taxon_abund / total if total > 0 else taxon_abund).to_frame(name="sylph_sample").T
            out.index.name = "sample"
            return out

        # ------------------------------------------------------------------ #
        # Run the full comparison for each taxonomic rank
        # ------------------------------------------------------------------ #
        rank_items = list(TAXON_RANKS.items())  # [("phylum","p__"), ("class","c__"), ("genus","g__")]

        for rank_name, rank_prefix in rank_items:
            wide_path = os.path.join(RESULTS_DIR, f"taxonomy_validation/{rank_name}_comparison_wide.tsv")
            corr_path = os.path.join(RESULTS_DIR, f"taxonomy_validation/{rank_name}_correlations.tsv")

            dfs = {
                "Kraken2_reads":   parse_combined_bracken(str(input.bracken_reads),   rank_prefix),
                "Kraken2_contigs": parse_combined_bracken(str(input.bracken_contigs), rank_prefix),
                "SingleM":         parse_singlem(str(input.singlem), rank_prefix),
                "Sylph":           parse_sylph(str(input.sylph), rank_prefix),
            }

            # Align all DataFrames to the union of taxon names
            all_taxa = sorted({t for df in dfs.values() if not df.empty for t in df.columns})
            for method in list(dfs.keys()):
                if not dfs[method].empty:
                    dfs[method] = dfs[method].reindex(columns=all_taxa, fill_value=0)

            # Write wide table: multi-index (method, sample) × taxon
            rows = []
            for method, df in dfs.items():
                if df.empty:
                    continue
                df2 = df.copy()
                df2.index = pd.MultiIndex.from_tuples(
                    [(method, s) for s in df2.index], names=["method", "sample"]
                )
                rows.append(df2)
            (pd.concat(rows, axis=0) if rows else pd.DataFrame()).to_csv(wide_path, sep="\t")

            # Compute pairwise Spearman correlations per taxon
            method_names = [m for m, df in dfs.items() if not df.empty]
            pairs = [(m1, m2) for i, m1 in enumerate(method_names) for m2 in method_names[i+1:]]
            corr_records = []
            for taxon in all_taxa:
                for (m1, m2) in pairs:
                    df1, df2 = dfs.get(m1), dfs.get(m2)
                    if df1 is None or df1.empty or df2 is None or df2.empty:
                        continue
                    if taxon not in df1.columns or taxon not in df2.columns:
                        continue
                    common = df1.index.intersection(df2.index)
                    if len(common) < 3:
                        continue
                    v1 = df1.loc[common, taxon].values.astype(float)
                    v2 = df2.loc[common, taxon].values.astype(float)
                    if v1.sum() == 0 and v2.sum() == 0:
                        continue
                    try:
                        r, pval = spearmanr(v1, v2)
                    except Exception:
                        r, pval = np.nan, np.nan
                    corr_records.append({
                        "rank": rank_name,
                        "taxon": taxon,
                        "method_A": m1,
                        "method_B": m2,
                        "spearman_r": round(float(r), 4) if not np.isnan(r) else np.nan,
                        "pvalue": round(float(pval), 6) if not np.isnan(pval) else np.nan,
                        "n_samples": len(common)
                    })
            pd.DataFrame(
                corr_records,
                columns=["rank", "taxon", "method_A", "method_B", "spearman_r", "pvalue", "n_samples"]
            ).to_csv(corr_path, sep="\t", index=False)


############################################
# Rule 2: write reviewer-facing summary TSV across all three ranks
rule taxonomy_validation_report:
    input:
        corrs=expand(os.path.join(RESULTS_DIR, "taxonomy_validation/{rank}_correlations.tsv"), rank=TAXON_RANKS.keys())
    output:
        summary=os.path.join(RESULTS_DIR, "taxonomy_validation/validation_summary.tsv")
    message:
        "Step: write taxonomy validation summary at phylum, class, and genus levels"
    run:
        import pandas as pd
        import numpy as np

        all_rows = []
        for corr_path in input.corrs:
            try:
                corr_df = pd.read_csv(corr_path, sep="\t")
            except Exception:
                continue

            # rank column was added in compare_taxonomy_methods; fall back to filename
            if "rank" not in corr_df.columns:
                rank_name = os.path.basename(corr_path).replace("_correlations.tsv", "")
                corr_df["rank"] = rank_name

            taxon_col = "taxon" if "taxon" in corr_df.columns else "phylum"

            for (rank_name, taxon), grp in corr_df.groupby(["rank", taxon_col]):
                row = {"rank": rank_name, "taxon": taxon}

                def _mean_r(ma, mb):
                    mask = (
                        ((grp["method_A"] == ma) & (grp["method_B"] == mb)) |
                        ((grp["method_A"] == mb) & (grp["method_B"] == ma))
                    )
                    vals = grp.loc[mask, "spearman_r"].dropna()
                    return round(float(vals.mean()), 4) if len(vals) > 0 else np.nan

                row["spearman_r_contigs_vs_SingleM"] = _mean_r("Kraken2_contigs", "SingleM")
                row["spearman_r_contigs_vs_Sylph"]   = _mean_r("Kraken2_contigs", "Sylph")
                row["spearman_r_reads_vs_contigs"]   = _mean_r("Kraken2_reads",   "Kraken2_contigs")
                # high_confidence: both key comparisons r > 0.7
                r_sm = row["spearman_r_contigs_vs_SingleM"]
                r_sy = row["spearman_r_contigs_vs_Sylph"]
                row["high_confidence"] = (
                    not np.isnan(r_sm) and not np.isnan(r_sy) and r_sm > 0.7 and r_sy > 0.7
                )
                all_rows.append(row)

        summary_df = pd.DataFrame(all_rows, columns=[
            "rank", "taxon",
            "spearman_r_contigs_vs_SingleM",
            "spearman_r_contigs_vs_Sylph",
            "spearman_r_reads_vs_contigs",
            "high_confidence"
        ]).sort_values(["rank", "taxon"]).reset_index(drop=True)

        summary_df.to_csv(str(output.summary), sep="\t", index=False)
