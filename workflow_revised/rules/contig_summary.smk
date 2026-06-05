"""
Author: Susheel Bhanu BUSI
Affiliation: Molecular Ecology group, UKCEH
Date: [2024-06-04]
Run: snakemake -s workflow_revised/Snakefile --configfile config/config.yaml --use-conda -rp
Latest modification:
Purpose: Build an integrated per-contig master table joining:
           - RGI ARG annotations + per-gene coverage
           - BacMet2 metal resistance gene hits
           - geNomad contig classification (chromosome/plasmid/virus)
           - ISEScan insertion sequence elements
           - Kraken2 taxonomy at unfiltered, 2 kb, and 10 kb contig cutoffs
           - fetchMG-normalised ARG abundance (ARGs per genome equivalent)
         Then run a co-selection analysis to identify ARG + metal resistance
         gene pairs co-occurring on the same contig, with their taxonomy,
         abundance, and mobility context.

         Master table columns:
           sample, contig_id, contig_length,
           ARG_name, ARG_family, ARG_drug_class, ARG_resistance_mechanism,
           ARG_pct_identity, ARG_coverage, ARG_normalized_abundance,
           n_args_on_contig,
           metal_resistance_gene, metal_class, metal_pct_identity,
           n_metal_genes_on_contig,
           coselected,
           genomad_classification, genomad_plasmid_score, genomad_virus_score,
           IS_element_families,
           kraken_taxid_raw, kraken_status_raw,
           kraken_taxid_2kb, kraken_status_2kb,
           kraken_taxid_10kb, kraken_status_10kb
"""


############################################
rule contig_summary:
    input:
        os.path.join(RESULTS_DIR, "contig_summary/master_table.tsv"),
        os.path.join(RESULTS_DIR, "contig_summary/coselection/coselected_pairs.tsv"),
        os.path.join(RESULTS_DIR, "contig_summary/coselection/coselection_by_classification.tsv"),
        os.path.join(RESULTS_DIR, "contig_summary/coselection/coselection_statistics.tsv")
    output:
        touch("status/contig_summary.done")


############################################
localrules: build_contig_master_table, coselection_analysis


############################################
# ── Per-sample helper: parse Kraken2 .out file → contig_id → (status, taxid)
# The Kraken2 output format is:
#   C/U  <seq_id>  <taxid>  <length>  <kmer_data>
# ─────────────────────────────────────────────────────────────────────────────

rule build_contig_master_table:
    input:
        # RGI coverage (has Coverage + Total_Coverage + Sample + all RGI cols)
        rgi_cov=expand(
            os.path.join(RESULTS_DIR, "amr/coverage/{sid}/{sid}_rgi_coverage.txt"),
            sid=SAMPLES.index
        ),
        # Normalised ARG abundance
        normalized=os.path.join(RESULTS_DIR, "amr/normalized/arg_normalized_abundance.tsv"),
        # BacMet2 metal resistance hits
        bacmet=expand(
            os.path.join(RESULTS_DIR, "metal_resistance/{sid}/{sid}_bacmet.tsv"),
            sid=SAMPLES.index
        ),
        # geNomad classifications
        genomad=expand(
            os.path.join(
                RESULTS_DIR,
                "genomad/{sid}/{sid}_aggregated_classification/{sid}_aggregated_classification.tsv"
            ),
            sid=SAMPLES.index
        ),
        # ISEScan summaries
        isescan=expand(
            os.path.join(RESULTS_DIR, "isescan/{sid}/assembly/{sid}.fasta.sum"),
            sid=SAMPLES.index
        ),
        # Kraken2 contig-level .out files at three length cutoffs
        k2_raw=expand(
            os.path.join(RESULTS_DIR, "kraken2/contig/{sid}_kraken.out"),
            sid=SAMPLES.index
        ),
        k2_2kb=expand(
            os.path.join(RESULTS_DIR, "kraken2/contig_2kb/{sid}_kraken.out"),
            sid=SAMPLES.index
        ),
        k2_10kb=expand(
            os.path.join(RESULTS_DIR, "kraken2/contig_10kb/{sid}_kraken.out"),
            sid=SAMPLES.index
        )
    output:
        master=os.path.join(RESULTS_DIR, "contig_summary/master_table.tsv")
    log:
        os.path.join(RESULTS_DIR, "logs/contig_summary/build_master_table.log")
    message:
        "Contig summary: building integrated master table"
    run:
        import os
        import re
        import numpy as np
        import pandas as pd

        os.makedirs(os.path.join(RESULTS_DIR, "contig_summary"), exist_ok=True)

        # ── helpers ──────────────────────────────────────────────────────────

        BACMET_COLS = ["qseqid", "sseqid", "pident", "length", "mismatch",
                       "gap", "qs", "qe", "ss", "se", "evalue", "bitscore"]

        METAL_CLASS = {
            "mera": "mercury",  "merb": "mercury",  "merc": "mercury",
            "merd": "mercury",  "mere": "mercury",  "merp": "mercury",
            "mert": "mercury",
            "copa": "copper",   "copb": "copper",   "copc": "copper",
            "cusa": "copper",   "cusb": "copper",
            "arsa": "arsenic",  "arsb": "arsenic",  "arsc": "arsenic",
            "arsd": "arsenic",
            "cada": "cadmium",  "cadc": "cadmium",
            "chra": "chromate", "chrb": "chromate",
            "znta": "zinc",     "zntb": "zinc",      "znur": "zinc",
            "silp": "silver",   "sils": "silver",
            "terb": "tellurite","terc": "tellurite",
        }

        def assign_metal_class(sseqid):
            sl = sseqid.lower()
            for prefix, cls in METAL_CLASS.items():
                if prefix in sl:
                    return cls
            return "other"

        def parse_kraken_out(path):
            """Return dict contig_id -> (status, taxid)."""
            result = {}
            try:
                df = pd.read_csv(
                    path, sep="\t", header=None,
                    names=["status", "seq_id", "taxid", "length", "kmer"],
                    usecols=[0, 1, 2]
                )
                for _, row in df.iterrows():
                    result[str(row["seq_id"])] = (row["status"], int(row["taxid"]))
            except Exception:
                pass
            return result

        def parse_isescan(path):
            """Return dict contig_id -> comma-joined IS families."""
            result = {}
            try:
                df = pd.read_csv(path, sep=r"\s+", comment="#", engine="python")
                if "seqID" in df.columns and "isFamily" in df.columns:
                    for cid, grp in df.groupby("seqID"):
                        result[str(cid)] = ",".join(sorted(set(grp["isFamily"].astype(str))))
            except Exception:
                pass
            return result

        def parse_genomad(path):
            """Return dict contig_id -> {classification, plasmid_score, virus_score}."""
            result = {}
            try:
                df = pd.read_csv(path, sep="\t")
                # geNomad columns: seq_name, length, topology, coordinates,
                #   n_genes, genetic_code, plasmid_score, fdr, n_hallmarks,
                #   marker_enrichment, classification, taxonomy
                for _, row in df.iterrows():
                    cid = str(row.get("seq_name", row.get("sequence", "")))
                    result[cid] = {
                        "genomad_classification": row.get("classification", "chromosome"),
                        "genomad_plasmid_score": row.get("plasmid_score", np.nan),
                        "genomad_virus_score": row.get("virus_score", np.nan),
                    }
            except Exception:
                pass
            return result

        # ── load normalised abundances (sample + Contig + arg_per_genome_eq) ─
        try:
            norm_df = pd.read_csv(str(input.normalized), sep="\t")
            # Contig col in RGI has gene suffix like _1; strip it for joining
            norm_df["_contig_join"] = (
                norm_df["Contig"].astype(str).str.rsplit("_", n=1).str[0]
            )
            norm_lookup = norm_df.set_index(["Sample", "_contig_join"])[
                "arg_per_genome_equivalent"
            ].to_dict()
        except Exception:
            norm_lookup = {}

        # ── index input files by sample ──────────────────────────────────────
        def _sid_from(path, split_on, strip):
            return os.path.basename(path).replace(strip, "")

        rgi_by_sid     = {_sid_from(p, "/", "_rgi_coverage.txt"): p for p in input.rgi_cov}
        bacmet_by_sid  = {_sid_from(p, "/", "_bacmet.tsv"):        p for p in input.bacmet}
        genomad_by_sid = {
            p.split("/genomad/")[1].split("/")[0]: p for p in input.genomad
        }
        isescan_by_sid = {
            os.path.basename(p).split(".fasta")[0]: p for p in input.isescan
        }
        k2_raw_by_sid  = {_sid_from(p, "/", "_kraken.out"): p for p in input.k2_raw}
        k2_2kb_by_sid  = {_sid_from(p, "/", "_kraken.out"): p for p in input.k2_2kb}
        k2_10kb_by_sid = {_sid_from(p, "/", "_kraken.out"): p for p in input.k2_10kb}

        # ── build per-sample tables then concatenate ─────────────────────────
        all_tables = []

        for sid in list(SAMPLES.index):

            # RGI coverage — one row per ARG per contig
            rgi_path = rgi_by_sid.get(sid)
            if not rgi_path:
                continue
            try:
                rgi = pd.read_csv(rgi_path, sep="\t")
            except Exception:
                continue

            rgi["sample"] = sid
            rgi["contig_id"] = rgi["Contig"].astype(str).str.rsplit("_", n=1).str[0]

            # Count ARGs per contig
            arg_counts = rgi.groupby("contig_id")["Best_Hit_ARO"].transform("count")
            rgi["n_args_on_contig"] = arg_counts

            # Normalised abundance
            rgi["ARG_normalized_abundance"] = rgi.apply(
                lambda r: norm_lookup.get((sid, r["contig_id"]), np.nan), axis=1
            )

            # Select and rename RGI columns
            rgi_out = rgi.rename(columns={
                "Best_Hit_ARO":           "ARG_name",
                "AMR Gene Family":        "ARG_family",
                "Drug Class":             "ARG_drug_class",
                "Resistance Mechanism":   "ARG_resistance_mechanism",
                "Best_Identities":        "ARG_pct_identity",
                "Coverage":               "ARG_coverage",
            })
            keep_rgi = [
                "sample", "contig_id",
                "ARG_name", "ARG_family", "ARG_drug_class",
                "ARG_resistance_mechanism", "ARG_pct_identity",
                "ARG_coverage", "ARG_normalized_abundance",
                "n_args_on_contig"
            ]
            keep_rgi = [c for c in keep_rgi if c in rgi_out.columns]
            rgi_out = rgi_out[keep_rgi].copy()

            # ── BacMet2 — aggregate metal genes per contig ───────────────────
            bacmet_path = bacmet_by_sid.get(sid)
            metal_by_contig = {}
            if bacmet_path:
                try:
                    bm = pd.read_csv(bacmet_path, sep="\t", header=None,
                                     names=BACMET_COLS)
                    # qseqid is the Prodigal gene ID: contig_strip + _N
                    bm["contig_id"] = bm["qseqid"].astype(str).str.rsplit("_", n=1).str[0]
                    bm["metal_class"] = bm["sseqid"].apply(assign_metal_class)
                    for cid, grp in bm.groupby("contig_id"):
                        metal_by_contig[cid] = {
                            "metal_resistance_gene": ",".join(grp["sseqid"].unique()),
                            "metal_class":           ",".join(sorted(set(grp["metal_class"]))),
                            "metal_pct_identity":    round(grp["pident"].mean(), 1),
                            "n_metal_genes_on_contig": len(grp["sseqid"].unique()),
                        }
                except Exception:
                    pass

            # ── geNomad ──────────────────────────────────────────────────────
            gd_path = genomad_by_sid.get(sid)
            gd_dict = parse_genomad(gd_path) if gd_path else {}

            # ── ISEScan ──────────────────────────────────────────────────────
            is_path = isescan_by_sid.get(sid)
            is_dict = parse_isescan(is_path) if is_path else {}

            # ── Kraken2 taxonomy at three cutoffs ────────────────────────────
            k2r  = parse_kraken_out(k2_raw_by_sid.get(sid, ""))
            k2_2 = parse_kraken_out(k2_2kb_by_sid.get(sid, ""))
            k2_10= parse_kraken_out(k2_10kb_by_sid.get(sid, ""))

            # ── Join everything onto the RGI table ───────────────────────────
            def _join_col(contig_id, d, key, default=np.nan):
                return d.get(contig_id, {}).get(key, default) if isinstance(d.get(contig_id), dict) else default

            rgi_out["metal_resistance_gene"]   = rgi_out["contig_id"].map(lambda c: metal_by_contig.get(c, {}).get("metal_resistance_gene", ""))
            rgi_out["metal_class"]             = rgi_out["contig_id"].map(lambda c: metal_by_contig.get(c, {}).get("metal_class", ""))
            rgi_out["metal_pct_identity"]      = rgi_out["contig_id"].map(lambda c: metal_by_contig.get(c, {}).get("metal_pct_identity", np.nan))
            rgi_out["n_metal_genes_on_contig"] = rgi_out["contig_id"].map(lambda c: metal_by_contig.get(c, {}).get("n_metal_genes_on_contig", 0))
            rgi_out["coselected"]              = rgi_out["n_metal_genes_on_contig"] > 0

            rgi_out["genomad_classification"]  = rgi_out["contig_id"].map(lambda c: gd_dict.get(c, {}).get("genomad_classification", "chromosome"))
            rgi_out["genomad_plasmid_score"]   = rgi_out["contig_id"].map(lambda c: gd_dict.get(c, {}).get("genomad_plasmid_score", np.nan))
            rgi_out["genomad_virus_score"]     = rgi_out["contig_id"].map(lambda c: gd_dict.get(c, {}).get("genomad_virus_score", np.nan))

            rgi_out["IS_element_families"]     = rgi_out["contig_id"].map(lambda c: is_dict.get(c, ""))

            rgi_out["kraken_status_raw"]  = rgi_out["contig_id"].map(lambda c: k2r.get(c, ("U", 0))[0])
            rgi_out["kraken_taxid_raw"]   = rgi_out["contig_id"].map(lambda c: k2r.get(c, ("U", 0))[1])
            rgi_out["kraken_status_2kb"]  = rgi_out["contig_id"].map(lambda c: k2_2.get(c, ("U", 0))[0])
            rgi_out["kraken_taxid_2kb"]   = rgi_out["contig_id"].map(lambda c: k2_2.get(c, ("U", 0))[1])
            rgi_out["kraken_status_10kb"] = rgi_out["contig_id"].map(lambda c: k2_10.get(c, ("U", 0))[0])
            rgi_out["kraken_taxid_10kb"]  = rgi_out["contig_id"].map(lambda c: k2_10.get(c, ("U", 0))[1])

            all_tables.append(rgi_out)

        if all_tables:
            master = pd.concat(all_tables, ignore_index=True)
        else:
            master = pd.DataFrame()

        master.to_csv(str(output.master), sep="\t", index=False)


############################################
# Co-selection analysis: ARG + metal resistance gene on the same contig
rule coselection_analysis:
    input:
        master=os.path.join(RESULTS_DIR, "contig_summary/master_table.tsv")
    output:
        pairs=os.path.join(RESULTS_DIR, "contig_summary/coselection/coselected_pairs.tsv"),
        by_class=os.path.join(RESULTS_DIR, "contig_summary/coselection/coselection_by_classification.tsv"),
        stats=os.path.join(RESULTS_DIR, "contig_summary/coselection/coselection_statistics.tsv")
    log:
        os.path.join(RESULTS_DIR, "logs/contig_summary/coselection_analysis.log")
    message:
        "Co-selection analysis: ARG + metal resistance gene pairs on shared contigs"
    run:
        import os
        import numpy as np
        import pandas as pd
        from scipy.stats import fisher_exact, chi2_contingency

        os.makedirs(os.path.join(RESULTS_DIR, "contig_summary/coselection"), exist_ok=True)

        master = pd.read_csv(str(input.master), sep="\t")

        # ── 1. Coselected pairs table ─────────────────────────────────────────
        # One row per unique (sample, contig_id, ARG_name, metal_resistance_gene)
        co = master[master["coselected"] == True].copy()

        if co.empty:
            for o in [output.pairs, output.by_class, output.stats]:
                pd.DataFrame().to_csv(o, sep="\t", index=False)
        else:
            # Expand metal genes (may be comma-joined) into one row each
            co = co.assign(
                metal_resistance_gene=co["metal_resistance_gene"].str.split(",")
            ).explode("metal_resistance_gene")
            co["metal_resistance_gene"] = co["metal_resistance_gene"].str.strip()
            co["metal_class"] = co["metal_resistance_gene"].apply(
                lambda g: next(
                    (v for k, v in {
                        "mer": "mercury", "cop": "copper", "cus": "copper",
                        "ars": "arsenic", "cad": "cadmium", "chr": "chromate",
                        "znt": "zinc",   "znr": "zinc",    "sil": "silver",
                        "ter": "tellurite",
                    }.items() if g.lower().startswith(k)),
                    "other"
                )
            )

            # Pair-level summary: frequency + abundance + taxonomy + mobility
            pair_cols = [
                "sample", "contig_id",
                "ARG_name", "ARG_family", "ARG_drug_class",
                "ARG_resistance_mechanism", "ARG_pct_identity",
                "ARG_coverage", "ARG_normalized_abundance",
                "metal_resistance_gene", "metal_class", "metal_pct_identity",
                "genomad_classification", "genomad_plasmid_score", "genomad_virus_score",
                "IS_element_families",
                "kraken_taxid_raw", "kraken_taxid_2kb", "kraken_taxid_10kb",
                "n_args_on_contig", "n_metal_genes_on_contig"
            ]
            pair_cols = [c for c in pair_cols if c in co.columns]
            pairs_df = co[pair_cols].drop_duplicates()
            pairs_df.to_csv(str(output.pairs), sep="\t", index=False)

            # ── 2. Co-selection by contig classification ──────────────────────
            # For each genomad_classification, count unique coselected contigs
            # and compare to all ARG-bearing contigs (not coselected) to test
            # whether plasmids/viruses are enriched for co-selection.
            arg_contigs = master.drop_duplicates(subset=["sample", "contig_id"]).copy()
            arg_contigs["is_coselected"] = arg_contigs["coselected"].fillna(False)

            class_counts = (
                arg_contigs.groupby(["genomad_classification", "is_coselected"])
                .size()
                .unstack(fill_value=0)
                .reset_index()
            )
            class_counts.columns.name = None

            # Fisher's exact test: plasmid/virus vs chromosome for co-selection
            stat_rows = []
            classifications = arg_contigs["genomad_classification"].unique()
            for cls in classifications:
                a = int(arg_contigs[
                    (arg_contigs["genomad_classification"] == cls) &
                    (arg_contigs["is_coselected"])
                ].shape[0])
                b = int(arg_contigs[
                    (arg_contigs["genomad_classification"] == cls) &
                    (~arg_contigs["is_coselected"])
                ].shape[0])
                c = int(arg_contigs[
                    (arg_contigs["genomad_classification"] != cls) &
                    (arg_contigs["is_coselected"])
                ].shape[0])
                d = int(arg_contigs[
                    (arg_contigs["genomad_classification"] != cls) &
                    (~arg_contigs["is_coselected"])
                ].shape[0])
                if (a + b) > 0 and (c + d) > 0:
                    _, pval = fisher_exact([[a, b], [c, d]], alternative="greater")
                    odds_ratio = (a / b) / (c / d) if b > 0 and c > 0 else np.nan
                else:
                    pval, odds_ratio = np.nan, np.nan
                stat_rows.append({
                    "classification": cls,
                    "n_coselected": a,
                    "n_not_coselected": b,
                    "odds_ratio": round(odds_ratio, 3) if not np.isnan(odds_ratio) else np.nan,
                    "fisher_pvalue": round(pval, 6) if not np.isnan(pval) else np.nan,
                })

            by_class_df = pd.DataFrame(stat_rows).sort_values("fisher_pvalue")
            by_class_df.to_csv(str(output.by_class), sep="\t", index=False)

            # ── 3. Statistics: top ARG–metal pairs + taxonomy breakdown ───────
            # Top ARG family × metal class pairs by count
            pair_freq = (
                co.groupby(["ARG_family", "metal_class"])
                .agg(
                    n_contigs=("contig_id", "nunique"),
                    n_samples=("sample", "nunique"),
                    mean_ARG_normalized_abundance=("ARG_normalized_abundance", "mean"),
                    mean_metal_pct_identity=("metal_pct_identity", "mean"),
                    genomad_plasmid_frac=(
                        "genomad_classification",
                        lambda x: (x == "Plasmid").mean()
                    ),
                    genomad_virus_frac=(
                        "genomad_classification",
                        lambda x: (x == "Virus").mean()
                    ),
                    top_kraken_taxid_raw=(
                        "kraken_taxid_raw",
                        lambda x: x.value_counts().index[0] if len(x) > 0 else np.nan
                    ),
                )
                .reset_index()
                .sort_values("n_contigs", ascending=False)
            )

            # Add specific ARG-name level detail (top 3 ARGs per pair)
            top_args = (
                co.groupby(["ARG_family", "metal_class"])["ARG_name"]
                .apply(lambda x: ",".join(x.value_counts().head(3).index))
                .reset_index()
                .rename(columns={"ARG_name": "top_ARG_names"})
            )
            top_metals = (
                co.groupby(["ARG_family", "metal_class"])["metal_resistance_gene"]
                .apply(lambda x: ",".join(x.value_counts().head(3).index))
                .reset_index()
                .rename(columns={"metal_resistance_gene": "top_metal_genes"})
            )

            pair_freq = pair_freq.merge(top_args,  on=["ARG_family", "metal_class"], how="left")
            pair_freq = pair_freq.merge(top_metals, on=["ARG_family", "metal_class"], how="left")

            pair_freq.to_csv(str(output.stats), sep="\t", index=False)
