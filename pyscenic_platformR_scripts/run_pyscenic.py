#!/usr/bin/env python
"""
pySCENIC Full Pipeline for CTR9 KO vs WT Analysis
==================================================
Runs the complete SCENIC pipeline on combined WT + KO single-cell data:
  Phase I   — GRNBoost2: Infer TF-target co-expression modules
  Phase II  — RcisTarget: Prune modules using cis-regulatory motif enrichment
  Phase III — AUCell: Score regulon activity per cell
  Phase IV  — Differential regulon analysis between WT and KO
"""

# =============================================================================
# NUMPY COMPATIBILITY PATCH — must come BEFORE any other imports
# pySCENIC / ctxcore uses deprecated np.object, np.bool, np.int, np.float, etc.
# These aliases were removed in NumPy 1.24. Restore them so pySCENIC works.
# =============================================================================
import numpy as np

_NP_COMPAT_ATTRS = {
    "object": object,
    "bool": bool,
    "int": int,
    "float": float,
    "complex": complex,
    "str": str,
}
for _attr, _builtin in _NP_COMPAT_ATTRS.items():
    if not hasattr(np, _attr):
        setattr(np, _attr, _builtin)

# ---------------------------------------------------------------------------

import os
import sys
import glob
import argparse
import logging
import pickle
import time
from datetime import timedelta

import pandas as pd
import scanpy as sc
import anndata as ad
from scipy import sparse
from scipy.stats import mannwhitneyu
from statsmodels.stats.multitest import multipletests

from arboreto.utils import load_tf_names
from arboreto.algo import grnboost2

from ctxcore.rnkdb import FeatherRankingDatabase as RankingDatabase
from pyscenic.utils import modules_from_adjacencies
from pyscenic.prune import prune2df, df2regulons
from pyscenic.aucell import aucell

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns

# ---------------------------------------------------------------------------
# TFs of Interest — Add your TFs here for deep cisTarget diagnostics
# ---------------------------------------------------------------------------
TFS_OF_INTEREST = ["Tfap2b", "Foxa1", "Esr1"]

# ---------------------------------------------------------------------------
# Setup
# ---------------------------------------------------------------------------

def setup_logging(output_dir):
    """Configure logging to both file and stdout."""
    log_file = os.path.join(output_dir, "pyscenic_pipeline.log")
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] %(message)s",
        handlers=[
            logging.FileHandler(log_file),
            logging.StreamHandler(sys.stdout),
        ],
    )
    return logging.getLogger(__name__)


def parse_args():
    parser = argparse.ArgumentParser(
        description="Run pySCENIC pipeline on CTR9 WT vs KO data"
    )
    parser.add_argument("--h5ad", required=True,
                        help="Path to combined .h5ad file (must contain 'sample' column in obs)")
    parser.add_argument("--condition_col", default="sample",
                        help="Column in adata.obs that distinguishes WT vs KO (default: 'sample')")
    parser.add_argument("--wt_label", default="WT",
                        help="Label for WT cells in --condition_col (default: 'WT')")
    parser.add_argument("--ko_label", default="KO",
                        help="Label for KO cells in --condition_col (default: 'KO')")
    parser.add_argument("--tf_list", required=True, help="Path to TF list (allTFs_mm.txt)")
    parser.add_argument("--db_dir", required=True,
                        help="Directory containing .feather ranking databases")
    parser.add_argument("--motif_annotations", required=True,
                        help="Path to motif annotations .tbl file")
    parser.add_argument("--output_dir", default="pyscenic_output", help="Output directory")
    parser.add_argument("--n_workers", type=int, default=64,
                        help="Number of parallel workers for GRNBoost2")
    parser.add_argument("--min_genes", type=int, default=200,
                        help="Minimum genes per cell for filtering (default: 200)")
    parser.add_argument(
        "--resume_from", type=str, default=None,
        choices=["modules", "regulons", "aucell"],
        help="Resume pipeline from a checkpoint (loads saved intermediate files)"
    )
    parser.add_argument(
        "--tfs_of_interest", type=str, nargs="+", default=None,
        help="TF names for deep cisTarget diagnostics (default: Tfap2b Foxa1 Esr1)"
    )
    return parser.parse_args()


# ---------------------------------------------------------------------------
# Phase 0: Load data from h5ad
# ---------------------------------------------------------------------------

def load_data(h5ad_path, condition_col, wt_label, ko_label, min_genes, logger):
    """Load a combined .h5ad file containing WT and KO cells and extract raw count matrix."""

    logger.info(f"Loading data from: {h5ad_path}")
    adata = sc.read_h5ad(h5ad_path)
    logger.info(f"  Loaded: {adata.shape[0]} cells x {adata.shape[1]} genes")

    # --- Validate condition column ---
    if condition_col not in adata.obs.columns:
        available = list(adata.obs.columns)
        logger.error(f"Column '{condition_col}' not found in adata.obs. Available: {available}")
        sys.exit(1)

    conditions = adata.obs[condition_col].astype(str).unique()
    logger.info(f"  Condition column '{condition_col}' has values: {conditions}")

    if wt_label not in conditions or ko_label not in conditions:
        logger.error(
            f"Expected labels '{wt_label}' and '{ko_label}' in column '{condition_col}', "
            f"but found: {conditions}"
        )
        sys.exit(1)

    # Keep only WT and KO cells (in case there are other labels)
    mask = adata.obs[condition_col].isin([wt_label, ko_label])
    adata = adata[mask].copy()
    logger.info(f"  Cells after filtering to WT/KO: {adata.shape[0]}")

    # Rename into a unified 'condition' column for downstream use
    adata.obs["condition"] = adata.obs[condition_col].astype(str).values

    # --- Extract raw count matrix ---
    # Prefer .raw if it exists and has the same cells, otherwise use .X
    if adata.raw is not None:
        logger.info("  Using adata.raw.X as expression matrix (raw counts)")
        X = adata.raw.X
        gene_names = adata.raw.var_names
    else:
        logger.info("  Using adata.X as expression matrix")
        X = adata.X
        gene_names = adata.var_names

    # Convert to dense if sparse
    if sparse.issparse(X):
        logger.info("  Converting sparse matrix to dense...")
        X = X.toarray()

    # Build DataFrame (cells × genes)
    ex_matrix = pd.DataFrame(X, index=adata.obs_names, columns=gene_names)

    # Remove duplicate gene columns
    ex_matrix = ex_matrix.loc[:, ~ex_matrix.columns.duplicated()]
    logger.info(f"  Expression matrix shape (after dedup): {ex_matrix.shape}")

    n_wt = (adata.obs["condition"] == wt_label).sum()
    n_ko = (adata.obs["condition"] == ko_label).sum()
    logger.info(f"  WT cells: {n_wt}, KO cells: {n_ko}")

    return adata, ex_matrix


# ---------------------------------------------------------------------------
# Phase I: GRNBoost2 — infer co-expression modules
# ---------------------------------------------------------------------------

def run_grnboost2(ex_matrix, tf_names, n_workers, output_dir, logger):
    """Run GRNBoost2 to infer TF-target adjacencies."""

    logger.info("=" * 60)
    logger.info("PHASE I: Running GRNBoost2 to infer co-expression modules")
    logger.info(f"  Expression matrix: {ex_matrix.shape[0]} cells x {ex_matrix.shape[1]} genes")
    logger.info(f"  TFs in expression data: {len(set(tf_names) & set(ex_matrix.columns))}")
    logger.info(f"  Workers: {n_workers}")
    logger.info("=" * 60)

    t0 = time.time()
    adjacencies = grnboost2(
        expression_data=ex_matrix,
        tf_names=tf_names,
        verbose=True,
        seed=42,
    )
    elapsed = timedelta(seconds=int(time.time() - t0))
    logger.info(f"  GRNBoost2 completed in {elapsed}")
    logger.info(f"  Adjacencies shape: {adjacencies.shape}")

    # Save checkpoint
    adj_path = os.path.join(output_dir, "adjacencies.tsv")
    adjacencies.to_csv(adj_path, index=False, sep="\t")
    logger.info(f"  Saved adjacencies to: {adj_path}")

    return adjacencies


def derive_modules(adjacencies, ex_matrix, output_dir, logger):
    """Derive potential regulon modules from adjacencies."""

    logger.info("Deriving co-expression modules from adjacencies...")
    modules = list(modules_from_adjacencies(adjacencies, ex_matrix))
    logger.info(f"  Number of modules: {len(modules)}")

    # --- Diagnostic: inspect modules for TFs of interest ---
    _log_modules_for_tfs_of_interest(modules, ex_matrix, adjacencies, output_dir, logger)

    # Save checkpoint
    mod_path = os.path.join(output_dir, "modules.pkl")
    with open(mod_path, "wb") as f:
        pickle.dump(modules, f)
    logger.info(f"  Saved modules to: {mod_path}")

    return modules


def _log_modules_for_tfs_of_interest(modules, ex_matrix, adjacencies, output_dir, logger):
    """Log detailed module info for TFs of interest before cisTarget pruning."""

    logger.info("=" * 60)
    logger.info("DIAGNOSTIC: Modules for TFs of Interest (pre-cisTarget)")
    logger.info("=" * 60)

    for tf in TFS_OF_INTEREST:
        logger.info(f"\n{'─' * 50}")
        logger.info(f"  TF: {tf}")
        logger.info(f"{'─' * 50}")

        # Check expression of the TF itself
        if tf in ex_matrix.columns:
            tf_expr = ex_matrix[tf]
            n_expressing = (tf_expr > 0).sum()
            mean_expr = tf_expr[tf_expr > 0].mean() if n_expressing > 0 else 0
            logger.info(f"  Expression: {n_expressing}/{ex_matrix.shape[0]} cells "
                         f"({n_expressing/ex_matrix.shape[0]*100:.1f}%), "
                         f"mean(nonzero)={mean_expr:.3f}")
        else:
            logger.warning(f"  *** TF '{tf}' NOT FOUND in expression matrix columns! ***")

        # Find all modules for this TF
        tf_modules = [m for m in modules if m.transcription_factor == tf]
        logger.info(f"  Modules generated: {len(tf_modules)}")

        if len(tf_modules) == 0:
            # Check adjacencies to see if there were edges
            tf_adj = adjacencies[adjacencies["TF"] == tf]
            logger.warning(f"  No modules generated for {tf}!")
            logger.info(f"  GRNBoost2 edges for {tf}: {len(tf_adj)}")
            if len(tf_adj) > 0:
                logger.info(f"  Top 10 edges by importance:")
                for _, row in tf_adj.nlargest(10, "importance").iterrows():
                    logger.info(f"    {row['target']:20s}  importance={row['importance']:.4f}")
            continue

        for i, mod in enumerate(tf_modules):
            gene_list = list(mod.genes)
            logger.info(f"\n  Module {i+1}/{len(tf_modules)}: "
                         f"name={mod.name}, genes={len(gene_list)}, "
                         f"context={mod.context}")

            # Log the module's target genes and their correlation with TF
            if tf in ex_matrix.columns and len(gene_list) > 0:
                corr_data = []
                for g in gene_list[:50]:  # first 50
                    if g in ex_matrix.columns:
                        corr_val = ex_matrix[tf].corr(ex_matrix[g])
                        g_expr_pct = (ex_matrix[g] > 0).sum() / ex_matrix.shape[0] * 100
                        corr_data.append((g, corr_val, g_expr_pct))
                corr_data.sort(key=lambda x: abs(x[1]) if not np.isnan(x[1]) else 0, reverse=True)
                logger.info(f"  Top target genes by |correlation| with {tf}:")
                logger.info(f"    {'Gene':<20s} {'Pearson r':>10s} {'% cells expr':>12s}")
                for g, r, pct in corr_data[:20]:
                    logger.info(f"    {g:<20s} {r:>10.4f} {pct:>11.1f}%")

    logger.info("=" * 60)


# ---------------------------------------------------------------------------
# Repair: fix regulons with string gene2weight (serialization bug)
# ---------------------------------------------------------------------------

def _repair_regulons(regulons, logger):
    """If gene2weight was built from a string (single-char keys), parse it back."""
    import ast
    from ctxcore.genesig import Regulon as Reg

    repaired = 0
    repaired_regulons = []
    for reg in regulons:
        g2w = reg.gene2weight
        needs_repair = False

        if isinstance(g2w, str):
            needs_repair = True
            raw = g2w
        else:
            # Detect frozendict/dict built from string iteration (single-char keys)
            gene_list = list(reg.genes)
            if len(gene_list) > 5 and all(len(g) == 1 for g in gene_list[:20]):
                needs_repair = True
                # We can't recover gene names from single chars — skip
                logger.warning(
                    f"  Regulon {reg.name} has single-char gene names (corrupted frozendict). "
                    f"Cannot repair from pickle — must re-run Phase II."
                )
                repaired_regulons.append(reg)
                continue

        if needs_repair and isinstance(g2w, str):
            try:
                parsed = ast.literal_eval(raw)
                if isinstance(parsed, list):
                    parsed = dict(parsed)
                new_reg = Reg(
                    name=reg.name,
                    score=reg.score,
                    context=reg.context,
                    transcription_factor=reg.transcription_factor,
                    gene2weight=parsed,
                    gene2occurrence=[],
                )
                repaired_regulons.append(new_reg)
                repaired += 1
            except Exception as e:
                logger.warning(f"  Could not repair regulon {reg.name}: {e}")
                repaired_regulons.append(reg)
        else:
            repaired_regulons.append(reg)

    if repaired > 0:
        logger.info(f"  [REPAIR] Fixed {repaired}/{len(regulons)} regulons with string gene2weight")
        r0 = repaired_regulons[0]
        logger.info(f"  [REPAIR] After fix — first regulon genes (first 10): {list(r0.genes)[:10]}")
        return repaired_regulons

    return regulons


# ---------------------------------------------------------------------------
# Phase II: RcisTarget — prune modules for cis-regulatory motifs
# ---------------------------------------------------------------------------

def run_cistarget(modules, dbs, motif_annotations_fname, output_dir, logger,
                  adjacencies=None, ex_matrix=None):
    """Prune modules using cisTarget motif enrichment — with deep diagnostics."""

    logger.info("=" * 60)
    logger.info("PHASE II: Running RcisTarget (motif pruning)")
    logger.info(f"  Number of modules: {len(modules)}")
    logger.info(f"  Number of databases: {len(dbs)}")
    logger.info("=" * 60)

    # =====================================================================
    # PRE-CISTARGET DIAGNOSTIC: Check module gene overlap with ranking DBs
    # =====================================================================
    _diagnose_db_gene_overlap(modules, dbs, motif_annotations_fname, logger,
                              adjacencies=adjacencies, ex_matrix=ex_matrix,
                              output_dir=output_dir)

    t0 = time.time()
    df = prune2df(dbs, modules, motif_annotations_fname)
    elapsed = timedelta(seconds=int(time.time() - t0))
    logger.info(f"  RcisTarget completed in {elapsed}")

    # Save motifs dataframe
    motifs_path = os.path.join(output_dir, "motifs.csv")
    df.to_csv(motifs_path)
    logger.info(f"  Saved motifs to: {motifs_path}")

    # =====================================================================
    # POST-CISTARGET DIAGNOSTIC: Inspect the motifs DataFrame for TFs of interest
    # =====================================================================
    _diagnose_cistarget_results(df, modules, output_dir, logger)

    # --- Debug: inspect TargetGenes column before converting to regulons ---
    df_flat = df.copy()
    if df_flat.columns.nlevels == 2:
        df_flat.columns = df_flat.columns.droplevel(0)
    if "TargetGenes" in df_flat.columns:
        sample_val = df_flat["TargetGenes"].iloc[0]
        logger.info(f"  [DEBUG] TargetGenes dtype: {type(sample_val)}")
        logger.info(f"  [DEBUG] TargetGenes sample (first 200 chars): {str(sample_val)[:200]}")

    # --- FIX: parse TargetGenes from string back to list of tuples ---
    import ast

    def _fix_target_genes_col(dataframe):
        """Parse TargetGenes strings back to list-of-tuples in-place."""
        col_name = "TargetGenes"
        if dataframe.columns.nlevels == 2:
            matching = [c for c in dataframe.columns if c[1] == col_name]
            if matching:
                col_key = matching[0]
            else:
                return dataframe
        else:
            col_key = col_name
            if col_key not in dataframe.columns:
                return dataframe

        def parse_if_str(val):
            if isinstance(val, str):
                return ast.literal_eval(val)
            return val

        dataframe[col_key] = dataframe[col_key].apply(parse_if_str)
        return dataframe

    df = _fix_target_genes_col(df)

    # Verify fix
    df_check = df.copy()
    if df_check.columns.nlevels == 2:
        df_check.columns = df_check.columns.droplevel(0)
    if "TargetGenes" in df_check.columns:
        sample_val = df_check["TargetGenes"].iloc[0]
        logger.info(f"  [DEBUG] After fix — TargetGenes dtype: {type(sample_val)}")
        if isinstance(sample_val, list) and len(sample_val) > 0:
            logger.info(f"  [DEBUG] After fix — first 3 entries: {sample_val[:3]}")

    # Convert to regulons
    regulons = df2regulons(df)
    logger.info(f"  Number of regulons: {len(regulons)}")

    # --- Post-conversion check: which TFs of interest survived? ---
    regulon_tfs = set()
    for r in regulons:
        tf_name = r.transcription_factor if hasattr(r, 'transcription_factor') else r.name.split("(")[0].strip()
        regulon_tfs.add(tf_name)

    logger.info("\n  === TFs of Interest — Final Regulon Status ===")
    for tf in TFS_OF_INTEREST:
        survived = tf in regulon_tfs
        matching_regs = [r for r in regulons
                         if (hasattr(r, 'transcription_factor') and r.transcription_factor == tf)
                         or r.name.startswith(tf)]
        if survived and matching_regs:
            for r in matching_regs:
                logger.info(f"  ✓ {tf} → regulon '{r.name}' with {len(r.genes)} target genes")
        else:
            logger.warning(f"  ✗ {tf} → PRUNED (no regulon survived cisTarget)")

    # --- Debug: inspect first regulon's gene structure ---
    if regulons:
        r0 = regulons[0]
        logger.info(f"  [DEBUG] First regulon: {r0.name}")
        logger.info(f"  [DEBUG] type(gene2weight): {type(r0.gene2weight)}")
        genes_preview = list(r0.genes)[:10]
        logger.info(f"  [DEBUG] First 10 genes: {genes_preview}")

    # Save regulons
    reg_path = os.path.join(output_dir, "regulons.pkl")
    with open(reg_path, "wb") as f:
        pickle.dump(regulons, f)
    logger.info(f"  Saved regulons to: {reg_path}")

    return regulons


def _diagnose_db_gene_overlap(modules, dbs, motif_annotations_fname, logger,
                              adjacencies=None, ex_matrix=None, output_dir=None):
    """Check whether TF-of-interest module genes exist in the ranking databases,
    cross-referenced with GRNBoost2 importance scores to quantify signal loss."""

    logger.info("\n" + "=" * 60)
    logger.info("DIAGNOSTIC: DB Gene Overlap for TFs of Interest")
    logger.info("=" * 60)

    # ------------------------------------------------------------------
    # 1. Motif annotation check — does the TF have ANY annotated motifs?
    # ------------------------------------------------------------------
    logger.info("\n--- Section 1: Motif Annotation Lookup ---")
    try:
        motif_annot = pd.read_csv(motif_annotations_fname, sep="\t", comment="#")
        logger.info(f"  Motif annotation table: {motif_annot.shape[0]} rows")
        annot_cols = motif_annot.columns.tolist()
        logger.info(f"  Annotation columns: {annot_cols}")

        for tf in TFS_OF_INTEREST:
            found_in = []
            for col in annot_cols:
                if motif_annot[col].dtype == object:
                    matches = motif_annot[motif_annot[col].astype(str).str.contains(
                        tf, case=False, na=False)]
                    if len(matches) > 0:
                        found_in.append((col, len(matches)))
            if found_in:
                logger.info(f"  TF '{tf}' found in annotation table: "
                             f"{', '.join(f'{col}({n} rows)' for col, n in found_in)}")
                for col, _ in found_in[:2]:
                    subset = motif_annot[motif_annot[col].astype(str).str.contains(
                        tf, case=False, na=False)]
                    for _, row in subset.head(5).iterrows():
                        logger.info(f"    Motif row: {dict(row)}")
            else:
                logger.warning(f"  *** TF '{tf}' NOT FOUND in motif annotation table! ***")
                logger.warning(f"      This means cisTarget cannot link any motif to {tf} "
                               f"→ all modules will be pruned regardless of enrichment.")
    except Exception as e:
        logger.warning(f"  Could not parse motif annotations for diagnostics: {e}")

    # ------------------------------------------------------------------
    # 2. Build a union of all DB genes (across all databases)
    # ------------------------------------------------------------------
    all_db_genes_per_db = {}
    db_gene_union = set()
    for db in dbs:
        try:
            genes = set(db.genes)
            all_db_genes_per_db[db.name] = genes
            db_gene_union |= genes
        except Exception as e:
            logger.warning(f"  Could not read genes from DB '{db.name}': {e}")
    logger.info(f"\n  Total unique genes across all {len(dbs)} databases: {len(db_gene_union)}")

    # ------------------------------------------------------------------
    # 3. Per-TF, per-module: importance-weighted overlap analysis
    # ------------------------------------------------------------------
    logger.info("\n--- Section 2: Importance-Weighted DB Overlap Analysis ---")

    # Collector for the detailed CSV
    all_gene_rows = []

    for tf in TFS_OF_INTEREST:
        logger.info(f"\n{'━' * 60}")
        logger.info(f"  TF: {tf}")
        logger.info(f"{'━' * 60}")

        tf_modules = [m for m in modules if m.transcription_factor == tf]
        if not tf_modules:
            logger.warning(f"  No modules for {tf} — skipping")
            continue

        # Get adjacencies for this TF (importance scores from GRNBoost2)
        tf_adj = None
        if adjacencies is not None:
            tf_adj = adjacencies[adjacencies["TF"] == tf].copy()
            if len(tf_adj) > 0:
                tf_adj = tf_adj.set_index("target")["importance"]
                logger.info(f"  GRNBoost2 edges: {len(tf_adj)}, "
                             f"total importance sum: {tf_adj.sum():.2f}")
            else:
                tf_adj = None
                logger.warning(f"  No GRNBoost2 edges found for {tf}")

        for mod_idx, mod in enumerate(tf_modules):
            mod_genes = list(mod.genes)
            logger.info(f"\n  Module {mod_idx+1}/{len(tf_modules)}: "
                         f"{mod.name}  |  {len(mod_genes)} genes  |  context={mod.context}")

            # Per-database overlap
            for db_name, db_genes in all_db_genes_per_db.items():
                overlap = set(mod_genes) & db_genes
                missing = set(mod_genes) - db_genes
                pct = len(overlap) / len(mod_genes) * 100 if mod_genes else 0
                logger.info(f"    DB '{db_name}': "
                             f"{len(overlap)}/{len(mod_genes)} in DB ({pct:.1f}%), "
                             f"{len(missing)} missing")

            # Union overlap
            overlap_union = set(mod_genes) & db_gene_union
            missing_union = set(mod_genes) - db_gene_union
            pct_union = len(overlap_union) / len(mod_genes) * 100 if mod_genes else 0
            logger.info(f"    DB union: "
                         f"{len(overlap_union)}/{len(mod_genes)} ({pct_union:.1f}%), "
                         f"{len(missing_union)} missing from ALL databases")

            # ----------------------------------------------------------
            # Importance-ranked gene table
            # ----------------------------------------------------------
            if tf_adj is not None:
                # Build per-gene table: importance, DB membership, expression stats
                gene_records = []
                for g in mod_genes:
                    importance = tf_adj.get(g, 0.0)
                    in_db = g in db_gene_union
                    # Expression stats
                    expr_pct = None
                    mean_expr = None
                    if ex_matrix is not None and g in ex_matrix.columns:
                        g_vals = ex_matrix[g]
                        expr_pct = (g_vals > 0).sum() / len(g_vals) * 100
                        mean_expr = g_vals[g_vals > 0].mean() if (g_vals > 0).any() else 0.0
                    gene_records.append({
                        "TF": tf,
                        "module": mod.name,
                        "module_idx": mod_idx + 1,
                        "context": str(mod.context),
                        "gene": g,
                        "importance": importance,
                        "in_db_union": in_db,
                        "pct_cells_expressing": expr_pct,
                        "mean_nonzero_expr": mean_expr,
                        "gene_name_pattern": _classify_gene_name(g),
                    })

                gene_df = pd.DataFrame(gene_records).sort_values("importance", ascending=False)
                all_gene_rows.extend(gene_records)

                # --- Summary statistics ---
                total_imp = gene_df["importance"].sum()
                in_db_imp = gene_df.loc[gene_df["in_db_union"], "importance"].sum()
                missing_imp = gene_df.loc[~gene_df["in_db_union"], "importance"].sum()
                pct_imp_visible = in_db_imp / total_imp * 100 if total_imp > 0 else 0
                pct_imp_lost = missing_imp / total_imp * 100 if total_imp > 0 else 0

                logger.info(f"\n    ┌─ Importance-Weighted Summary ─────────────────────┐")
                logger.info(f"    │  Total importance (all module genes):  {total_imp:>10.2f}  │")
                logger.info(f"    │  Importance VISIBLE to cisTarget:      {in_db_imp:>10.2f}  "
                             f"({pct_imp_visible:.1f}%)  │")
                logger.info(f"    │  Importance INVISIBLE (missing genes): {missing_imp:>10.2f}  "
                             f"({pct_imp_lost:.1f}%)  │")
                logger.info(f"    └────────────────────────────────────────────────────┘")

                # --- Top genes ranked by importance, annotated with DB status ---
                logger.info(f"\n    Top 30 genes by GRNBoost2 importance (▓=in DB, ░=MISSING):")
                logger.info(f"    {'Rank':<5s} {'Gene':<20s} {'Importance':>10s} {'DB?':>5s} "
                             f"{'%cells':>7s} {'MeanExpr':>9s} {'Pattern':<12s}")
                logger.info(f"    {'─'*5} {'─'*20} {'─'*10} {'─'*5} {'─'*7} {'─'*9} {'─'*12}")

                for rank, (_, row) in enumerate(gene_df.head(30).iterrows(), 1):
                    marker = "▓" if row["in_db_union"] else "░"
                    db_str = "YES" if row["in_db_union"] else "NO"
                    pct_str = f"{row['pct_cells_expressing']:.1f}" if row['pct_cells_expressing'] is not None else "?"
                    mean_str = f"{row['mean_nonzero_expr']:.3f}" if row['mean_nonzero_expr'] is not None else "?"
                    logger.info(f"    {rank:<5d} {marker} {row['gene']:<18s} {row['importance']:>10.4f} "
                                 f"{db_str:>5s} {pct_str:>7s} {mean_str:>9s} {row['gene_name_pattern']:<12s}")

                # --- Missing genes breakdown by naming pattern ---
                missing_df = gene_df[~gene_df["in_db_union"]]
                if len(missing_df) > 0:
                    pattern_counts = missing_df["gene_name_pattern"].value_counts()
                    pattern_imp = missing_df.groupby("gene_name_pattern")["importance"].sum()

                    logger.info(f"\n    Missing genes by naming pattern:")
                    logger.info(f"    {'Pattern':<20s} {'Count':>6s} {'Σ Importance':>12s} {'% of lost imp':>14s}")
                    logger.info(f"    {'─'*20} {'─'*6} {'─'*12} {'─'*14}")
                    for pat in pattern_counts.index:
                        cnt = pattern_counts[pat]
                        imp = pattern_imp.get(pat, 0)
                        pct_of_lost = imp / missing_imp * 100 if missing_imp > 0 else 0
                        logger.info(f"    {pat:<20s} {cnt:>6d} {imp:>12.4f} {pct_of_lost:>13.1f}%")

                    # --- List all missing genes sorted by importance ---
                    logger.info(f"\n    All {len(missing_df)} missing genes (sorted by importance):")
                    for _, row in missing_df.iterrows():
                        logger.info(f"      {row['gene']:<25s}  importance={row['importance']:.4f}  "
                                     f"pattern={row['gene_name_pattern']}")

                # --- Cumulative importance curve ---
                gene_df_sorted = gene_df.sort_values("importance", ascending=False).reset_index(drop=True)
                gene_df_sorted["cum_importance"] = gene_df_sorted["importance"].cumsum()
                gene_df_sorted["cum_pct"] = gene_df_sorted["cum_importance"] / total_imp * 100

                # Find where we lose the most — what rank do missing genes appear at?
                first_missing_ranks = gene_df_sorted[~gene_df_sorted["in_db_union"]].index.tolist()
                if first_missing_ranks:
                    logger.info(f"\n    Cumulative importance at key ranks:")
                    for check_rank in [5, 10, 15, 20, 25, 30, len(gene_df_sorted)]:
                        if check_rank <= len(gene_df_sorted):
                            top_n = gene_df_sorted.head(check_rank)
                            n_missing_in_top = (~top_n["in_db_union"]).sum()
                            cum_imp = top_n["importance"].sum()
                            cum_pct = cum_imp / total_imp * 100
                            missing_imp_in_top = top_n.loc[~top_n["in_db_union"], "importance"].sum()
                            logger.info(f"      Top {check_rank:>3d} genes: "
                                         f"cum_importance={cum_pct:.1f}% of total, "
                                         f"{n_missing_in_top} missing from DB "
                                         f"(lost importance={missing_imp_in_top:.4f})")

            else:
                # No adjacencies — just list missing genes
                logger.info(f"    Missing genes (no importance data available):")
                for g in sorted(missing_union):
                    logger.info(f"      {g}  pattern={_classify_gene_name(g)}")

    # ------------------------------------------------------------------
    # 4. Save detailed CSV
    # ------------------------------------------------------------------
    if all_gene_rows and output_dir:
        detail_df = pd.DataFrame(all_gene_rows).sort_values(
            ["TF", "module_idx", "importance"], ascending=[True, True, False]
        )
        detail_path = os.path.join(output_dir, "cistarget_gene_overlap_detail.csv")
        detail_df.to_csv(detail_path, index=False)
        logger.info(f"\n  Saved detailed gene-level overlap table to: {detail_path}")
        logger.info(f"  ({len(detail_df)} rows across {detail_df['TF'].nunique()} TFs, "
                     f"{detail_df['module'].nunique()} modules)")

    logger.info("=" * 60)


def _classify_gene_name(gene_name):
    """Classify a gene name into a naming pattern category."""
    import re
    if re.match(r'^Gm\d+$', gene_name):
        return "Gm_predicted"
    elif gene_name.endswith("Rik"):
        return "RIKEN_clone"
    elif re.match(r'^\d{4,}', gene_name):
        return "numeric_id"
    elif re.match(r'^[A-Z][a-z]+\d*$', gene_name):
        return "standard"
    elif re.match(r'^Mir\d+', gene_name, re.IGNORECASE):
        return "miRNA"
    elif re.match(r'^Linc', gene_name, re.IGNORECASE):
        return "lincRNA"
    elif re.match(r'^Snor', gene_name, re.IGNORECASE):
        return "snoRNA"
    elif "-" in gene_name:
        return "hyphenated"
    elif re.match(r'^[A-Z][a-z]', gene_name):
        return "standard"
    else:
        return "other"


def _diagnose_cistarget_results(df, modules, output_dir, logger):
    """Inspect the prune2df output DataFrame for TFs of interest."""

    logger.info("\n" + "=" * 60)
    logger.info("DIAGNOSTIC: cisTarget Results for TFs of Interest")
    logger.info("=" * 60)

    # Flatten multi-level columns if needed
    df_diag = df.copy()
    if df_diag.columns.nlevels == 2:
        df_diag.columns = [f"{a}_{b}" if a else b for a, b in df_diag.columns]

    logger.info(f"  Motifs DataFrame shape: {df_diag.shape}")
    logger.info(f"  Columns: {list(df_diag.columns)}")

    # Save full diagnostic CSV for TFs of interest
    diag_path = os.path.join(output_dir, "cistarget_diagnostics_tfs_of_interest.csv")

    # Try to identify the TF column
    tf_col = None
    for candidate in ["TF", "Enrichment_TF", "TF_TF", "gene_name"]:
        if candidate in df_diag.columns:
            tf_col = candidate
            break

    if tf_col is None:
        # Search for column containing TF names
        for col in df_diag.columns:
            if df_diag[col].dtype == object:
                sample_vals = df_diag[col].dropna().head(20).tolist()
                # Check if any TF of interest appears
                if any(tf in str(v) for v in sample_vals for tf in TFS_OF_INTEREST):
                    tf_col = col
                    break

    if tf_col:
        logger.info(f"  Using TF column: '{tf_col}'")
        logger.info(f"  Unique TFs in motifs df: {df_diag[tf_col].nunique()}")

        for tf in TFS_OF_INTEREST:
            tf_rows = df_diag[df_diag[tf_col].astype(str).str.contains(tf, case=False, na=False)]
            logger.info(f"\n  {'─' * 50}")
            logger.info(f"  TF: {tf} — {len(tf_rows)} rows in motifs DataFrame")
            logger.info(f"  {'─' * 50}")

            if len(tf_rows) == 0:
                logger.warning(f"  *** {tf} has ZERO rows in cisTarget output ***")
                logger.warning(f"      Possible reasons:")
                logger.warning(f"      1. No motif annotated for {tf} in the annotation table")
                logger.warning(f"      2. Module genes had insufficient overlap with DB genes")
                logger.warning(f"      3. NES score was below threshold (default: 3.0)")
                logger.warning(f"      4. AUC threshold not met (default: 0.05)")

                # Check if ANY motif was enriched for these modules (even without TF annotation)
                tf_modules = [m for m in modules if m.transcription_factor == tf]
                if tf_modules:
                    mod_names = [m.name for m in tf_modules]
                    logger.info(f"      Module names to search in motifs df: {mod_names}")
                    # Search by module name pattern
                    for mname in mod_names:
                        name_matches = df_diag[
                            df_diag.apply(lambda row: mname in str(row.values), axis=1)
                        ]
                        if len(name_matches) > 0:
                            logger.info(f"      Found {len(name_matches)} motif hits for module '{mname}'")
                continue

            # Show NES scores and other key columns
            nes_cols = [c for c in df_diag.columns if "NES" in c.upper()]
            auc_cols = [c for c in df_diag.columns if "AUC" in c.upper() and "AUCELL" not in c.upper()]
            motif_cols = [c for c in df_diag.columns if "MOTIF" in c.upper() or "MotifID" in c]
            context_cols = [c for c in df_diag.columns if "CONTEXT" in c.upper() or "Context" in c]
            target_cols = [c for c in df_diag.columns if "TARGET" in c.upper() or "TargetGenes" in c]

            display_cols = [tf_col] + nes_cols + auc_cols + motif_cols + context_cols
            display_cols = [c for c in display_cols if c in df_diag.columns]
            # Remove duplicates preserving order
            seen = set()
            display_cols = [c for c in display_cols if c not in seen and not seen.add(c)]

            logger.info(f"  Key columns for {tf}:")
            for _, row in tf_rows.iterrows():
                row_info = {c: row[c] for c in display_cols if c in row.index}
                logger.info(f"    {row_info}")

            # Show NES distribution
            for nes_col in nes_cols:
                if nes_col in tf_rows.columns:
                    nes_vals = pd.to_numeric(tf_rows[nes_col], errors="coerce").dropna()
                    if len(nes_vals) > 0:
                        logger.info(f"  NES scores ({nes_col}): "
                                     f"min={nes_vals.min():.3f}, "
                                     f"max={nes_vals.max():.3f}, "
                                     f"mean={nes_vals.mean():.3f}, "
                                     f"n_above_3.0={sum(nes_vals >= 3.0)}")

            # Show target gene counts
            for tgt_col in target_cols:
                if tgt_col in tf_rows.columns:
                    for idx, row in tf_rows.iterrows():
                        tgt_val = row[tgt_col]
                        if isinstance(tgt_val, list):
                            n_targets = len(tgt_val)
                        elif isinstance(tgt_val, str):
                            # Count tuples approximately
                            n_targets = tgt_val.count("(")
                        else:
                            n_targets = "?"
                        logger.info(f"  TargetGenes count for row {idx}: {n_targets}")

        # Save TF-of-interest rows
        tf_interest_rows = df_diag[
            df_diag[tf_col].astype(str).str.contains(
                "|".join(TFS_OF_INTEREST), case=False, na=False
            )
        ]
        if len(tf_interest_rows) > 0:
            tf_interest_rows.to_csv(diag_path, index=True)
            logger.info(f"\n  Saved TF-of-interest cisTarget rows to: {diag_path}")
    else:
        logger.warning("  Could not identify TF column in motifs DataFrame!")
        logger.info(f"  DataFrame head:\n{df_diag.head(2).to_string()}")
        # Save the full thing for manual inspection
        df_diag.head(50).to_csv(diag_path, index=True)
        logger.info(f"  Saved first 50 rows to: {diag_path}")

    # Also log a summary of ALL TFs that DID get regulons
    if tf_col:
        all_tfs = df_diag[tf_col].dropna().unique()
        logger.info(f"\n  Total unique TFs with enriched motifs: {len(all_tfs)}")
        logger.info(f"  First 30 TFs: {sorted(all_tfs)[:30]}")

    logger.info("=" * 60)


# ---------------------------------------------------------------------------
# Phase III: AUCell — score regulon activity per cell
# ---------------------------------------------------------------------------

def run_aucell(ex_matrix, regulons, n_workers, output_dir, logger):
    """Score regulon enrichment per cell using AUCell."""

    logger.info("=" * 60)
    logger.info("PHASE III: Running AUCell (regulon activity scoring)")
    logger.info(f"  Cells: {ex_matrix.shape[0]}")
    logger.info(f"  Regulons: {len(regulons)}")
    logger.info("=" * 60)

    # --- Diagnostic: regulon gene overlap with expression matrix ---
    ex_genes = set(ex_matrix.columns)
    logger.info(f"  Expression matrix has {len(ex_genes)} genes")
    logger.info("-" * 60)
    logger.info("  Regulon gene-overlap diagnostics (vs expression matrix):")
    logger.info(f"  {'Regulon':<30s} {'Total':>6s} {'Present':>7s} {'Missing':>7s} {'Pct':>6s}")
    logger.info(f"  {'-'*30} {'-'*6} {'-'*7} {'-'*7} {'-'*6}")

    for reg in regulons:
        reg_genes = set(reg.genes)
        present = reg_genes & ex_genes
        missing = reg_genes - ex_genes
        pct = (len(present) / len(reg_genes) * 100) if len(reg_genes) > 0 else 0.0

        logger.info(
            f"  {reg.name:<30s} {len(reg_genes):>6d} {len(present):>7d} "
            f"{len(missing):>7d} {pct:>5.1f}%"
        )

        # For regulons below 80% overlap, list present AND missing genes
        if pct < 80.0:
            present_sorted = sorted(present)
            missing_sorted = sorted(missing)
            logger.info(f"    ^ BELOW 80% — genes PRESENT in expr matrix ({len(present)}):")
            # Print in chunks of 10 for readability
            for i in range(0, len(present_sorted), 10):
                chunk = present_sorted[i:i+10]
                logger.info(f"      {', '.join(chunk)}")
            logger.info(f"    ^ genes MISSING from expr matrix ({len(missing)}):")
            for i in range(0, len(missing_sorted), 10):
                chunk = missing_sorted[i:i+10]
                logger.info(f"      {', '.join(chunk)}")

    logger.info("-" * 60)

    t0 = time.time()
    auc_mtx = aucell(ex_matrix, regulons, num_workers=n_workers)
    elapsed = timedelta(seconds=int(time.time() - t0))
    logger.info(f"  AUCell completed in {elapsed}")

    # Save AUCell matrix
    auc_path = os.path.join(output_dir, "auc_matrix.csv")
    auc_mtx.to_csv(auc_path)
    logger.info(f"  Saved AUCell matrix to: {auc_path}")

    return auc_mtx


# ---------------------------------------------------------------------------
# Phase IV: Differential regulon activity (WT vs KO)
# ---------------------------------------------------------------------------

def differential_regulon_activity(auc_mtx, adata, output_dir, logger,
                                   wt_label="WT", ko_label="KO"):
    """Compare regulon activity between WT and KO conditions."""

    logger.info("=" * 60)
    logger.info("PHASE IV: Differential regulon activity (WT vs KO)")
    logger.info("=" * 60)

    # Align cell barcodes
    common_cells = auc_mtx.index.intersection(adata.obs_names)
    auc_mtx = auc_mtx.loc[common_cells]
    conditions = adata.obs.loc[common_cells, "condition"]

    wt_mask = conditions == wt_label
    ko_mask = conditions == ko_label
    logger.info(f"  WT cells: {wt_mask.sum()}, KO cells: {ko_mask.sum()}")

    results = []
    for regulon in auc_mtx.columns:
        wt_scores = auc_mtx.loc[wt_mask, regulon].values
        ko_scores = auc_mtx.loc[ko_mask, regulon].values

        # Mann-Whitney U test
        try:
            stat, pval = mannwhitneyu(wt_scores, ko_scores, alternative="two-sided")
        except ValueError:
            stat, pval = np.nan, 1.0

        mean_wt = np.mean(wt_scores)
        mean_ko = np.mean(ko_scores)
        log2fc = np.log2((mean_ko + 1e-9) / (mean_wt + 1e-9))

        results.append({
            "regulon": regulon,
            "mean_AUC_WT": mean_wt,
            "mean_AUC_KO": mean_ko,
            "log2FC_KO_vs_WT": log2fc,
            "U_statistic": stat,
            "pvalue": pval,
        })

    diff_df = pd.DataFrame(results)

    # Multiple testing correction
    diff_df["padj"] = multipletests(diff_df["pvalue"], method="fdr_bh")[1]
    diff_df = diff_df.sort_values("padj")

    # Save
    diff_path = os.path.join(output_dir, "differential_regulons_WT_vs_KO.csv")
    diff_df.to_csv(diff_path, index=False)
    logger.info(f"  Saved differential regulon results to: {diff_path}")
    logger.info(f"  Significant regulons (padj < 0.05): {(diff_df['padj'] < 0.05).sum()}")

    # Log top hits
    top = diff_df.head(20)
    logger.info("\n  Top 20 differentially active regulons:")
    for _, row in top.iterrows():
        direction = "UP in KO" if row["log2FC_KO_vs_WT"] > 0 else "DOWN in KO"
        logger.info(
            f"    {row['regulon']:30s}  log2FC={row['log2FC_KO_vs_WT']:+.3f}  "
            f"padj={row['padj']:.2e}  ({direction})"
        )

    return diff_df


# ---------------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------------

def plot_results(auc_mtx, adata, diff_df, output_dir, logger,
                 wt_label="WT", ko_label="KO"):
    """Generate summary plots."""

    logger.info("Generating plots...")
    plot_dir = os.path.join(output_dir, "plots")
    os.makedirs(plot_dir, exist_ok=True)

    common_cells = auc_mtx.index.intersection(adata.obs_names)
    auc_mtx = auc_mtx.loc[common_cells]
    conditions = adata.obs.loc[common_cells, "condition"]

    color_map = {wt_label: "#1f77b4", ko_label: "#d62728"}

    # --- 1. Clustermap of top variable regulons ---
    logger.info("  Plotting regulon clustermap...")
    top_var = auc_mtx.var().nlargest(50).index
    if len(top_var) > 0:
        plot_mtx = auc_mtx[top_var]
        # Drop zero-variance columns to avoid NaN from z_score normalization
        nonzero_var = plot_mtx.var() > 0
        plot_mtx = plot_mtx.loc[:, nonzero_var]
        if plot_mtx.shape[1] == 0:
            logger.warning("  No regulons with non-zero variance; skipping clustermap.")
        else:
            # Subsample for visualization if too many cells
            if plot_mtx.shape[0] > 2000:
                idx = np.random.choice(plot_mtx.index, 2000, replace=False)
                plot_mtx = plot_mtx.loc[idx]
                cond_colors = conditions.loc[idx]
            else:
                cond_colors = conditions

            row_colors = cond_colors.map(color_map)

            try:
                g = sns.clustermap(
                    plot_mtx,
                    row_colors=row_colors,
                    figsize=(16, 12),
                    cmap="viridis",
                    z_score=1,
                    xticklabels=True,
                    yticklabels=False,
                    dendrogram_ratio=(0.1, 0.15),
                )
                g.savefig(os.path.join(plot_dir, "regulon_clustermap.png"), dpi=150, bbox_inches="tight")
                plt.close()
            except ValueError as e:
                logger.warning(f"  Clustermap failed ({e}); skipping.")
                plt.close("all")

    # --- 2. Volcano plot of differential regulons ---
    logger.info("  Plotting volcano plot...")
    fig, ax = plt.subplots(figsize=(10, 8))
    sig = diff_df["padj"] < 0.05
    ax.scatter(
        diff_df.loc[~sig, "log2FC_KO_vs_WT"],
        -np.log10(diff_df.loc[~sig, "padj"]),
        c="grey", alpha=0.5, s=30, label="Not significant",
    )
    ax.scatter(
        diff_df.loc[sig, "log2FC_KO_vs_WT"],
        -np.log10(diff_df.loc[sig, "padj"]),
        c="red", alpha=0.7, s=50, label="padj < 0.05",
    )
    # Label top hits
    for _, row in diff_df.head(15).iterrows():
        if row["padj"] < 0.05:
            ax.annotate(
                row["regulon"],
                (row["log2FC_KO_vs_WT"], -np.log10(row["padj"])),
                fontsize=7, ha="center", va="bottom",
            )
    ax.set_xlabel("log2FC (KO / WT)")
    ax.set_ylabel("-log10(adjusted p-value)")
    ax.set_title("Differential Regulon Activity: CTR9 KO vs WT")
    ax.axhline(-np.log10(0.05), ls="--", c="black", alpha=0.3)
    ax.legend()
    fig.savefig(os.path.join(plot_dir, "volcano_regulons.png"), dpi=150, bbox_inches="tight")
    plt.close()

    # --- 3. Top regulon boxplots ---
    logger.info("  Plotting top regulon boxplots...")
    top_regulons = diff_df[diff_df["padj"] < 0.05].head(12)["regulon"].values
    if len(top_regulons) > 0:
        n_plots = len(top_regulons)
        ncols = 4
        nrows = int(np.ceil(n_plots / ncols))
        fig, axes = plt.subplots(nrows, ncols, figsize=(4 * ncols, 4 * nrows))
        axes = axes.flatten() if n_plots > 1 else [axes]

        for i, reg in enumerate(top_regulons):
            plot_data = pd.DataFrame({
                "AUC": auc_mtx.loc[common_cells, reg].values,
                "Condition": conditions.values,
            })
            sns.boxplot(data=plot_data, x="Condition", y="AUC", ax=axes[i], palette=color_map)
            padj_val = diff_df.loc[diff_df["regulon"] == reg, "padj"].values[0]
            axes[i].set_title(f"{reg}\npadj={padj_val:.2e}", fontsize=9)

        # Turn off unused axes
        for j in range(i + 1, len(axes)):
            axes[j].set_visible(False)

        fig.suptitle("Top Differentially Active Regulons (KO vs WT)", fontsize=14)
        fig.tight_layout()
        fig.savefig(os.path.join(plot_dir, "top_regulon_boxplots.png"), dpi=150, bbox_inches="tight")
        plt.close()

    logger.info(f"  All plots saved to: {plot_dir}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    args = parse_args()
    os.makedirs(args.output_dir, exist_ok=True)
    logger = setup_logging(args.output_dir)

    logger.info("=" * 60)
    logger.info("pySCENIC Pipeline — CTR9 KO vs WT")
    logger.info("=" * 60)
    logger.info(f"Arguments: {vars(args)}")

    # Override TFs of interest if provided via CLI
    global TFS_OF_INTEREST
    if args.tfs_of_interest:
        TFS_OF_INTEREST = args.tfs_of_interest
    logger.info(f"TFs of interest for diagnostics: {TFS_OF_INTEREST}")

    # ------------------------------------------------------------------
    # Load data
    # ------------------------------------------------------------------
    adata, ex_matrix = load_data(
        args.h5ad, args.condition_col, args.wt_label, args.ko_label,
        args.min_genes, logger
    )

    # ------------------------------------------------------------------
    # Load TF list
    # ------------------------------------------------------------------
    logger.info(f"Loading TF names from: {args.tf_list}")
    tf_names = load_tf_names(args.tf_list)
    tf_in_data = list(set(tf_names) & set(ex_matrix.columns))
    logger.info(f"  Total TFs in list: {len(tf_names)}")
    logger.info(f"  TFs found in expression data: {len(tf_in_data)}")

    # ------------------------------------------------------------------
    # Load ranking databases
    # ------------------------------------------------------------------
    db_fnames = glob.glob(os.path.join(args.db_dir, "*.feather"))
    if not db_fnames:
        logger.error(f"No .feather databases found in {args.db_dir}!")
        sys.exit(1)

    dbs = [
        RankingDatabase(fname=f, name=os.path.splitext(os.path.basename(f))[0])
        for f in db_fnames
    ]
    logger.info(f"Loaded {len(dbs)} ranking databases:")
    for db in dbs:
        logger.info(f"  - {db.name}")

    # ------------------------------------------------------------------
    # Phase I: GRNBoost2
    # ------------------------------------------------------------------
    if args.resume_from in ("modules", "regulons", "aucell"):
        adj_path = os.path.join(args.output_dir, "adjacencies.tsv")
        logger.info(f"Resuming: loading adjacencies from {adj_path}")
        adjacencies = pd.read_csv(adj_path, sep="\t")
    else:
        adjacencies = run_grnboost2(ex_matrix, tf_names, args.n_workers, args.output_dir, logger)

    # --- Check for specific TFs of interest in the base GRN ---
    for tf in TFS_OF_INTEREST:
        tf_edges = adjacencies[adjacencies["TF"] == tf]
        if len(tf_edges) > 0:
            logger.info(f"  TF '{tf}' found in GRN adjacencies: {len(tf_edges)} target edges")
            top_targets = tf_edges.nlargest(10, "importance")[["target", "importance"]]
            logger.info(f"  Top 10 targets for {tf}:\n{top_targets.to_string(index=False)}")
        else:
            logger.warning(f"  TF '{tf}' NOT found in GRN adjacencies! "
                           f"Check if it is in TF list and expressed in data.")

    # Derive modules
    if args.resume_from in ("regulons", "aucell"):
        mod_path = os.path.join(args.output_dir, "modules.pkl")
        logger.info(f"Resuming: loading modules from {mod_path}")
        with open(mod_path, "rb") as f:
            modules = pickle.load(f)
    else:
        modules = derive_modules(adjacencies, ex_matrix, args.output_dir, logger)

    # ------------------------------------------------------------------
    # Phase II: RcisTarget
    # ------------------------------------------------------------------
    if args.resume_from == "aucell":
        reg_path = os.path.join(args.output_dir, "regulons.pkl")
        logger.info(f"Resuming: loading regulons from {reg_path}")
        with open(reg_path, "rb") as f:
            regulons = pickle.load(f)
        regulons = _repair_regulons(regulons, logger)
    else:
        regulons = run_cistarget(modules, dbs, args.motif_annotations, args.output_dir, logger,
                                adjacencies=adjacencies, ex_matrix=ex_matrix)

    # ------------------------------------------------------------------
    # Phase III: AUCell
    # ------------------------------------------------------------------
    auc_mtx = run_aucell(ex_matrix, regulons, args.n_workers, args.output_dir, logger)

    # ------------------------------------------------------------------
    # Phase IV: Differential analysis + plots
    # ------------------------------------------------------------------
    diff_df = differential_regulon_activity(auc_mtx, adata, args.output_dir, logger,
                                             wt_label=args.wt_label, ko_label=args.ko_label)
    plot_results(auc_mtx, adata, diff_df, args.output_dir, logger,
                 wt_label=args.wt_label, ko_label=args.ko_label)

    # ------------------------------------------------------------------
    # Save AUCell matrix into the AnnData object
    # ------------------------------------------------------------------
    logger.info("Saving final AnnData with regulon scores...")
    common = auc_mtx.index.intersection(adata.obs_names)
    adata_out = adata[common].copy()
    adata_out.obsm["X_aucell"] = auc_mtx.loc[common].values
    aucell_var = pd.DataFrame(index=auc_mtx.columns)
    adata_out.uns["aucell_regulon_names"] = list(auc_mtx.columns)
    adata_out.write_h5ad(os.path.join(args.output_dir, "adata_with_aucell.h5ad"))

    logger.info("=" * 60)
    logger.info("Pipeline complete!")
    logger.info(f"All outputs saved to: {args.output_dir}")
    logger.info("=" * 60)


if __name__ == "__main__":
    main()
