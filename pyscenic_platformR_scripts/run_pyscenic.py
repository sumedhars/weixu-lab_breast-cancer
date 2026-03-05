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

    # Save checkpoint
    mod_path = os.path.join(output_dir, "modules.pkl")
    with open(mod_path, "wb") as f:
        pickle.dump(modules, f)
    logger.info(f"  Saved modules to: {mod_path}")

    return modules


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

def run_cistarget(modules, dbs, motif_annotations_fname, output_dir, logger):
    """Prune modules using cisTarget motif enrichment."""

    logger.info("=" * 60)
    logger.info("PHASE II: Running RcisTarget (motif pruning)")
    logger.info(f"  Number of modules: {len(modules)}")
    logger.info(f"  Number of databases: {len(dbs)}")
    logger.info("=" * 60)

    t0 = time.time()
    df = prune2df(dbs, modules, motif_annotations_fname)
    elapsed = timedelta(seconds=int(time.time() - t0))
    logger.info(f"  RcisTarget completed in {elapsed}")

    # Save motifs dataframe
    motifs_path = os.path.join(output_dir, "motifs.csv")
    df.to_csv(motifs_path)
    logger.info(f"  Saved motifs to: {motifs_path}")

    # --- Debug: inspect TargetGenes column before converting to regulons ---
    df_flat = df.copy()
    if df_flat.columns.nlevels == 2:
        df_flat.columns = df_flat.columns.droplevel(0)
    if "TargetGenes" in df_flat.columns:
        sample_val = df_flat["TargetGenes"].iloc[0]
        logger.info(f"  [DEBUG] TargetGenes dtype: {type(sample_val)}")
        logger.info(f"  [DEBUG] TargetGenes sample (first 200 chars): {str(sample_val)[:200]}")

    # --- FIX: parse TargetGenes from string back to list of tuples ---
    # prune2df can return TargetGenes as string representations instead of
    # actual lists, which causes df2regulons to build frozendicts from
    # individual characters. Fix this by parsing string columns.
    import ast

    def _fix_target_genes_col(dataframe):
        """Parse TargetGenes strings back to list-of-tuples in-place."""
        col_name = "TargetGenes"
        # Handle multi-level columns
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
    tfs_of_interest = ["Tfap2b"]
    for tf in tfs_of_interest:
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
        regulons = run_cistarget(modules, dbs, args.motif_annotations, args.output_dir, logger)

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
