#!/usr/bin/env python
"""
pySCENIC Pipeline — DB-Gene-Filtered Version
=============================================
Identical to the original pipeline EXCEPT that the expression matrix is
filtered to only genes present in the cisTarget ranking databases BEFORE
GRNBoost2 runs.  This guarantees 100 % module–DB overlap and prevents
TFs from being pruned simply because their top co-expression targets
(Gm-predicted genes, RIKEN clones, etc.) are absent from the DB.

Phases:
  Phase 0   — Load data + load DBs + filter expression matrix to DB genes
  Phase I   — GRNBoost2 on filtered matrix
  Phase II  — RcisTarget (motif pruning)
  Phase III — AUCell (regulon activity scoring)
  Phase IV  — Differential regulon analysis (WT vs KO) + plots
"""

# =============================================================================
# NUMPY COMPATIBILITY PATCH — must come BEFORE any other imports
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
import ast
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
# TFs of Interest — for diagnostic logging
# ---------------------------------------------------------------------------
TFS_OF_INTEREST = ["Tfap2b", "Foxa1", "Esr1"]


# ===========================================================================
# Setup
# ===========================================================================

def setup_logging(output_dir):
    log_file = os.path.join(output_dir, "pyscenic_filtered_pipeline.log")
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
        description="Run DB-gene-filtered pySCENIC pipeline on CTR9 WT vs KO data"
    )
    parser.add_argument("--h5ad", required=True,
                        help="Path to combined .h5ad file")
    parser.add_argument("--condition_col", default="sample",
                        help="Column in adata.obs for WT/KO labels (default: 'sample')")
    parser.add_argument("--wt_label", default="WT",
                        help="WT label in --condition_col")
    parser.add_argument("--ko_label", default="KO",
                        help="KO label in --condition_col")
    parser.add_argument("--tf_list", required=True,
                        help="Path to TF list (allTFs_mm.txt)")
    parser.add_argument("--db_dir", required=True,
                        help="Directory containing .feather ranking databases")
    parser.add_argument("--motif_annotations", required=True,
                        help="Path to motif annotations .tbl file")
    parser.add_argument("--output_dir", default="pyscenic_filtered_output",
                        help="Output directory")
    parser.add_argument("--n_workers", type=int, default=64,
                        help="Parallel workers for GRNBoost2 / AUCell")
    parser.add_argument("--min_genes", type=int, default=200,
                        help="Minimum genes per cell (default: 200)")
    parser.add_argument(
        "--resume_from", type=str, default=None,
        choices=["modules", "regulons", "aucell"],
        help="Resume from checkpoint (loads saved intermediate files)"
    )
    parser.add_argument(
        "--tfs_of_interest", type=str, nargs="+", default=None,
        help="TF names for diagnostics (default: Tfap2b Foxa1 Esr1)"
    )
    return parser.parse_args()


# ===========================================================================
# Phase 0: Load data
# ===========================================================================

def load_data(h5ad_path, condition_col, wt_label, ko_label, min_genes, logger):
    """Load combined .h5ad and extract raw count matrix."""
    logger.info(f"Loading data from: {h5ad_path}")
    adata = sc.read_h5ad(h5ad_path)
    logger.info(f"  Loaded: {adata.shape[0]} cells x {adata.shape[1]} genes")

    # Validate condition column
    if condition_col not in adata.obs.columns:
        logger.error(f"Column '{condition_col}' not found. Available: {list(adata.obs.columns)}")
        sys.exit(1)

    conditions = adata.obs[condition_col].astype(str).unique()
    logger.info(f"  Condition column '{condition_col}' has values: {conditions}")

    if wt_label not in conditions or ko_label not in conditions:
        logger.error(f"Expected '{wt_label}' and '{ko_label}', found: {conditions}")
        sys.exit(1)

    mask = adata.obs[condition_col].isin([wt_label, ko_label])
    adata = adata[mask].copy()
    logger.info(f"  Cells after WT/KO filter: {adata.shape[0]}")

    adata.obs["condition"] = adata.obs[condition_col].astype(str).values

    # Extract raw counts
    if adata.raw is not None:
        logger.info("  Using adata.raw.X as expression matrix")
        X = adata.raw.X
        gene_names = adata.raw.var_names
    else:
        logger.info("  Using adata.X as expression matrix")
        X = adata.X
        gene_names = adata.var_names

    if sparse.issparse(X):
        logger.info("  Converting sparse → dense...")
        X = X.toarray()

    ex_matrix = pd.DataFrame(X, index=adata.obs_names, columns=gene_names)
    ex_matrix = ex_matrix.loc[:, ~ex_matrix.columns.duplicated()]
    logger.info(f"  Expression matrix: {ex_matrix.shape}")

    n_wt = (adata.obs["condition"] == wt_label).sum()
    n_ko = (adata.obs["condition"] == ko_label).sum()
    logger.info(f"  WT cells: {n_wt}, KO cells: {n_ko}")

    return adata, ex_matrix


def load_databases(db_dir, logger):
    """Load cisTarget ranking databases and return (list[DB], set of all DB genes)."""
    db_fnames = glob.glob(os.path.join(db_dir, "*.feather"))
    if not db_fnames:
        logger.error(f"No .feather databases found in {db_dir}!")
        sys.exit(1)

    dbs = [
        RankingDatabase(fname=f, name=os.path.splitext(os.path.basename(f))[0])
        for f in db_fnames
    ]
    logger.info(f"Loaded {len(dbs)} ranking databases:")

    db_gene_union = set()
    for db in dbs:
        genes = set(db.genes)
        db_gene_union |= genes
        logger.info(f"  - {db.name}: {len(genes)} genes")

    logger.info(f"  Union of all DB genes: {len(db_gene_union)}")
    return dbs, db_gene_union


def filter_expression_to_db_genes(ex_matrix, db_gene_union, tf_names, logger):
    """
    Filter expression matrix columns to only genes present in the DB union.

    TF genes are ALWAYS retained even if absent from the DB, because GRNBoost2
    needs them as predictors (they are not ranked targets in cisTarget anyway).
    """
    all_expr_genes = set(ex_matrix.columns)
    tf_set = set(tf_names) & all_expr_genes

    # Keep: (genes in DB) ∪ (TF genes in expression data)
    keep_genes = (all_expr_genes & db_gene_union) | tf_set
    removed_genes = all_expr_genes - keep_genes

    ex_filtered = ex_matrix[sorted(keep_genes)]

    logger.info("=" * 60)
    logger.info("DB-GENE FILTERING (pre-GRNBoost2)")
    logger.info("=" * 60)
    logger.info(f"  Original gene count:       {len(all_expr_genes)}")
    logger.info(f"  Genes in DB union:         {len(all_expr_genes & db_gene_union)}")
    logger.info(f"  TF genes (always kept):    {len(tf_set)}")
    logger.info(f"  Genes KEPT (DB ∪ TF):      {ex_filtered.shape[1]}")
    logger.info(f"  Genes REMOVED:             {len(removed_genes)}")

    # Log removed genes for TFs of interest
    for tf in TFS_OF_INTEREST:
        if tf in all_expr_genes and tf not in db_gene_union:
            logger.info(f"  Note: TF '{tf}' is NOT in DB but kept as predictor (TF)")
        elif tf in all_expr_genes:
            logger.info(f"  TF '{tf}': in expression data AND in DB ✓")
        else:
            logger.warning(f"  TF '{tf}': NOT in expression data!")

    # Save the list of removed genes for reference
    logger.info(f"\n  Sample of removed genes (first 50):")
    for g in sorted(removed_genes)[:50]:
        logger.info(f"    {g}")

    logger.info("=" * 60)
    return ex_filtered, removed_genes


# ===========================================================================
# Phase I: GRNBoost2
# ===========================================================================

def run_grnboost2(ex_matrix, tf_names, n_workers, output_dir, logger):
    """Run GRNBoost2 on the (filtered) expression matrix."""
    logger.info("=" * 60)
    logger.info("PHASE I: GRNBoost2 (co-expression inference)")
    logger.info(f"  Expression matrix: {ex_matrix.shape[0]} cells x {ex_matrix.shape[1]} genes")
    logger.info(f"  TFs in data: {len(set(tf_names) & set(ex_matrix.columns))}")
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

    adj_path = os.path.join(output_dir, "adjacencies.tsv")
    adjacencies.to_csv(adj_path, index=False, sep="\t")
    logger.info(f"  Saved adjacencies → {adj_path}")

    return adjacencies


def derive_modules(adjacencies, ex_matrix, output_dir, logger):
    """Derive co-expression modules from adjacencies."""
    logger.info("Deriving co-expression modules from adjacencies...")
    modules = list(modules_from_adjacencies(adjacencies, ex_matrix))
    logger.info(f"  Number of modules: {len(modules)}")

    # Quick diagnostic for TFs of interest
    for tf in TFS_OF_INTEREST:
        tf_mods = [m for m in modules if m.transcription_factor == tf]
        if tf_mods:
            for i, m in enumerate(tf_mods):
                logger.info(f"  Module for {tf} [{i+1}/{len(tf_mods)}]: "
                            f"{len(m.genes)} genes, context={m.context}")
        else:
            logger.warning(f"  No modules for {tf}")

    mod_path = os.path.join(output_dir, "modules.pkl")
    with open(mod_path, "wb") as f:
        pickle.dump(modules, f)
    logger.info(f"  Saved modules → {mod_path}")

    return modules


# ===========================================================================
# Phase II: RcisTarget
# ===========================================================================

def _fix_target_genes_col(dataframe):
    """Parse TargetGenes strings → list-of-tuples (pySCENIC serialization quirk)."""
    col_name = "TargetGenes"
    if dataframe.columns.nlevels == 2:
        matching = [c for c in dataframe.columns if c[1] == col_name]
        col_key = matching[0] if matching else None
    else:
        col_key = col_name if col_name in dataframe.columns else None

    if col_key is None:
        return dataframe

    def parse_if_str(val):
        if isinstance(val, str):
            return ast.literal_eval(val)
        return val

    dataframe[col_key] = dataframe[col_key].apply(parse_if_str)
    return dataframe


def _repair_regulons(regulons, logger):
    """Fix regulons with string gene2weight (serialization bug)."""
    from ctxcore.genesig import Regulon as Reg
    repaired = 0
    out = []
    for reg in regulons:
        g2w = reg.gene2weight
        if isinstance(g2w, str):
            try:
                parsed = ast.literal_eval(g2w)
                if isinstance(parsed, list):
                    parsed = dict(parsed)
                new_reg = Reg(
                    name=reg.name, score=reg.score, context=reg.context,
                    transcription_factor=reg.transcription_factor,
                    gene2weight=parsed, gene2occurrence=[],
                )
                out.append(new_reg)
                repaired += 1
            except Exception as e:
                logger.warning(f"  Could not repair regulon {reg.name}: {e}")
                out.append(reg)
        else:
            gene_list = list(reg.genes)
            if len(gene_list) > 5 and all(len(g) == 1 for g in gene_list[:20]):
                logger.warning(f"  Regulon {reg.name} has single-char gene names (corrupted). "
                               f"Must re-run Phase II.")
            out.append(reg)

    if repaired:
        logger.info(f"  [REPAIR] Fixed {repaired}/{len(regulons)} regulons")
    return out


def run_cistarget(modules, dbs, motif_annotations_fname, output_dir, logger):
    """Prune modules using cisTarget motif enrichment."""
    logger.info("=" * 60)
    logger.info("PHASE II: RcisTarget (motif pruning)")
    logger.info(f"  Modules: {len(modules)}")
    logger.info(f"  Databases: {len(dbs)}")
    logger.info("=" * 60)

    # Quick pre-check: DB overlap should now be ~100 %
    db_gene_union = set()
    for db in dbs:
        db_gene_union |= set(db.genes)

    for tf in TFS_OF_INTEREST:
        tf_mods = [m for m in modules if m.transcription_factor == tf]
        for m in tf_mods:
            mod_genes = set(m.genes)
            overlap = mod_genes & db_gene_union
            pct = len(overlap) / len(mod_genes) * 100 if mod_genes else 0
            logger.info(f"  Pre-cisTarget check — {tf} module ({len(mod_genes)} genes): "
                        f"{len(overlap)} in DB ({pct:.1f}%)")

    t0 = time.time()
    df = prune2df(dbs, modules, motif_annotations_fname)
    elapsed = timedelta(seconds=int(time.time() - t0))
    logger.info(f"  RcisTarget completed in {elapsed}")

    motifs_path = os.path.join(output_dir, "motifs.csv")
    df.to_csv(motifs_path)
    logger.info(f"  Saved motifs → {motifs_path}")

    # Fix TargetGenes serialization
    df = _fix_target_genes_col(df)

    # Convert to regulons
    regulons = df2regulons(df)
    logger.info(f"  Number of regulons: {len(regulons)}")

    # Check which TFs of interest survived
    regulon_tfs = set()
    for r in regulons:
        tf_name = r.transcription_factor if hasattr(r, 'transcription_factor') \
            else r.name.split("(")[0].strip()
        regulon_tfs.add(tf_name)

    logger.info("\n  === TFs of Interest — Regulon Status ===")
    for tf in TFS_OF_INTEREST:
        matching = [r for r in regulons
                    if (hasattr(r, 'transcription_factor') and r.transcription_factor == tf)
                    or r.name.startswith(tf)]
        if matching:
            for r in matching:
                logger.info(f"  ✓ {tf} → regulon '{r.name}' with {len(r.genes)} target genes")
        else:
            logger.warning(f"  ✗ {tf} → PRUNED (no regulon survived cisTarget)")

    # Debug first regulon
    if regulons:
        r0 = regulons[0]
        logger.info(f"  [DEBUG] First regulon: {r0.name}, "
                    f"type(gene2weight)={type(r0.gene2weight)}, "
                    f"first 10 genes={list(r0.genes)[:10]}")

    reg_path = os.path.join(output_dir, "regulons.pkl")
    with open(reg_path, "wb") as f:
        pickle.dump(regulons, f)
    logger.info(f"  Saved regulons → {reg_path}")

    return regulons


# ===========================================================================
# Phase III: AUCell
# ===========================================================================

def run_aucell(ex_matrix, regulons, n_workers, output_dir, logger):
    """Score regulon activity per cell using AUCell."""
    logger.info("=" * 60)
    logger.info("PHASE III: AUCell (regulon activity scoring)")
    logger.info(f"  Cells: {ex_matrix.shape[0]}, Regulons: {len(regulons)}")
    logger.info("=" * 60)

    # Gene overlap diagnostics
    ex_genes = set(ex_matrix.columns)
    logger.info(f"  {'Regulon':<30s} {'Total':>6s} {'Present':>7s} {'Missing':>7s} {'Pct':>6s}")
    for reg in regulons:
        rg = set(reg.genes)
        present = rg & ex_genes
        missing = rg - ex_genes
        pct = len(present) / len(rg) * 100 if rg else 0
        logger.info(f"  {reg.name:<30s} {len(rg):>6d} {len(present):>7d} "
                    f"{len(missing):>7d} {pct:>5.1f}%")

    t0 = time.time()
    auc_mtx = aucell(ex_matrix, regulons, num_workers=n_workers)
    elapsed = timedelta(seconds=int(time.time() - t0))
    logger.info(f"  AUCell completed in {elapsed}")

    auc_path = os.path.join(output_dir, "auc_matrix.csv")
    auc_mtx.to_csv(auc_path)
    logger.info(f"  Saved AUCell matrix → {auc_path}")

    return auc_mtx


# ===========================================================================
# Phase IV: Differential regulon activity
# ===========================================================================

def differential_regulon_activity(auc_mtx, adata, output_dir, logger,
                                  wt_label="WT", ko_label="KO"):
    """Mann-Whitney U test on regulon AUC scores between WT and KO."""
    logger.info("=" * 60)
    logger.info("PHASE IV: Differential regulon activity (WT vs KO)")
    logger.info("=" * 60)

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
    diff_df["padj"] = multipletests(diff_df["pvalue"], method="fdr_bh")[1]
    diff_df = diff_df.sort_values("padj")

    diff_path = os.path.join(output_dir, "differential_regulons_WT_vs_KO.csv")
    diff_df.to_csv(diff_path, index=False)
    logger.info(f"  Saved differential results → {diff_path}")
    logger.info(f"  Significant (padj < 0.05): {(diff_df['padj'] < 0.05).sum()}")

    top = diff_df.head(20)
    logger.info("\n  Top 20 differentially active regulons:")
    for _, row in top.iterrows():
        direction = "UP in KO" if row["log2FC_KO_vs_WT"] > 0 else "DOWN in KO"
        logger.info(f"    {row['regulon']:30s}  log2FC={row['log2FC_KO_vs_WT']:+.3f}  "
                    f"padj={row['padj']:.2e}  ({direction})")

    return diff_df


# ===========================================================================
# Plotting
# ===========================================================================

def plot_results(auc_mtx, adata, diff_df, output_dir, logger,
                 wt_label="WT", ko_label="KO"):
    """Generate summary plots (clustermap, volcano, boxplots)."""
    logger.info("Generating plots...")
    plot_dir = os.path.join(output_dir, "plots")
    os.makedirs(plot_dir, exist_ok=True)

    common_cells = auc_mtx.index.intersection(adata.obs_names)
    auc_mtx = auc_mtx.loc[common_cells]
    conditions = adata.obs.loc[common_cells, "condition"]
    color_map = {wt_label: "#1f77b4", ko_label: "#d62728"}

    # --- 1. Clustermap ---
    logger.info("  Plotting regulon clustermap...")
    top_var = auc_mtx.var().nlargest(50).index
    if len(top_var) > 0:
        plot_mtx = auc_mtx[top_var]
        nonzero_var = plot_mtx.var() > 0
        plot_mtx = plot_mtx.loc[:, nonzero_var]
        if plot_mtx.shape[1] > 0:
            if plot_mtx.shape[0] > 2000:
                idx = np.random.choice(plot_mtx.index, 2000, replace=False)
                plot_mtx = plot_mtx.loc[idx]
                cond_colors = conditions.loc[idx]
            else:
                cond_colors = conditions
            row_colors = cond_colors.map(color_map)
            try:
                g = sns.clustermap(
                    plot_mtx, row_colors=row_colors, figsize=(16, 12),
                    cmap="viridis", z_score=1, xticklabels=True, yticklabels=False,
                    dendrogram_ratio=(0.1, 0.15),
                )
                g.savefig(os.path.join(plot_dir, "regulon_clustermap.png"),
                          dpi=150, bbox_inches="tight")
                plt.close()
            except ValueError as e:
                logger.warning(f"  Clustermap failed ({e}); skipping.")
                plt.close("all")

    # --- 2. Volcano plot ---
    logger.info("  Plotting volcano plot...")
    fig, ax = plt.subplots(figsize=(10, 8))
    sig = diff_df["padj"] < 0.05
    ax.scatter(diff_df.loc[~sig, "log2FC_KO_vs_WT"],
               -np.log10(diff_df.loc[~sig, "padj"]),
               c="grey", alpha=0.5, s=30, label="Not significant")
    ax.scatter(diff_df.loc[sig, "log2FC_KO_vs_WT"],
               -np.log10(diff_df.loc[sig, "padj"]),
               c="red", alpha=0.7, s=50, label="padj < 0.05")
    for _, row in diff_df.head(15).iterrows():
        if row["padj"] < 0.05:
            ax.annotate(row["regulon"],
                        (row["log2FC_KO_vs_WT"], -np.log10(row["padj"])),
                        fontsize=7, ha="center", va="bottom")
    ax.set_xlabel("log2FC (KO / WT)")
    ax.set_ylabel("-log10(adjusted p-value)")
    ax.set_title("Differential Regulon Activity: CTR9 KO vs WT (DB-filtered)")
    ax.axhline(-np.log10(0.05), ls="--", c="black", alpha=0.3)
    ax.legend()
    fig.savefig(os.path.join(plot_dir, "volcano_regulons.png"), dpi=150, bbox_inches="tight")
    plt.close()

    # --- 3. Top regulon boxplots ---
    logger.info("  Plotting top regulon boxplots...")
    top_regulons = diff_df[diff_df["padj"] < 0.05].head(12)["regulon"].values
    if len(top_regulons) > 0:
        ncols = 4
        nrows = int(np.ceil(len(top_regulons) / ncols))
        fig, axes = plt.subplots(nrows, ncols, figsize=(4 * ncols, 4 * nrows))
        axes = axes.flatten() if len(top_regulons) > 1 else [axes]
        for i, reg in enumerate(top_regulons):
            plot_data = pd.DataFrame({
                "AUC": auc_mtx.loc[common_cells, reg].values,
                "Condition": conditions.values,
            })
            sns.boxplot(data=plot_data, x="Condition", y="AUC",
                        ax=axes[i], palette=color_map)
            padj_val = diff_df.loc[diff_df["regulon"] == reg, "padj"].values[0]
            axes[i].set_title(f"{reg}\npadj={padj_val:.2e}", fontsize=9)
        for j in range(i + 1, len(axes)):
            axes[j].set_visible(False)
        fig.suptitle("Top Differentially Active Regulons — KO vs WT (DB-filtered)", fontsize=14)
        fig.tight_layout()
        fig.savefig(os.path.join(plot_dir, "top_regulon_boxplots.png"),
                    dpi=150, bbox_inches="tight")
        plt.close()

    logger.info(f"  All plots saved → {plot_dir}")


# ===========================================================================
# Main
# ===========================================================================

def main():
    args = parse_args()
    os.makedirs(args.output_dir, exist_ok=True)
    logger = setup_logging(args.output_dir)

    logger.info("=" * 60)
    logger.info("pySCENIC Pipeline — DB-Gene-Filtered Version")
    logger.info("=" * 60)
    logger.info(f"Arguments: {vars(args)}")

    global TFS_OF_INTEREST
    if args.tfs_of_interest:
        TFS_OF_INTEREST = args.tfs_of_interest
    logger.info(f"TFs of interest: {TFS_OF_INTEREST}")

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
    # Load ranking databases (BEFORE GRNBoost2)
    # ------------------------------------------------------------------
    dbs, db_gene_union = load_databases(args.db_dir, logger)

    # ------------------------------------------------------------------
    # ★ KEY STEP: Filter expression matrix to DB genes
    # ------------------------------------------------------------------
    ex_matrix_filtered, removed_genes = filter_expression_to_db_genes(
        ex_matrix, db_gene_union, tf_names, logger
    )

    # Save the removed gene list for reference
    removed_path = os.path.join(args.output_dir, "removed_non_db_genes.txt")
    with open(removed_path, "w") as f:
        for g in sorted(removed_genes):
            f.write(g + "\n")
    logger.info(f"  Saved removed gene list → {removed_path}")

    # ------------------------------------------------------------------
    # Phase I: GRNBoost2 (on filtered matrix)
    # ------------------------------------------------------------------
    if args.resume_from in ("modules", "regulons", "aucell"):
        adj_path = os.path.join(args.output_dir, "adjacencies.tsv")
        logger.info(f"Resuming: loading adjacencies from {adj_path}")
        adjacencies = pd.read_csv(adj_path, sep="\t")
    else:
        adjacencies = run_grnboost2(
            ex_matrix_filtered, tf_names, args.n_workers, args.output_dir, logger
        )

    # Log TFs of interest in adjacencies
    for tf in TFS_OF_INTEREST:
        tf_edges = adjacencies[adjacencies["TF"] == tf]
        if len(tf_edges) > 0:
            logger.info(f"  TF '{tf}': {len(tf_edges)} target edges in GRN")
            top = tf_edges.nlargest(10, "importance")[["target", "importance"]]
            logger.info(f"  Top 10 targets:\n{top.to_string(index=False)}")
        else:
            logger.warning(f"  TF '{tf}' NOT found in adjacencies!")

    # Derive modules
    if args.resume_from in ("regulons", "aucell"):
        mod_path = os.path.join(args.output_dir, "modules.pkl")
        logger.info(f"Resuming: loading modules from {mod_path}")
        with open(mod_path, "rb") as f:
            modules = pickle.load(f)
    else:
        modules = derive_modules(adjacencies, ex_matrix_filtered, args.output_dir, logger)

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
        regulons = run_cistarget(modules, dbs, args.motif_annotations,
                                 args.output_dir, logger)

    # ------------------------------------------------------------------
    # Phase III: AUCell
    #   NOTE: We use the FILTERED expression matrix for AUCell so that
    #   the gene universe matches the regulon gene sets.
    # ------------------------------------------------------------------
    auc_mtx = run_aucell(ex_matrix_filtered, regulons, args.n_workers,
                         args.output_dir, logger)

    # ------------------------------------------------------------------
    # Phase IV: Differential analysis + plots
    # ------------------------------------------------------------------
    diff_df = differential_regulon_activity(
        auc_mtx, adata, args.output_dir, logger,
        wt_label=args.wt_label, ko_label=args.ko_label
    )
    plot_results(
        auc_mtx, adata, diff_df, args.output_dir, logger,
        wt_label=args.wt_label, ko_label=args.ko_label
    )

    # ------------------------------------------------------------------
    # Save final AnnData
    # ------------------------------------------------------------------
    logger.info("Saving final AnnData with regulon scores...")
    common = auc_mtx.index.intersection(adata.obs_names)
    adata_out = adata[common].copy()
    adata_out.obsm["X_aucell"] = auc_mtx.loc[common].values
    adata_out.uns["aucell_regulon_names"] = list(auc_mtx.columns)
    adata_out.write_h5ad(os.path.join(args.output_dir, "adata_with_aucell.h5ad"))

    logger.info("=" * 60)
    logger.info("Pipeline complete!")
    logger.info(f"All outputs saved to: {args.output_dir}")
    logger.info("=" * 60)


if __name__ == "__main__":
    main()
