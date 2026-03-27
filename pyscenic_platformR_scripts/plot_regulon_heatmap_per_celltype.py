#!/usr/bin/env python
"""
Per-cell-type regulon heatmaps comparing WT vs KO.

For each cell type:
  1. Subset to cells of that type only.
  2. For each regulon, fit:  AUCell_score ~ condition  (WT=0, KO=1)
  3. Extract the t-value for the condition coefficient (beta).
  4. Select top N regulons by |t-value|.
  5. Plot per cell type:
       (a) A diverging t-value heatmap of the condition effect (top regulons).
       (b) A two-row heatmap of mean z-scored AUCell for WT vs KO across
           those regulons, for direct visual comparison.
       (c) Violin plots of per-cell AUCell score distributions, split by
           WT vs KO, for each top regulon.
       (d) DotPlot of TF gene expression (from adata.X) grouped by condition,
           showing fraction expressing and mean expression.

Inputs:
  --h5ad           : Path to AnnData (.h5ad) with cell type annotations
  --auc_csv        : Path to AUCell matrix CSV (cells x regulons), from pySCENIC
  --cell_type_col  : Column in adata.obs with cell type labels
  --condition_col  : Column in adata.obs with WT/KO labels
  --n_top          : Number of top regulons per cell type (default: 15)
  --output_dir     : Output directory for heatmaps and CSVs
"""

import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
import statsmodels.api as sm
import argparse
import os
import warnings
warnings.filterwarnings("ignore")


def regulon_to_tf(regulon_name):
    """Extract TF gene name from regulon name, e.g. 'Stat5a(+)' -> 'Stat5a'."""
    name = str(regulon_name)
    for suffix in ["(+)", "(-)", "(~)"]:
        if name.endswith(suffix):
            name = name[: -len(suffix)]
            break
    return name.strip()


def parse_args():
    parser = argparse.ArgumentParser(
        description="Per-cell-type regulon heatmaps: WT vs KO condition effect"
    )
    parser.add_argument("--h5ad", required=True,
                        help="Path to AnnData with cell type annotations")
    parser.add_argument("--auc_csv", default=None,
                        help="Path to AUCell matrix CSV (cells x regulons). "
                             "If not provided, uses adata.obsm['X_aucell'].")
    parser.add_argument("--cell_type_col", required=True,
                        help="Column in adata.obs with cell type labels")
    parser.add_argument("--condition_col", required=True,
                        help="Column in adata.obs for WT/KO condition")
    parser.add_argument("--n_top", type=int, default=15,
                        help="Number of top regulons per cell type (default: 15)")
    parser.add_argument("--output_dir", default="regulon_heatmap_per_celltype",
                        help="Output directory")
    return parser.parse_args()


def load_aucell_matrix(h5ad_path, auc_csv_path):
    """Load AnnData and AUCell matrix, return aligned (adata, auc_df)."""

    print(f"Loading AnnData from: {h5ad_path}")
    adata = sc.read_h5ad(h5ad_path)
    print(f"  Shape: {adata.shape[0]} cells x {adata.shape[1]} genes")
    print(f"  obs columns: {list(adata.obs.columns)}")

    if auc_csv_path is not None:
        print(f"Loading AUCell matrix from: {auc_csv_path}")
        auc_df = pd.read_csv(auc_csv_path, index_col=0)
        print(f"  AUCell shape: {auc_df.shape}")
    elif "X_aucell" in adata.obsm and "aucell_regulon_names" in adata.uns:
        print("Using AUCell scores from adata.obsm['X_aucell']")
        regulon_names = adata.uns["aucell_regulon_names"]
        auc_df = pd.DataFrame(
            adata.obsm["X_aucell"],
            index=adata.obs_names,
            columns=regulon_names,
        )
        print(f"  AUCell shape: {auc_df.shape}")
    else:
        raise ValueError(
            "No AUCell data found. Provide --auc_csv or ensure adata has "
            "obsm['X_aucell'] and uns['aucell_regulon_names']."
        )

    # Align cells
    common = auc_df.index.intersection(adata.obs_names)
    print(f"  Common cells between adata and AUCell: {len(common)}")
    adata = adata[common].copy()
    auc_df = auc_df.loc[common]

    return adata, auc_df


def compute_condition_tvalues_within_celltype(auc_df_ct, conditions_ct):
    """
    Given AUCell scores and condition labels for cells of ONE cell type,
    fit per-regulon:  AUCell ~ condition   (OLS with intercept).

    Returns: Series of t-values (index = regulon names), one per regulon.
             Positive t-value means higher AUCell in the second condition
             (alphabetically), e.g. WT > KO if WT sorts after KO.
    """
    unique_conds = sorted(conditions_ct.unique())
    # Encode: first alphabetically = 0, second = 1
    cond_binary = (conditions_ct == unique_conds[1]).astype(int).values

    X = np.column_stack([
        np.ones(len(cond_binary)),  # intercept
        cond_binary,                 # condition indicator
    ])

    tvals = pd.Series(index=auc_df_ct.columns, dtype=float)

    for reg in auc_df_ct.columns:
        y = auc_df_ct[reg].values
        try:
            model = sm.OLS(y, X).fit()
            tvals[reg] = model.tvalues[1]  # condition coefficient
        except Exception:
            tvals[reg] = 0.0

    return tvals


def plot_tvalue_heatmap(tvals, ct_name, n_top, cond_labels, output_path):
    """
    Plot a single-column diverging heatmap of condition t-values for the
    top N regulons (by |t-value|) within one cell type.
    """
    top_regs = tvals.abs().nlargest(n_top).index
    plot_data = tvals.loc[top_regs].to_frame(name="condition t-value")

    fig_height = max(4, n_top * 0.4 + 1.5)
    fig, ax = plt.subplots(figsize=(3.5, fig_height))

    vmax = plot_data.values.__abs__().max()
    sns.heatmap(
        plot_data,
        cmap="RdBu_r",
        center=0,
        vmin=-vmax,
        vmax=vmax,
        linewidths=0.5,
        linecolor="white",
        annot=True,
        fmt=".1f",
        annot_kws={"fontsize": 8},
        cbar_kws={"label": "t-value", "shrink": 0.6},
        ax=ax,
    )
    ax.set_ylabel("")
    ax.set_xlabel("")
    ax.tick_params(axis="y", rotation=0, labelsize=9)
    ax.tick_params(axis="x", rotation=0, labelsize=9)
    ax.set_title(
        f"{ct_name}\nCondition effect (top {n_top} regulons)\n"
        f"← {cond_labels[0]}   |   {cond_labels[1]} →",
        fontsize=11,
    )

    fig.savefig(output_path, dpi=200, bbox_inches="tight")
    plt.close()


def plot_wt_ko_mean_heatmap(auc_df_ct, conditions_ct, tvals, ct_name,
                             n_top, output_path):
    """
    Two-row heatmap: rows = conditions (WT, KO), columns = top N regulons.
    Values are z-scored mean AUCell within this cell type.
    """
    unique_conds = sorted(conditions_ct.unique())
    top_regs = tvals.abs().nlargest(n_top).index

    # Compute mean AUCell per condition for selected regulons
    means = auc_df_ct[top_regs].copy()
    means["condition"] = conditions_ct.values
    mean_per_cond = means.groupby("condition")[top_regs].mean()
    mean_per_cond = mean_per_cond.loc[unique_conds]  # enforce order

    # Z-score across the two conditions (per regulon column)
    zscore_df = (mean_per_cond - mean_per_cond.mean(axis=0)) / (
        mean_per_cond.std(axis=0) + 1e-12
    )

    # Sort regulons by t-value (most positive on left = enriched in 2nd cond)
    sorted_regs = tvals.loc[top_regs].sort_values(ascending=False).index
    zscore_df = zscore_df[sorted_regs]

    fig_width = max(8, n_top * 0.55 + 2)
    fig, ax = plt.subplots(figsize=(fig_width, 2.8))

    sns.heatmap(
        zscore_df,
        cmap="RdBu_r",
        center=0,
        linewidths=0.5,
        linecolor="white",
        xticklabels=True,
        yticklabels=True,
        cbar_kws={"label": "z-scored mean AUCell", "shrink": 0.6},
        ax=ax,
    )
    ax.set_xlabel("Regulons", fontsize=11)
    ax.set_ylabel("")
    ax.tick_params(axis="x", rotation=90, labelsize=8)
    ax.tick_params(axis="y", rotation=0, labelsize=10)
    ax.set_title(
        f"{ct_name} — Mean Regulon Activity: WT vs KO\n"
        f"(top {n_top} regulons by condition effect)",
        fontsize=11,
    )

    fig.savefig(output_path, dpi=200, bbox_inches="tight")
    plt.close()


def plot_violin_aucell(auc_df_ct, conditions_ct, tvals, ct_name, n_top,
                       output_path):
    """
    Violin plot of per-cell AUCell score distributions for top N regulons,
    split by condition (WT vs KO).  Each regulon gets its own panel.
    """
    top_regs = tvals.abs().nlargest(n_top).index.tolist()

    # Build long-form DataFrame
    plot_data = auc_df_ct[top_regs].copy()
    plot_data["condition"] = conditions_ct.values
    long_df = plot_data.melt(
        id_vars="condition",
        var_name="regulon",
        value_name="AUCell",
    )
    # Preserve order by |t-value|
    long_df["regulon"] = pd.Categorical(
        long_df["regulon"], categories=top_regs, ordered=True,
    )

    unique_conds = sorted(conditions_ct.unique())

    # Layout: up to 5 columns
    n_cols = min(5, n_top)
    n_rows = int(np.ceil(n_top / n_cols))
    fig, axes = plt.subplots(
        n_rows, n_cols,
        figsize=(n_cols * 2.8, n_rows * 3.2),
        squeeze=False,
    )

    palette = {unique_conds[0]: "#4878D0", unique_conds[1]: "#D65F5F"}

    for idx, reg in enumerate(top_regs):
        row_i, col_i = divmod(idx, n_cols)
        ax = axes[row_i][col_i]
        sub = long_df[long_df["regulon"] == reg]

        sns.violinplot(
            data=sub,
            x="condition",
            y="AUCell",
            order=unique_conds,
            palette=palette,
            inner="box",
            linewidth=0.8,
            cut=0,
            ax=ax,
        )
        t = tvals[reg]
        ax.set_title(f"{reg}\nt={t:.1f}", fontsize=9)
        ax.set_xlabel("")
        ax.set_ylabel("AUCell" if col_i == 0 else "", fontsize=8)
        ax.tick_params(labelsize=8)

    # Hide empty panels
    for idx in range(n_top, n_rows * n_cols):
        row_i, col_i = divmod(idx, n_cols)
        axes[row_i][col_i].set_visible(False)

    fig.suptitle(
        f"{ct_name} — AUCell Score Distributions (top {n_top} regulons)",
        fontsize=13, y=1.01,
    )
    fig.tight_layout()
    fig.savefig(output_path, dpi=200, bbox_inches="tight")
    plt.close()


def plot_dotplot_tf_expression(adata_ct, conditions_ct, tvals, condition_col,
                                ct_name, n_top, output_path):
    """
    Scanpy-style DotPlot of TF gene expression (from adata.X) for the top N
    regulon TFs, grouped by condition (WT vs KO).

    Dot size  = fraction of cells expressing the gene (> 0).
    Dot color = mean expression in expressing cells.
    """
    top_regs = tvals.abs().nlargest(n_top).index.tolist()
    tf_genes = [regulon_to_tf(r) for r in top_regs]

    # Map regulon -> gene, keep only genes present in adata
    reg_gene_pairs = []
    for reg, gene in zip(top_regs, tf_genes):
        if gene in adata_ct.var_names:
            reg_gene_pairs.append((reg, gene))
        else:
            # Try case-insensitive match
            matches = [g for g in adata_ct.var_names if g.lower() == gene.lower()]
            if matches:
                reg_gene_pairs.append((reg, matches[0]))

    if not reg_gene_pairs:
        print(f"  WARNING: No TF genes found in adata.var_names for {ct_name}, "
              f"skipping DotPlot.")
        return

    regulon_labels = [rg[0] for rg in reg_gene_pairs]
    gene_list = [rg[1] for rg in reg_gene_pairs]
    n_found = len(gene_list)
    n_missing = n_top - n_found
    if n_missing > 0:
        missing = set(tf_genes) - {rg[1] for rg in reg_gene_pairs}
        print(f"  DotPlot: {n_missing} TF gene(s) not in adata: {missing}")

    # Ensure condition column is set for groupby
    adata_sub = adata_ct[:, gene_list].copy()
    adata_sub.obs[condition_col] = conditions_ct.values

    # Rename var_names to regulon labels for display
    adata_sub.var_names = pd.Index(regulon_labels)

    try:
        dp = sc.pl.dotplot(
            adata_sub,
            var_names=regulon_labels,
            groupby=condition_col,
            standard_scale="var",
            show=False,
            return_fig=True,
        )
        dp.savefig(output_path, dpi=200, bbox_inches="tight")
        plt.close()
    except Exception as e:
        # Fallback: manual dotplot with seaborn
        print(f"  scanpy dotplot failed ({e}), using manual fallback")
        _manual_dotplot(adata_sub, regulon_labels, condition_col, ct_name,
                        n_top, output_path)


def _manual_dotplot(adata_sub, regulon_labels, condition_col, ct_name,
                    n_top, output_path):
    """Manual dot plot fallback if sc.pl.dotplot doesn't cooperate."""
    from scipy.sparse import issparse

    unique_conds = sorted(adata_sub.obs[condition_col].unique())
    records = []

    X = adata_sub.X
    if issparse(X):
        X = X.toarray()

    for cond in unique_conds:
        mask = (adata_sub.obs[condition_col] == cond).values
        X_sub = X[mask]
        for j, gene in enumerate(regulon_labels):
            vals = X_sub[:, j]
            frac_expr = (vals > 0).mean()
            mean_expr = vals[vals > 0].mean() if (vals > 0).any() else 0.0
            records.append({
                "condition": cond,
                "regulon": gene,
                "frac_expressing": frac_expr,
                "mean_expression": mean_expr,
            })

    df = pd.DataFrame(records)
    df["regulon"] = pd.Categorical(
        df["regulon"], categories=regulon_labels, ordered=True,
    )

    fig, ax = plt.subplots(
        figsize=(max(8, len(regulon_labels) * 0.6 + 2), 3.5)
    )

    for i, cond in enumerate(unique_conds):
        sub = df[df["condition"] == cond]
        x_pos = np.arange(len(regulon_labels)) + (i - 0.5) * 0.3
        scatter = ax.scatter(
            x_pos,
            [cond] * len(sub),
            s=sub["frac_expressing"].values * 300,
            c=sub["mean_expression"].values,
            cmap="Reds",
            edgecolors="black",
            linewidths=0.5,
            vmin=0,
        )

    ax.set_xticks(range(len(regulon_labels)))
    ax.set_xticklabels(regulon_labels, rotation=90, fontsize=8)
    ax.set_yticks(range(len(unique_conds)))
    ax.set_yticklabels(unique_conds, fontsize=10)
    ax.set_title(
        f"{ct_name} — TF Gene Expression (top {n_top} regulons)",
        fontsize=11,
    )
    plt.colorbar(scatter, ax=ax, label="Mean expression\n(expressing cells)",
                 shrink=0.6)
    # Add a size legend
    for frac in [0.25, 0.50, 0.75]:
        ax.scatter([], [], s=frac * 300, c="gray", edgecolors="black",
                   linewidths=0.5, label=f"{frac:.0%}")
    ax.legend(title="Fraction\nexpressing", loc="upper right",
              fontsize=7, title_fontsize=8, framealpha=0.8)

    fig.savefig(output_path, dpi=200, bbox_inches="tight")
    plt.close()


def plot_combined_summary(all_tvals, n_top, output_path):
    """
    Summary heatmap across ALL cell types: rows = cell types,
    columns = union of top N regulons per cell type,
    values = condition t-value within that cell type.
    """
    selected = []
    for ct, tvals in all_tvals.items():
        top_regs = tvals.abs().nlargest(n_top).index.tolist()
        selected.extend(top_regs)
    # Deduplicate preserving order
    seen = set()
    unique_selected = []
    for r in selected:
        if r not in seen:
            seen.add(r)
            unique_selected.append(r)

    # Build matrix: cell_types x regulons
    cell_types_ordered = sorted(all_tvals.keys())
    mat = pd.DataFrame(index=cell_types_ordered, columns=unique_selected, dtype=float)
    for ct in cell_types_ordered:
        for reg in unique_selected:
            mat.loc[ct, reg] = all_tvals[ct].get(reg, 0.0)

    fig_width = max(10, len(unique_selected) * 0.45 + 3)
    fig_height = max(4, len(cell_types_ordered) * 0.6 + 2)

    g = sns.clustermap(
        mat.astype(float),
        cmap="RdBu_r",
        center=0,
        figsize=(fig_width, fig_height),
        row_cluster=False,
        col_cluster=True,
        linewidths=0.5,
        linecolor="white",
        xticklabels=True,
        yticklabels=True,
        dendrogram_ratio=(0.08, 0.12),
        cbar_kws={"label": "condition t-value", "shrink": 0.5},
        cbar_pos=(0.02, 0.82, 0.03, 0.15),
    )
    g.ax_heatmap.set_xlabel("Regulons", fontsize=12)
    g.ax_heatmap.set_ylabel("")
    g.ax_heatmap.tick_params(axis="x", rotation=90, labelsize=8)
    g.ax_heatmap.tick_params(axis="y", rotation=0, labelsize=10)
    g.fig.suptitle(
        f"Within-Cell-Type Condition Effect (WT vs KO)\n"
        f"(top {n_top} regulons per cell type, union)",
        fontsize=13, y=1.02,
    )
    g.savefig(output_path, dpi=200, bbox_inches="tight")
    plt.close()
    print(f"\nSummary heatmap saved to: {output_path}")


def main():
    args = parse_args()
    os.makedirs(args.output_dir, exist_ok=True)

    # Load data
    adata, auc_df = load_aucell_matrix(args.h5ad, args.auc_csv)

    # Validate columns
    for col_name, col_label in [(args.cell_type_col, "cell_type_col"),
                                 (args.condition_col, "condition_col")]:
        if col_name not in adata.obs.columns:
            print(f"\nERROR: Column '{col_name}' not found in adata.obs.")
            print(f"Available columns: {list(adata.obs.columns)}")
            return

    cell_types = adata.obs[args.cell_type_col].astype(str)
    conditions = adata.obs[args.condition_col].astype(str)

    print(f"\nCell types in '{args.cell_type_col}':")
    print(cell_types.value_counts().to_string())
    print(f"\nConditions in '{args.condition_col}':")
    print(conditions.value_counts().to_string())

    unique_conds = sorted(conditions.unique())
    if len(unique_conds) != 2:
        print(f"\nWARNING: Expected 2 conditions, found {len(unique_conds)}: "
              f"{unique_conds}")
    print(f"  Condition encoding: {unique_conds[0]} = 0, {unique_conds[1]} = 1")

    # Drop zero-variance regulons
    nonzero_var = auc_df.var() > 0
    n_dropped = (~nonzero_var).sum()
    if n_dropped > 0:
        print(f"\nDropping {n_dropped} zero-variance regulons (global)")
        auc_df = auc_df.loc[:, nonzero_var]

    # =========================================================================
    # Per-cell-type loop
    # =========================================================================
    unique_types = sorted(cell_types.unique())
    all_tvals = {}  # ct -> Series of t-values

    for ct in unique_types:
        print(f"\n{'='*60}")
        print(f"  Cell type: {ct}")
        print(f"{'='*60}")

        # Subset to this cell type
        mask = (cell_types == ct).values
        auc_ct = auc_df.loc[mask]
        cond_ct = conditions.loc[mask]

        n_wt = (cond_ct == unique_conds[0]).sum()
        n_ko = (cond_ct == unique_conds[1]).sum()
        print(f"  {n_wt} cells in {unique_conds[0]}, "
              f"{n_ko} cells in {unique_conds[1]}")

        if n_wt < 3 or n_ko < 3:
            print(f"  SKIPPING — too few cells in one condition "
                  f"(need >= 3 per condition)")
            continue

        # Drop regulons with zero variance within this cell type
        ct_var = auc_ct.var()
        auc_ct_filt = auc_ct.loc[:, ct_var > 0]
        n_dropped_ct = auc_ct.shape[1] - auc_ct_filt.shape[1]
        if n_dropped_ct > 0:
            print(f"  Dropped {n_dropped_ct} zero-variance regulons "
                  f"within {ct}")
        print(f"  Fitting OLS: AUCell ~ condition for "
              f"{auc_ct_filt.shape[1]} regulons ...")

        # Compute condition t-values
        tvals = compute_condition_tvalues_within_celltype(auc_ct_filt, cond_ct)
        all_tvals[ct] = tvals

        # Make safe directory name
        ct_safe = ct.replace("/", "_").replace(" ", "_")
        ct_dir = os.path.join(args.output_dir, ct_safe)
        os.makedirs(ct_dir, exist_ok=True)

        # Save t-values CSV
        tval_df = tvals.to_frame(name="condition_tvalue")
        tval_df["abs_tvalue"] = tvals.abs()
        tval_df = tval_df.sort_values("abs_tvalue", ascending=False)
        tval_path = os.path.join(ct_dir, "condition_tvalues.csv")
        tval_df.to_csv(tval_path)
        print(f"  t-values saved to: {tval_path}")

        # Plot t-value heatmap
        plot_tvalue_heatmap(
            tvals, ct, args.n_top, unique_conds,
            os.path.join(ct_dir, "regulon_tvalue_heatmap.png"),
        )
        print(f"  t-value heatmap saved")

        # Plot WT vs KO mean AUCell heatmap
        plot_wt_ko_mean_heatmap(
            auc_ct_filt, cond_ct, tvals, ct, args.n_top,
            os.path.join(ct_dir, "regulon_wt_vs_ko_heatmap.png"),
        )
        print(f"  WT vs KO mean heatmap saved")

        # Violin plot of AUCell score distributions split by condition
        plot_violin_aucell(
            auc_ct_filt, cond_ct, tvals, ct, args.n_top,
            os.path.join(ct_dir, "regulon_violin_aucell.png"),
        )
        print(f"  AUCell violin plot saved")

        # DotPlot of TF gene expression split by condition
        adata_ct = adata[mask].copy()
        plot_dotplot_tf_expression(
            adata_ct, cond_ct, tvals, args.condition_col, ct, args.n_top,
            os.path.join(ct_dir, "regulon_dotplot_tf_expression.png"),
        )
        print(f"  TF expression dotplot saved")

    # =========================================================================
    # Combined summary heatmap across all cell types
    # =========================================================================
    if all_tvals:
        # Save combined t-value matrix
        combined = pd.DataFrame(all_tvals).fillna(0.0)
        combined.index.name = "regulon"
        combined_path = os.path.join(
            args.output_dir, "condition_tvalues_all_celltypes.csv"
        )
        combined.to_csv(combined_path)
        print(f"\nCombined t-value matrix saved to: {combined_path}")

        plot_combined_summary(
            all_tvals, args.n_top,
            os.path.join(args.output_dir, "summary_condition_heatmap.png"),
        )

    print("\nDone!")


if __name__ == "__main__":
    main()
