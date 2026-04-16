#!/usr/bin/env python
"""
Create a Figure 2D-style heatmap of SCENIC regulon enrichment per cell type.

Method (adapted from Brand et al. 2024 STAR Methods):
  For each regulon, fit a linear model with TWO covariates:

      AUCell_score ~ is_celltype + condition

  where is_celltype is a binary one-vs-rest indicator and condition is the
  WT/KO binary covariate. Extract the t-value for the is_celltype coefficient.
  This gives cell-type-specific regulon enrichment while controlling for
  condition (WT vs KO) effects.

  Select top N regulons per cell type by highest t-value and plot as a heatmap.

Inputs:
  --h5ad           : Path to your AnnData (.h5ad) with cell type annotations
  --auc_csv        : Path to AUCell matrix CSV (cells x regulons), from pySCENIC
  --cell_type_col  : Column in adata.obs with your 8 cell type labels
  --condition_col  : Column in adata.obs with WT/KO labels (binary covariate)
  --n_top          : Number of top regulons to select per cell type (default: 5)
  --output_dir     : Output directory for heatmaps and t-value CSV
"""

import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
import argparse
import os
import warnings
warnings.filterwarnings("ignore")


def parse_args():
    parser = argparse.ArgumentParser(
        description="Plot SCENIC regulon heatmap per cell type (Figure 2D style)"
    )
    parser.add_argument("--h5ad", required=True,
                        help="Path to AnnData with cell type annotations")
    parser.add_argument("--auc_csv", default=None,
                        help="Path to AUCell matrix CSV (cells x regulons). "
                             "If not provided, uses adata.obsm['X_aucell'].")
    parser.add_argument("--cell_type_col", required=True,
                        help="Column in adata.obs with cell type labels (8 types)")
    parser.add_argument("--condition_col", required=True,
                        help="Column in adata.obs for WT/KO condition (included as covariate)")
    parser.add_argument("--n_top", type=int, default=5,
                        help="Number of top regulons per cell type (default: 5)")
    parser.add_argument("--output_dir", default="regulon_heatmap_output",
                        help="Output directory")
    parser.add_argument("--regulons_pkl", default=None,
                        help="Path to regulons.pkl produced by run_pyscenic.py "
                             "(pyscenic_output/regulons.pkl). Used to look up the exact "
                             "number of target genes for each regulon from the Regulon "
                             "objects themselves (r.name → len(r.genes)). "
                             "If omitted, the script falls back to parsing target-gene "
                             "counts from regulon names (pattern: TF(+)(Ng)), or omits "
                             "the count entirely if that pattern is also absent.")
    return parser.parse_args()


# ---------------------------------------------------------------------------
# Helper: parse number of target genes from pySCENIC regulon names
# ---------------------------------------------------------------------------

def parse_n_targets_from_regulon_names(regulon_names):
    """
    Extract TF name and target-gene count from pySCENIC regulon name strings.

    Supported formats
    -----------------
    - ``Tfap2b(+)(52g)``  → tf_name='Tfap2b', n_targets=52
    - ``Tfap2b(-)``       → tf_name='Tfap2b', n_targets=None
    - ``Tfap2b(+)``       → tf_name='Tfap2b', n_targets=None
    - anything else       → tf_name=<original name>, n_targets=None

    Returns
    -------
    dict  :  regulon_name  →  (tf_name: str, n_targets: int | None)
    """
    import re
    # Matches e.g.  Tfap2b(+)(52g)  or  Tfap2b(-)(7g)
    _with_count    = re.compile(r'^(.+?)\([+-]\)\((\d+)g\)$')
    # Matches e.g.  Tfap2b(+)  or  Tfap2b(-)
    _without_count = re.compile(r'^(.+?)\([+-]\)$')

    result = {}
    for name in regulon_names:
        m = _with_count.match(name)
        if m:
            result[name] = (m.group(1), int(m.group(2)))
            continue
        m = _without_count.match(name)
        if m:
            result[name] = (m.group(1), None)
            continue
        result[name] = (name, None)
    return result


# ---------------------------------------------------------------------------
# Load n_targets directly from pySCENIC regulons.pkl
# ---------------------------------------------------------------------------

def load_n_targets_from_pkl(pkl_path):
    """
    Read the ``regulons.pkl`` file written by ``run_pyscenic.py`` (Phase II)
    and return a mapping of regulon name → number of target genes.

    Each entry in the pickle is a pySCENIC ``Regulon`` object with:
      - ``r.name``   — the regulon label used as AUCell column names
                       (e.g. ``'Tfap2b(+)'``)
      - ``r.genes``  — frozenset of target gene names

    Parameters
    ----------
    pkl_path : str
        Absolute or relative path to ``pyscenic_output/regulons.pkl``.

    Returns
    -------
    dict
        ``{ regulon_name: n_targets }``  e.g. ``{'Tfap2b(+)': 52, ...}``
        Returns an empty dict if the file cannot be loaded.
    """
    import pickle

    if not os.path.exists(pkl_path):
        print(f"  WARNING: --regulons_pkl path not found: {pkl_path}")
        return {}

    print(f"Loading regulon target-gene counts from: {pkl_path}")
    try:
        with open(pkl_path, "rb") as fh:
            regulons = pickle.load(fh)
    except Exception as e:
        print(f"  WARNING: Could not load regulons.pkl ({e}). "
              f"Falling back to name-parsing.")
        return {}

    n_targets_map = {}
    for r in regulons:
        try:
            genes = r.genes                   # frozenset
            n     = len(genes)
            name  = r.name                    # e.g. 'Tfap2b(+)'
            n_targets_map[name] = n
        except Exception:
            pass                              # malformed entry — skip

    print(f"  Loaded n_targets for {len(n_targets_map)} regulons from pkl")
    if n_targets_map:
        sample = list(n_targets_map.items())[:4]
        print(f"  Sample mapping: { {k: v for k, v in sample} }")

    return n_targets_map


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


def compute_tvalues_per_celltype(auc_df, cell_types, conditions):
    """
    For each regulon, fit a linear model:

        AUCell_score ~ is_celltype + condition

    where:
      - is_celltype = 1 if cell belongs to the focal cluster, 0 otherwise
      - condition   = 1 for KO, 0 for WT (binary covariate)

    Extract the t-value for the is_celltype coefficient. This gives the
    cell-type-specific enrichment of the regulon *controlling for* the
    WT/KO condition effect.

    Returns: DataFrame (regulons x cell_types) of t-values.
    """
    import statsmodels.api as sm

    unique_types = sorted(cell_types.unique())
    unique_conds = sorted(conditions.unique())
    print(f"\nComputing t-values for {len(unique_types)} cell types "
          f"across {auc_df.shape[1]} regulons...")
    print(f"  Conditions in model: {unique_conds}")

    # Encode condition as binary (alphabetical: first = 0, second = 1)
    # e.g. KO_DM = 0, WT_DM = 1  or  KO = 0, WT = 1
    cond_binary = (conditions == unique_conds[1]).astype(int).values
    print(f"  Condition encoding: {unique_conds[0]} = 0, {unique_conds[1]} = 1")

    tvalue_matrix = pd.DataFrame(
        index=auc_df.columns,  # regulons
        columns=unique_types,
        dtype=float,
    )

    n_regulons = auc_df.shape[1]

    for ct in unique_types:
        # Binary indicator: 1 if cell is this type, 0 otherwise
        is_ct = (cell_types == ct).astype(int).values
        n_in = is_ct.sum()
        n_out = len(is_ct) - n_in

        # Design matrix: [intercept, is_celltype, condition]
        X = np.column_stack([
            np.ones(len(is_ct)),   # intercept
            is_ct,                  # cell type indicator
            cond_binary,            # condition (WT/KO)
        ])

        for i, reg in enumerate(auc_df.columns):
            y = auc_df[reg].values

            try:
                model = sm.OLS(y, X).fit()
                # Coefficient index 1 = is_celltype
                tvalue_matrix.loc[reg, ct] = model.tvalues[1]
            except Exception:
                tvalue_matrix.loc[reg, ct] = 0.0

        print(f"  {ct}: {n_in} cells (vs {n_out} rest) — done")

    tvalue_matrix = tvalue_matrix.astype(float)
    return tvalue_matrix


def select_top_regulons(tvalue_matrix, n_top=5):
    """
    Select top N regulons per cell type by highest (most positive) t-value,
    keeping the union across all cell types (no duplicates).
    """

    selected = []
    print(f"\n--- Top {n_top} regulons per cell type (by highest t-value) ---")
    for ct in tvalue_matrix.columns:
        top_series = tvalue_matrix[ct].nlargest(n_top)
        top_idx = top_series.index.tolist()
        selected.extend(top_idx)
        print(f"\n  {ct}:")
        for rank, reg in enumerate(top_idx, 1):
            tval = tvalue_matrix.loc[reg, ct]
            print(f"    {rank:2d}. {reg:30s}  t = {tval:+.4f}")

    # Remove duplicates while preserving order
    seen = set()
    unique_selected = []
    for r in selected:
        if r not in seen:
            seen.add(r)
            unique_selected.append(r)

    print(f"\nSelected {len(unique_selected)} unique regulons "
          f"({n_top} per cell type, union)")
    return unique_selected


def plot_heatmap_sorted(tvalue_matrix, selected_regulons, output_path, n_top,
                        regulons_pkl=None):
    """
    Publication-quality sorted heatmap:

    Layout
    ------
    - Y-axis (rows)  : TF / regulon, labelled as ``TF_name (N_targets)``
    - X-axis (cols)  : cell types, with a hierarchical dendrogram on top
    - Row ordering   : regulons grouped by the cell type in which they score
                       highest; within each group sorted by descending t-value
    - Separators     : horizontal dotted lines between each cell-type group

    Target-gene counts
    ------------------
    Priority order:
      1. ``regulons.pkl`` from run_pyscenic.py — authoritative source; reads
         ``len(r.genes)`` directly from the pySCENIC Regulon objects.
      2. Parsed from the regulon name itself  (pattern ``TF(+)(Ng)``).
      3. TF name only (no count) if neither source yields a number.
    """
    from scipy.cluster.hierarchy import linkage
    from scipy.spatial.distance import pdist

    # ------------------------------------------------------------------
    # 1. Sort regulons by their top cell type, track group boundaries
    # ------------------------------------------------------------------
    max_ct     = tvalue_matrix.loc[selected_regulons].idxmax(axis=1)
    ct_order   = list(tvalue_matrix.columns)          # original cell-type order

    sorted_regulons = []
    block_sizes     = []
    for ct in ct_order:
        regs = max_ct[max_ct == ct].index.tolist()
        regs = sorted(regs,
                      key=lambda r: tvalue_matrix.loc[r, ct],
                      reverse=True)
        sorted_regulons.extend(regs)
        block_sizes.append(len(regs))

    # ------------------------------------------------------------------
    # 2. Build plot matrix: regulons × cell_types
    # ------------------------------------------------------------------
    plot_df = tvalue_matrix.loc[sorted_regulons, ct_order].copy()

    # ------------------------------------------------------------------
    # 3. Build Y-axis labels:  TF_name (n_targets)
    #    Priority 1: regulons.pkl  →  len(r.genes)
    #    Priority 2: parse count from regulon name  (pattern TF(+)(Ng))
    #    Priority 3: TF name only
    # ------------------------------------------------------------------
    pkl_map = load_n_targets_from_pkl(regulons_pkl) if regulons_pkl else {}

    # Name-parsing gives us clean TF names even when the count is absent
    parsed = parse_n_targets_from_regulon_names(sorted_regulons)

    def make_label(reg):
        tf_name, n_parsed = parsed.get(reg, (reg, None))
        # pkl is authoritative; fall back to whatever the name encodes
        n = pkl_map.get(reg, n_parsed)
        return f"{tf_name} ({n})" if n is not None else tf_name

    y_labels = [make_label(r) for r in sorted_regulons]
    plot_df.index = y_labels

    # ------------------------------------------------------------------
    # 4. Cluster cell types (columns) for the X-axis dendrogram
    # ------------------------------------------------------------------
    col_linkage = linkage(
        pdist(plot_df.T.values, metric="euclidean"),
        method="average",
    )

    # ------------------------------------------------------------------
    # 5. Build figure with clustermap
    # ------------------------------------------------------------------
    n_rows = len(sorted_regulons)
    n_cols = len(ct_order)

    fig_height = max(6, n_rows * 0.32 + 2.5)
    fig_width  = max(5, n_cols * 1.1  + 3.0)

    g = sns.clustermap(
        plot_df,
        cmap="RdBu_r",
        center=0,
        figsize=(fig_width, fig_height),
        row_cluster=False,          # fixed regulon order
        col_cluster=True,           # cluster cell types → dendrogram on top
        col_linkage=col_linkage,
        z_score=None,
        linewidths=0.4,
        linecolor="white",
        xticklabels=True,
        yticklabels=True,
        dendrogram_ratio=0.12,      # space for the column (cell-type) dendrogram
        cbar_kws={"label": "t-value (one-vs-rest)", "shrink": 0.55,
                  "orientation": "vertical"},
    )

    # ------------------------------------------------------------------
    # 6. Horizontal dotted separator lines between cell-type blocks
    # ------------------------------------------------------------------
    cumulative  = 0
    boundaries  = []
    for size in block_sizes:
        cumulative += size
        if size > 0:
            boundaries.append(cumulative)

    # All boundaries except the very last (end of heatmap)
    for b in boundaries[:-1]:
        g.ax_heatmap.axhline(
            y=b, color="black", linewidth=1.1,
            linestyle="--", alpha=0.75,
        )

    # ------------------------------------------------------------------
    # 7. Axis labels & tick formatting
    # ------------------------------------------------------------------
    g.ax_heatmap.set_xlabel("Cell Types", fontsize=12, labelpad=8)
    g.ax_heatmap.set_ylabel("TF / Regulon  (# target genes)", fontsize=11, labelpad=8)
    g.ax_heatmap.tick_params(axis="x", rotation=45, labelsize=9)
    g.ax_heatmap.tick_params(axis="y", rotation=0,  labelsize=8)

    g.fig.suptitle(
        f"SCENIC Regulon Enrichment per Cell Type\n"
        f"(top {n_top} regulons per type · one-vs-rest t-values)",
        fontsize=13, y=1.01,
    )

    g.savefig(output_path, dpi=200, bbox_inches="tight")
    plt.close()
    print(f"Sorted heatmap saved to: {output_path}")


def compute_tvalues_per_group(auc_df, cell_types, conditions):
    """
    For each regulon and each (cell_type, condition) group, fit:

        AUCell_score ~ is_group

    where is_group = 1 if the cell belongs to e.g. Fibroblast_KO, 0 otherwise.
    No separate condition covariate — condition is already baked into the group.

    Returns: DataFrame (regulons x groups) of t-values, where group labels
             are 'celltype_condition' (e.g. 'Fibroblast_KO_DM').
    """
    import statsmodels.api as sm

    unique_types = sorted(cell_types.unique())
    unique_conds = sorted(conditions.unique())

    # Build group labels
    group_labels = cell_types.astype(str) + "_" + conditions.astype(str)

    # Ordered: cell types first, conditions nested within
    group_order = []
    for ct in unique_types:
        for cond in unique_conds:
            group_order.append(f"{ct}_{cond}")

    print(f"\nComputing t-values for {len(group_order)} "
          f"(cell_type × condition) groups "
          f"across {auc_df.shape[1]} regulons...")

    tvalue_matrix = pd.DataFrame(
        index=auc_df.columns,  # regulons
        columns=group_order,
        dtype=float,
    )

    for grp in group_order:
        is_grp = (group_labels == grp).astype(int).values
        n_in = is_grp.sum()
        n_out = len(is_grp) - n_in

        if n_in == 0:
            print(f"  {grp}: 0 cells — skipping")
            tvalue_matrix[grp] = 0.0
            continue

        # Design matrix: [intercept, is_group]
        X = np.column_stack([
            np.ones(len(is_grp)),  # intercept
            is_grp,                 # group indicator
        ])

        for i, reg in enumerate(auc_df.columns):
            y = auc_df[reg].values
            try:
                model = sm.OLS(y, X).fit()
                # Coefficient index 1 = is_group
                tvalue_matrix.loc[reg, grp] = model.tvalues[1]
            except Exception:
                tvalue_matrix.loc[reg, grp] = 0.0

        print(f"  {grp}: {n_in} cells (vs {n_out} rest) — done")

    tvalue_matrix = tvalue_matrix.astype(float)
    return tvalue_matrix


def plot_heatmap_by_condition(
    tvalue_group_matrix,
    selected_regulons,
    tvalue_celltype_matrix,
    cell_types,
    conditions,
    output_path,
    n_top,
):
    """
    Plot heatmap with y-axis = celltype_condition and x-axis = regulons.
    Values are one-vs-rest t-values from the per-group model.
    Regulon ordering reuses the cell-type-level t-value selection.
    """

    plot_df = tvalue_group_matrix.loc[selected_regulons].T  # groups x regulons

    # Sort regulons by which cell type they are most enriched in
    # (using the original cell-type-level t-values)
    max_ct = tvalue_celltype_matrix.loc[selected_regulons].idxmax(axis=1)
    unique_types = sorted(cell_types.unique())
    unique_conds = sorted(conditions.unique())
    sorted_regulons = []
    for ct in unique_types:
        regs_for_ct = max_ct[max_ct == ct].index.tolist()
        regs_for_ct = sorted(
            regs_for_ct,
            key=lambda r: tvalue_celltype_matrix.loc[r, ct],
            reverse=True,
        )
        sorted_regulons.extend(regs_for_ct)

    # Order rows: group by cell type, conditions side-by-side
    row_order = []
    for ct in unique_types:
        for cond in unique_conds:
            label = f"{ct}_{cond}"
            if label in plot_df.index:
                row_order.append(label)

    plot_df = plot_df.loc[row_order, sorted_regulons]

    # --- Plot ---
    fig_width = max(10, len(sorted_regulons) * 0.45 + 3)
    fig_height = max(5, len(row_order) * 0.45 + 2)

    fig, ax = plt.subplots(figsize=(fig_width, fig_height))
    sns.heatmap(
        plot_df,
        cmap="RdBu_r",
        center=0,
        linewidths=0.5,
        linecolor="white",
        xticklabels=True,
        yticklabels=True,
        cbar_kws={"label": "t-value", "shrink": 0.5},
        ax=ax,
    )
    ax.set_xlabel("Regulons", fontsize=12)
    ax.set_ylabel("")
    ax.tick_params(axis="x", rotation=90, labelsize=8)
    ax.tick_params(axis="y", rotation=0, labelsize=9)
    ax.set_title(
        f"SCENIC Regulon Enrichment by Cell Type × Condition\n"
        f"(top {n_top} regulons per type, one-vs-rest t-values per group)",
        fontsize=13,
    )

    # Add horizontal lines to separate cell type blocks
    n_conds = len(unique_conds)
    for i in range(1, len(unique_types)):
        ax.axhline(y=i * n_conds, color="black", linewidth=1.5)

    fig.savefig(output_path, dpi=200, bbox_inches="tight")
    plt.close()
    print(f"Condition-split heatmap saved to: {output_path}")





def main():
    args = parse_args()
    os.makedirs(args.output_dir, exist_ok=True)

    # Load data
    adata, auc_df = load_aucell_matrix(args.h5ad, args.auc_csv)

    # Validate cell type column
    if args.cell_type_col not in adata.obs.columns:
        print(f"\nERROR: Column '{args.cell_type_col}' not found in adata.obs.")
        print(f"Available columns: {list(adata.obs.columns)}")
        return

    # Validate condition column
    if args.condition_col not in adata.obs.columns:
        print(f"\nERROR: Column '{args.condition_col}' not found in adata.obs.")
        print(f"Available columns: {list(adata.obs.columns)}")
        return

    cell_types = adata.obs[args.cell_type_col].astype(str)
    conditions = adata.obs[args.condition_col].astype(str)

    print(f"\nCell types in '{args.cell_type_col}':")
    print(cell_types.value_counts().to_string())
    print(f"\nConditions in '{args.condition_col}':")
    print(conditions.value_counts().to_string())

    if conditions.nunique() != 2:
        print(f"\nWARNING: Expected exactly 2 conditions for the binary covariate, "
              f"found {conditions.nunique()}: {sorted(conditions.unique())}")
        print("The model will encode these alphabetically as 0/1.")

    # Drop regulons with zero variance (uninformative)
    nonzero_var = auc_df.var() > 0
    n_dropped = (~nonzero_var).sum()
    if n_dropped > 0:
        print(f"\nDropping {n_dropped} zero-variance regulons")
        auc_df = auc_df.loc[:, nonzero_var]

    # Compute t-values with linear model: AUCell ~ is_celltype + condition
    tvalue_matrix = compute_tvalues_per_celltype(auc_df, cell_types, conditions)

    # Save full t-value matrix
    tval_path = os.path.join(args.output_dir, "regulon_tvalues_per_celltype.csv")
    tvalue_matrix.to_csv(tval_path)
    print(f"Full t-value matrix saved to: {tval_path}")

    # Select top regulons
    selected = select_top_regulons(tvalue_matrix, n_top=args.n_top)

    # Plot sorted heatmap:
    #   - regulons on Y-axis, cell types on X-axis
    #   - dendrogram clusters cell types on X
    #   - dotted lines separate regulon groups per cell type
    #   - TF labels include number of target genes
    plot_heatmap_sorted(
        tvalue_matrix, selected,
        os.path.join(args.output_dir, "regulon_heatmap_sorted.png"),
        args.n_top,
        regulons_pkl=args.regulons_pkl,
    )

    # ---- Condition-split heatmap (celltype_WT / celltype_KO rows) ----
    # Fit one-vs-rest model per (celltype, condition) group:
    #   AUCell ~ is_group  (no separate condition covariate)
    tvalue_group_matrix = compute_tvalues_per_group(auc_df, cell_types, conditions)

    # Save full group-level t-value matrix
    tval_grp_path = os.path.join(args.output_dir, "regulon_tvalues_per_group.csv")
    tvalue_group_matrix.to_csv(tval_grp_path)
    print(f"Group-level t-value matrix saved to: {tval_grp_path}")

    # Sorted condition-split heatmap
    plot_heatmap_by_condition(
        tvalue_group_matrix, selected, tvalue_matrix,
        cell_types, conditions,
        os.path.join(args.output_dir, "regulon_heatmap_condition_sorted.png"),
        args.n_top,
    )

    print("\nDone!")


if __name__ == "__main__":
    main()
