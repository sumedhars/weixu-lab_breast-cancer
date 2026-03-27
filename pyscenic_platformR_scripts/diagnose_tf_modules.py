#!/usr/bin/env python
"""
Diagnostic: Why are Tfap2b, Foxa1, Esr1 missing from SCENIC regulons?
======================================================================
Checks:
  1. Do modules exist for these TFs after modules_from_adjacencies?
  2. How well do module genes overlap with the ranking database genes?
  3. Are these TFs' genes present in the databases at all?
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
import pickle
import argparse

from ctxcore.rnkdb import FeatherRankingDatabase as RankingDatabase


def parse_args():
    parser = argparse.ArgumentParser(
        description="Diagnose why specific TFs are absent from SCENIC regulons"
    )
    parser.add_argument("--modules_pkl", required=True,
                        help="Path to saved modules.pkl from Phase I")
    parser.add_argument("--db_dir", required=True,
                        help="Directory containing .feather ranking databases")
    parser.add_argument("--motif_annotations", required=True,
                        help="Path to motif annotations .tbl file")
    parser.add_argument("--tfs", nargs="+", default=["Tfap2b", "Foxa1", "Esr1"],
                        help="TFs of interest to diagnose (default: Tfap2b Foxa1 Esr1)")
    parser.add_argument("--adjacencies_tsv", default=None,
                        help="(Optional) Path to adjacencies.tsv to also check raw GRN edges")
    return parser.parse_args()


def main():
    args = parse_args()
    tfs = args.tfs

    print("=" * 70)
    print("SCENIC TF Regulon Diagnostic")
    print(f"TFs of interest: {tfs}")
    print("=" * 70)

    # -----------------------------------------------------------------
    # 1. Load modules and check TF presence
    # -----------------------------------------------------------------
    print(f"\n{'='*70}")
    print("CHECK 1: Do modules exist for these TFs?")
    print(f"{'='*70}")

    with open(args.modules_pkl, "rb") as f:
        modules = pickle.load(f)

    print(f"Total modules loaded: {len(modules)}")

    for tf in tfs:
        tf_mods = [m for m in modules if m.transcription_factor == tf]
        print(f"\n  {tf}: {len(tf_mods)} module(s)")
        if len(tf_mods) == 0:
            # Check if the TF appears as transcription_factor under any name
            all_tfs_in_modules = set(m.transcription_factor for m in modules)
            close = [t for t in all_tfs_in_modules if tf.lower() in t.lower()]
            if close:
                print(f"    Possible name variants found in modules: {close}")
            else:
                print(f"    NOT FOUND in any module. This TF was likely filtered "
                      f"during modules_from_adjacencies.")
        else:
            for m in tf_mods:
                genes = list(m.genes)
                print(f"    Module: {m.name}")
                print(f"      Genes: {len(genes)}")
                print(f"      First 20 genes: {genes[:20]}")

    # -----------------------------------------------------------------
    # 2. Load ranking databases
    # -----------------------------------------------------------------
    print(f"\n{'='*70}")
    print("CHECK 2: Gene overlap between modules and ranking databases")
    print(f"{'='*70}")

    db_fnames = sorted(glob.glob(os.path.join(args.db_dir, "*.feather")))
    if not db_fnames:
        print(f"ERROR: No .feather databases found in {args.db_dir}")
        sys.exit(1)

    dbs = []
    db_gene_sets = []
    for f in db_fnames:
        db = RankingDatabase(fname=f, name=os.path.splitext(os.path.basename(f))[0])
        db_genes = set(db.genes)
        dbs.append(db)
        db_gene_sets.append(db_genes)
        print(f"\n  Database: {db.name}")
        print(f"    Genes in DB: {len(db_genes)}")

    # Union of all DB genes
    all_db_genes = set()
    for gs in db_gene_sets:
        all_db_genes |= gs
    print(f"\n  Union of all DB genes: {len(all_db_genes)}")

    # -----------------------------------------------------------------
    # 3. For each TF, check module gene overlap with each database
    # -----------------------------------------------------------------
    for tf in tfs:
        tf_mods = [m for m in modules if m.transcription_factor == tf]
        if not tf_mods:
            print(f"\n  {tf}: No modules — skipping overlap check")
            continue

        for m in tf_mods:
            mod_genes = set(m.genes)
            print(f"\n  {tf} | Module '{m.name}' ({len(mod_genes)} genes):")

            for db, db_genes in zip(dbs, db_gene_sets):
                overlap = mod_genes & db_genes
                pct = 100 * len(overlap) / len(mod_genes) if mod_genes else 0
                print(f"    vs {db.name}: {len(overlap)}/{len(mod_genes)} genes "
                      f"overlap ({pct:.1f}%)")

            # Union overlap
            union_overlap = mod_genes & all_db_genes
            union_pct = 100 * len(union_overlap) / len(mod_genes) if mod_genes else 0
            print(f"    vs ALL DBs (union): {len(union_overlap)}/{len(mod_genes)} "
                  f"({union_pct:.1f}%)")

            # Show genes NOT in any database
            missing = mod_genes - all_db_genes
            if missing:
                missing_sorted = sorted(missing)
                print(f"    Genes NOT in any DB ({len(missing)}): "
                      f"{missing_sorted[:30]}{'...' if len(missing) > 30 else ''}")

    # -----------------------------------------------------------------
    # 4. Check if the TFs themselves are in the databases as genes
    # -----------------------------------------------------------------
    print(f"\n{'='*70}")
    print("CHECK 3: Are the TFs themselves present as genes in the databases?")
    print(f"{'='*70}")

    for tf in tfs:
        for db, db_genes in zip(dbs, db_gene_sets):
            present = tf in db_genes
            print(f"  {tf} in {db.name}: {'YES' if present else 'NO'}")

    # -----------------------------------------------------------------
    # 5. Check motif annotation counts for these TFs
    # -----------------------------------------------------------------
    print(f"\n{'='*70}")
    print("CHECK 4: Motif annotation counts for these TFs")
    print(f"{'='*70}")

    with open(args.motif_annotations, "r") as f:
        header = f.readline()
        lines = f.readlines()

    for tf in tfs:
        # Case-sensitive whole-word search in each line
        hits = [l for l in lines if tf in l]
        print(f"  {tf}: {len(hits)} annotation rows")

    # -----------------------------------------------------------------
    # 6. (Optional) Check adjacencies for edge counts / importance
    # -----------------------------------------------------------------
    if args.adjacencies_tsv:
        import pandas as pd

        print(f"\n{'='*70}")
        print("CHECK 5: GRNBoost2 adjacency stats for these TFs")
        print(f"{'='*70}")

        adj = pd.read_csv(args.adjacencies_tsv, sep="\t")
        for tf in tfs:
            tf_edges = adj[adj["TF"] == tf]
            print(f"\n  {tf}: {len(tf_edges)} edges in adjacencies")
            if len(tf_edges) > 0:
                print(f"    Importance — mean: {tf_edges['importance'].mean():.4f}, "
                      f"median: {tf_edges['importance'].median():.4f}, "
                      f"max: {tf_edges['importance'].max():.4f}")
                top10 = tf_edges.nlargest(10, "importance")[["target", "importance"]]
                print(f"    Top 10 targets:\n{top10.to_string(index=False)}")

    # -----------------------------------------------------------------
    print(f"\n{'='*70}")
    print("DIAGNOSTIC COMPLETE")
    print(f"{'='*70}")


if __name__ == "__main__":
    main()
