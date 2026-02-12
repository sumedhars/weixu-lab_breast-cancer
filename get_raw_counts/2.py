#!/usr/bin/env python3
"""
Script to add raw counts layer to CTR9 h5ad file
"""

import scanpy as sc
import scipy.io
import pandas as pd
import numpy as np

# Load the existing h5ad file
print("Loading existing h5ad file...")
adata = sc.read_h5ad("CTR9_snRNASeq/CTR9_snRNASeq_full.h5ad")
print(f"Loaded: {adata.shape[0]} cells x {adata.shape[1]} genes")

# Load the raw counts matrix
print("\nLoading raw counts matrix...")
counts_mtx = scipy.io.mmread("CTR9_snRNASeq/CTR9_counts.mtx")

# Load gene names and cell barcodes
with open("CTR9_snRNASeq/CTR9_genes.txt", 'r') as f:
    genes = [line.strip() for line in f]

with open("CTR9_snRNASeq/CTR9_cells.txt", 'r') as f:
    cells = [line.strip() for line in f]

print(f"Raw counts matrix: {len(genes)} genes x {len(cells)} cells")

# Convert to CSR format (efficient for row operations)
counts_mtx = counts_mtx.T.tocsr()

# Check if the dimensions match
print(f"\nChecking dimensions...")
print(f"adata shape: {adata.shape}")
print(f"Raw counts shape: {counts_mtx.shape}")

print(f"\nMatching genes...")
print(f"Genes in adata: {adata.n_vars}")
print(f"Genes in raw counts: {len(genes)}")

# Create a mapping from gene names to indices
gene_to_idx = {gene: i for i, gene in enumerate(genes)}

# Find indices of adata genes in the raw counts
indices = []
matched_genes = []
for gene in adata.var_names:
    if gene in gene_to_idx:
        indices.append(gene_to_idx[gene])
        matched_genes.append(gene)

print(f"Matched {len(matched_genes)} genes out of {adata.n_vars}")

# Subset the raw counts matrix
counts_subset = counts_mtx[:, indices]


# Verify cell names match
if adata.obs_names.tolist() != cells:
    print("WARNING: Cell order doesn't match! Reordering...")
    # Create cell mapping
    cell_to_idx = {cell: i for i, cell in enumerate(cells)}
    cell_indices = [cell_to_idx[cell] for cell in adata.obs_names]
    counts_subset = counts_subset[cell_indices, :]
    print("Cells reordered to match adata")

# Add as a layer
adata.layers["raw_count"] = counts_subset
# Verify
print(f"adata.layers keys: {list(adata.layers.keys())}")
print(f"raw_count layer shape: {adata.layers['raw_count'].shape}")

# Save the updated h5ad file
output_file = "CTR9_snRNASeq/CTR9_snRNASeq_with_raw.h5ad"
adata.write_h5ad(output_file)

print("\n✓ Done! Raw counts layer added successfully.")
print(f"\nUpdated AnnData object:")
print(adata)
