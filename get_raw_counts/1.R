library(Seurat)
library(Matrix)

setwd("~/Documents/Sumedha/weixu-lab_breast-cancer/CTR9_snRNASeq/")
load("~/Documents/Sumedha/weixu-lab_breast-cancer/CTR9_snRNASeq/CTR9_snRNASeq_seurat_obj_annot.rda")

# raw counts as a separate file
counts_matrix <- GetAssayData(seurat, layer = "counts")

# market matrix (efficient for sparse data)
writeMM(counts_matrix, "CTR9_counts.mtx")

# Save gene names and cell barcodes
writeLines(rownames(counts_matrix), "CTR9_genes.txt")
writeLines(colnames(counts_matrix), "CTR9_cells.txt")

print(paste("Matrix dimensions:", nrow(counts_matrix), "genes x", ncol(counts_matrix), "cells"))