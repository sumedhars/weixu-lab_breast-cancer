library(Seurat)

wt <- Read10X_h5("CTR9_snRNASeq/WT/filtered_feature_bc_matrix.h5")
ko <- Read10X_h5("CTR9_snRNASeq/KO/filtered_feature_bc_matrix.h5")

wt <- CreateSeuratObject(wt, project="WT_DM")
ko <- CreateSeuratObject(ko, project="KO_DM")

wt$sample <- "WT_DM"
ko$sample <- "KO_DM"

obj <- merge(wt, y=ko) 

obj <- NormalizeData(obj)
obj <- FindVariableFeatures(obj)
obj <- ScaleData(obj)
obj <- RunPCA(obj)

obj <- FindNeighbors(obj, dims = 1:20)
obj <- FindClusters(obj, resolution = 0.5)
obj <- RunUMAP(obj, dims = 1:20)

DimPlot(obj)

DotPlot(obj, features = c("Ncr1", "Cd8a", "Cd3e", "Cd4", "Foxp3", "Cd79a", "Kit", "Csf1r", "Csf3r", "S100a8", "Epcam","Krt14"))

m <- FindAllMarkers(obj)
