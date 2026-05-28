


# ----------------------------------------------- #
# Doublet Markers removal strategy
# ----------------------------------------------- #
library(dplyr)


merged_seurat <- NormalizeData(merged_seurat)
merged_seurat <- FindVariableFeatures(merged_seurat, nfeatures = 3000)
merged_seurat <- ScaleData(merged_seurat)

merged_seurat <- RunPCA(merged_seurat, npcs = 30)

merged_seurat <- RunHarmony(
  merged_seurat,
  group.by = "orig.ident"
)

merged_seurat <- RunUMAP(merged_seurat, reduction = "harmony",dims = 1:30)
merged_seurat <- FindNeighbors(merged_seurat, reduction = "harmony", dims = 1:30)
merged_seurat <- FindClusters(merged_seurat, resolution = 1)

marker_sets <- list(
  Neuron = c("Snap25", "Rbfox3", "Pvalb"),
  SGC = c("Fabp7", "Ednrb"),
  mySC = c("Mpz", "Mbp"),
  nmSC = c("Scn7a"),
  Immune = c("Ptprc", "Ccr2"),
  Endothelial = c("Pecam1", "Flt1"),
  Pericyte = c("Notch3", "Kcnj8"),
  Fibroblast = c("Pdgfra", "Tbx18")
)

merged_seurat <- AddModuleScore(
  merged_seurat,
  features = marker_sets,
  name = names(marker_sets)
)

score_cols <- grep("1$", colnames(merged_seurat@meta.data), value = TRUE)
score_cols

score_matrix <- merged_seurat@meta.data[, score_cols]
colnames(score_matrix) <- names(marker_sets)

merged_seurat$celltype <- apply(score_matrix, 1, function(x) {
  names(which.max(x))
})

neural <- c("Snap25", "Rbfox3", "Pvalb")

non_neural <- unlist(marker_sets[names(marker_sets) != "Neuron"])

merged_seurat <- AddModuleScore(
  merged_seurat,
  features = list(
    neural = neural,
    non_neural = non_neural
  ),
  name = c("neural_score", "nonneural_score")
)
grep("neural_score|nonneural_score", colnames(merged_seurat@meta.data), value = TRUE)

merged_seurat$doublet_score <- 
  merged_seurat$neural_score1 + merged_seurat$nonneural_score1
threshold <- quantile(merged_seurat$doublet_score, 0.95)

merged_seurat$doublet_like <- merged_seurat$doublet_score > threshold
DimPlot(merged_seurat, group.by = "celltype", label = TRUE)
DimPlot(merged_seurat, group.by = "doublet_like")
FeaturePlot(merged_seurat, features = c("neural_score1", "nonneural_score1"))
merged_seurat_clean <- merged_seurat[, !merged_seurat$doublet_like]