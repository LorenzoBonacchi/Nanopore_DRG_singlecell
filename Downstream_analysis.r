library(Seurat)
library(harmony)
library(DoubletFinder)
library(tibble)
library(presto)
load("objects_post_filtering.RData")

merged_seurat <- merge(seurat_objects_filtered[[1]], 
                       y = seurat_objects_filtered[-1], 
                       add.cell.ids = names(seurat_objects_filtered), 
                       project = "IntegratedProject")

merged_seurat$condition <- ifelse(
  grepl("adeno", merged_seurat$condition),
  "adeno",
  "sham"
)

# Preprocess the merged object
merged_seurat <- NormalizeData(merged_seurat)
merged_seurat <- FindVariableFeatures(merged_seurat, selection.method = "vst", nfeatures = 2000)
merged_seurat <- ScaleData(merged_seurat, vars.to.regress = c("percent.mt", "percent.redcell"))
merged_seurat <- RunPCA(merged_seurat, npcs = 30)
# Run Harmony
merged_seurat <- RunHarmony(
  object = merged_seurat,
  group.by.vars = "condition", # Adjust this based on your batch metadata
  dims.use = 1:30
)

# Update embeddings for downstream use
resolutions <- c(0.1, 0.2, 0.3, 0.4, 0.5, 1.0) 
merged_seurat <- RunUMAP(merged_seurat, reduction = "harmony", dims = 1:30)
merged_seurat <- FindNeighbors(merged_seurat, reduction = "harmony", dims = 1:30)
merged_seurat <- FindClusters(merged_seurat, resolution = resolutions)

Idents(merged_seurat) <- merged_seurat$RNA_snn_res.0.5 # da decidere
merged_seurat$seurat_clusters <- merged_seurat$RNA_snn_res.0.5
merged_seurat = JoinLayers(merged_seurat)
markers <- FindAllMarkers(object = merged_seurat, logfc.threshold = 0.1, only.pos = TRUE, test.use="wilcox", min.pct = 0.01,assay="RNA")

save(markers,merged_seurat,file="data_markers.RData")