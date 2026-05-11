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
merged_seurat <- FindVariableFeatures(merged_seurat, selection.method = "vst", nfeatures = 3000)
merged_seurat <- ScaleData(merged_seurat, vars.to.regress = c("mitoRatio"))
merged_seurat <- RunPCA(merged_seurat, npcs = 30)
# Run Harmony
merged_seurat <- RunHarmony(
  object = merged_seurat,
  group.by.vars = "ident", # Adjust this based on your batch metadata
  dims.use = 1:30
)

# Update embeddings for downstream use
#resolutions <- c(0.1, 0.2, 0.3, 0.4, 0.5, 1.0) 
resolutions <- c(0.5) 
merged_seurat <- RunUMAP(merged_seurat, reduction = "harmony", dims = 1:30)
merged_seurat <- FindNeighbors(merged_seurat, reduction = "harmony", dims = 1:30)
merged_seurat <- FindClusters(merged_seurat, resolution = resolutions)

# Clusters polishing
cluster_sizes <- table(merged_seurat$seurat_clusters)
# scegli soglia
threshold <- 50 #20 for exploration 
small_clusters <- names(cluster_sizes[cluster_sizes < threshold])
merged_seurat_filtered <- subset(
  merged_seurat,
  subset = !(seurat_clusters %in% small_clusters)
)

# Re-clustering after filtering small clusters
merged_seurat_filtered <- FindNeighbors(merged_seurat_filtered, dims = 1:20)
merged_seurat_filtered <- FindClusters(merged_seurat_filtered, resolution = 0.5)


Idents(merged_seurat) <- merged_seurat$RNA_snn_res.0.5 # da decidere
#merged_seurat$seurat_clusters <- merged_seurat$RNA_snn_res.0.5
Idents(merged_seurat_filtered) <- merged_seurat_filtered$seurat_clusters # da decidere
merged_seurat_filtered = JoinLayers(merged_seurat_filtered)
markers <- FindAllMarkers(object = merged_seurat_filtered, logfc.threshold = 0.1, only.pos = TRUE, test.use="wilcox", min.pct = 0.01,assay="RNA")