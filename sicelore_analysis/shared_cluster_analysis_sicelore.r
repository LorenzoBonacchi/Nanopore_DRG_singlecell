# tabella cluster vs condition
tab <- table(merged_seurat$condition, merged_seurat$seurat_clusters)

# condizioni adeno e sham
adeno_rows <- grep("^adeno", rownames(tab))
sham_rows  <- grep("^sham", rownames(tab))

# cluster presenti in almeno un adeno E almeno uno sham
shared_clusters <- colnames(tab)[
  colSums(tab[adeno_rows, ] > 0) > 0 &
  colSums(tab[sham_rows, ] > 0) > 0
]

shared_clusters

merged_subset <- subset(
  merged_seurat,
  idents = shared_clusters
)


library(Seurat)
merged_subset <- NormalizeData(merged_subset)
merged_subset <- FindVariableFeatures(merged_subset, selection.method = "vst", nfeatures = 3000)
merged_subset <- ScaleData(merged_subset, vars.to.regress = c("mitoRatio"))
merged_subset <- RunPCA(merged_subset, npcs = 30)

merged_subset <- RunHarmony(
  object = merged_subset,
  group.by.vars = "orig.ident", # Adjust this based on your batch metadata
  dims.use = 1:30
)


# Update embeddings for downstream use
resolutions <- c(0.1, 0.2, 0.3, 0.4, 0.5, 1.0) 
#resolutions <- c(0.5) 
merged_subset <- RunUMAP(merged_subset, reduction = "harmony", dims = 1:30)
merged_subset <- FindNeighbors(merged_subset, reduction = "harmony", dims = 1:30)
merged_subset <- FindClusters(merged_subset, resolution = resolutions)



Idents(merged_subset) <- merged_subset$decontX_snn_res.0.3 # da decidere
merged_subset[["RNA"]] <- JoinLayers(merged_subset[["RNA"]])
markers <- FindAllMarkers(
  object = merged_subset,
  logfc.threshold = 0.1,
  only.pos = TRUE,
  test.use = "wilcox",
  min.pct = 0.01,
  assay = "RNA"
)


top20_markers <- markers %>%
  filter(p_val_adj < 0.05) %>%
  group_by(cluster) %>%
  arrange(desc(avg_log2FC), .by_group = TRUE) %>%
  slice_head(n = 20) %>%
  ungroup()