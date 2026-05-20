

library(clusteree)
clustree(merged_seurat, prefix = "decontX_snn_res.")
Idents(merged_seurat) <- merged_seurat$decontX_snn_res.0.3 

# Clusters polishing
cluster_table <- table(
  merged_seurat$decontX_snn_res.0.3 ,
  merged_seurat$condition
)
adeno_counts <- cluster_table[, "adeno"]
keep_clusters <- names(adeno_counts[adeno_counts >= 50])

merged_seurat_filtered <- subset(
  merged_seurat,
  subset = decontX_snn_res.0.3 %in% keep_clusters
)

# Re-clustering after filtering small cluster
merged_seurat_filtered <- FindNeighbors(merged_seurat_filtered, dims = 1:20)
merged_seurat_filtered <- FindClusters(merged_seurat_filtered, resolution = 03)


Idents(merged_seurat_filtered) <- merged_seurat_filtered$decontX_snn_res.0.3 # da decidere
merged_seurat_filtered$seurat_clusters <- merged_seurat_filtered$decontX_snn_res.0.3

# --- Annotation ----------- #
neurons = c("Snap25","Syt1","Avil","Gap43","Nefl","Nefm")
noci = c("Trpv1", "Scn10a", "Scn9a", "Calca", "Tac1", "Asic3")
itch = c("Il31ra", "Nppb", "Sst")
proprioceptive = c("Pvalb", "Parvalbumin", "Ntrk3", "Runx3")
mecha = c("Piezo2","Mef2c")
satellite = c("Gja1", "Fabp7", "Glul", "Aqp4", "Kcnj10")
mschwann = c("Mbp", "Krox20", "Mag", "Pmp22")
nomschwann = c("Ngfr", "Sox10", "Bnc1")
fibroblasts = c("Dcn","Apod","Col1a1","Col1a2","Lum","Pdgfra")
endothelial = c("Cldn5","Ly6c1","Pecam1","Kdr","Vwf")
smooth_muscle = c("Acta2","Tagln","Myh11","Cnn1")
macrophages = c("Mrc1","Cd68","Adgre1","Cx3cr1","Lyz2")
capillary = c("S100a8","Ly6c2")
immune = c("Cd74","Cd3d","Cd3e","Trac")
rbc = c("Hba-a1","Hbb-bs","Hbb-bt")
pericytes = c("Rgs5","Pdgfrb","Des")

all = c(neurons, noci, itch, proprioceptive, mecha, satellite, mschwann, nomschwann, fibroblasts, endothelial, smooth_muscle, macrophages, capillary, immune, rbc, pericytes)
neu_all = c(neurons, noci, itch, proprioceptive, mecha)
glia_all = c(satellite, mschwann, nomschwann)
immune = c(macrophages, capillary, immune)

library(Seurat)
library(dplyr)
library(DESeq2)

# =====================================================
# SUBSET CLUSTER 9
# =====================================================

cluster9 <- subset(
  merged_seurat_filtered,
  subset = seurat_clusters == 9
)

# =====================================================
# WILCOXON DE (single-cell level)
# =====================================================

Idents(cluster9) <- "condition"

wilcox_markers <- FindMarkers(
  cluster9,
  ident.1 = "adeno",
  ident.2 = "sham",
  test.use = "wilcox",
  logfc.threshold = 0.25,
  min.pct = 0.1
)

# view top genes
head(wilcox_markers)

# save
write.csv(
  wilcox_markers,
  file = "cluster9_wilcox_markers.csv"
)

# =====================================================
# PSEUDOBULK DESEQ2
# =====================================================

# IMPORTANT:
# You need biological replicates
# Example:
# adeno1, adeno2, sham1, sham2
#
# Here we assume:
# merged_seurat_filtered$sample
# contains biological replicate names

# =====================================================
# CREATE METADATA
# =====================================================

meta <- cluster9@meta.data

# check samples
table(meta$sample, meta$condition)

# =====================================================
# EXTRACT COUNTS
# =====================================================

counts <- GetAssayData(
  cluster9,
  slot = "counts"
)

# =====================================================
# PSEUDOBULK AGGREGATION
# Sum counts per sample
# =====================================================

sample_ids <- unique(meta$sample)

pseudobulk_counts <- sapply(sample_ids, function(s) {

  cells <- rownames(meta)[meta$sample == s]

  Matrix::rowSums(counts[, cells, drop = FALSE])

})

pseudobulk_counts <- as.matrix(pseudobulk_counts)

# =====================================================
# COLDATA FOR DESEQ2
# =====================================================

coldata <- data.frame(
  sample = sample_ids,
  condition = sapply(sample_ids, function(x) {
    unique(meta$condition[meta$sample == x])
  })
)

rownames(coldata) <- coldata$sample

# =====================================================
# RUN DESEQ2
# =====================================================

dds <- DESeqDataSetFromMatrix(
  countData = pseudobulk_counts,
  colData = coldata,
  design = ~ condition
)

# optional filtering
dds <- dds[rowSums(counts(dds)) > 10, ]

dds <- DESeq(dds)

# =====================================================
# RESULTS
# =====================================================

res <- results(
  dds,
  contrast = c("condition", "adeno", "sham")
)

res <- as.data.frame(res)

# order by adjusted pvalue
res <- res[order(res$padj), ]

head(res)

# save
write.csv(
  res,
  file = "cluster9_pseudobulk_DESeq2.csv"
)