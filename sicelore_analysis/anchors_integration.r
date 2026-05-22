library(Seurat)
library(SingleCellExperiment)
library(celda)
library(dplyr)
library(stringr)
library(harmony)
library(tibble)
library(presto)
library(celldex)
library(SingleR)


data_dir <- "/home/lab-user/data/seurat_sicelore_analysis"
files <- list.files(
  data_dir,
  pattern = "matrix_sicelore_post\\.txt$",
  full.names = TRUE
)
seurat_objects <- list()

for (file in files) {
  dataset_name <- gsub(
    "_matrix_sicelore\\.txt",
    "",
    basename(file)
  )
  matrix_data <- read.table(
    file = file,
    header = TRUE,
    row.names = 1,
    sep = "\t",
    check.names = FALSE
  )
  matrix_data <- as.matrix(matrix_data)
  seu <- CreateSeuratObject(
    counts = matrix_data,
    project = dataset_name
  )
  seu$condition <- dataset_name
  seu$log10GenesPerUMI <-
    log10(seu$nFeature_RNA) /
    log10(seu$nCount_RNA)
  seu$mitoRatio <-
    PercentageFeatureSet(
      seu,
      pattern = "^mt-"
    ) / 100
  # QC filtering
  seu <- subset(
    seu,
    subset =
      nCount_RNA > 500 &
      nFeature_RNA > 200 &
      log10GenesPerUMI > 0.80 &
      mitoRatio < 0.1
  )
  # -------------------------------------------------
  # Preliminary clustering for decontX
  # -------------------------------------------------
  seu <- NormalizeData(seu)
  seu <- FindVariableFeatures(seu)
  seu <- ScaleData(seu)
  seu <- RunPCA(seu)
  seu <- FindNeighbors(seu, dims = 1:20)
  seu <- FindClusters(seu, resolution = 0.5)
  # -------------------------------------------------
  # decontX
  # -------------------------------------------------
  sce <- as.SingleCellExperiment(seu)
  sce <- decontX(
    sce,
    z = seu$seurat_clusters
  )
  # contamination metadata
  seu$decontX_contamination <-
    colData(sce)$decontX_contamination
  # corrected counts
  seu[["decontX"]] <- CreateAssayObject(
    counts = decontXcounts(sce)
  )
  DefaultAssay(seu) <- "decontX"
  seurat_objects[[dataset_name]] <- seu
}

merged_seurat <- merge(
  x = seurat_objects[[1]],
  y = seurat_objects[-1],
  add.cell.ids = names(seurat_objects),
  project = "IntegratedProject"
)

DefaultAssay(merged_seurat) <- "decontX"

# This override the previous condition metadata to group adeno vs sham, batchs are in original ident
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

sham <- NormalizeData(sham)
sham <- FindVariableFeatures(sham, selection.method = "vst", nfeatures = 3000)
sham <- ScaleData(sham, vars.to.regress = c("mitoRatio"))
sham <- RunPCA(sham, npcs = 30)

# Anchors integration
adeno = merged_seurat
sham = seurat_objects[[1]]
adeno <- subset(merged_seurat, condition == "adeno")
sham <- subset(merged_seurat, condition == "sham")


anchors <- FindIntegrationAnchors(
  object.list = list(adeno, sham),
  dims = 1:30
)

combined <- IntegrateData(
  anchorset = anchors,
  dims = 1:30
)
DefaultAssay(combined) <- "integrated"
# Update embeddings for downstream use
resolutions <- c(0.1, 0.2, 0.3, 0.4, 0.5, 1.0) 
#resolutions <- c(0.5) 
combined <- ScaleData(combined, vars.to.regress = c("mitoRatio"))
combined <- RunPCA(combined, npcs = 30)
combined <- RunUMAP(combined, reduction = "pca", dims = 1:30)
combined <- FindNeighbors(combined, reduction = "pca", dims = 1:30)
combined <- FindClusters(combined, resolution = resolutions)

Idents(combined) <- combined$integrated_snn_res.0.5 



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
immune_all = c(macrophages, capillary, immune)

combined = JoinLayers(combined, assay="RNA")
markers <- FindAllMarkers(object = combined, logfc.threshold = 0.1, only.pos = TRUE, test.use="wilcox", min.pct = 0.01,assay="RNA")


library(dplyr)

top25 <- markers %>%
  group_by(cluster) %>%
  slice_max(order_by = p_val_adj, n = 25, with_ties = FALSE) %>%
  arrange(cluster, desc(avg_log2FC))
  
markers24 <- FindMarkers(combined, ident.1 = 24)

markers24 %>%
  arrange(desc(avg_log2FC)) %>%
  head(20)

save(merged_seurat,top25, file = "adeno_preANCHOR_sicelore.RData")