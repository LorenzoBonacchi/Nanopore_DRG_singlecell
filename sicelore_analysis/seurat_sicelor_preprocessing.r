library(Seurat)
library(SingleCellExperiment)
library(celda)
library(dplyr)
library(stringr)



data_dir <- "/home/lab-user/data/seurat_sicelore_analysis"
files <- list.files(
  data_dir,
  pattern = "matrix_sicelore\\.txt$",
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
      nCount_RNA > 200 &
      nFeature_RNA > 200 &
      log10GenesPerUMI > 0.80 &
      mitoRatio < 0.15
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