library(Seurat)
library(SoupX)
library(celda)
library(SingleCellExperiment)
library(dplyr)
library(ggplot2)
library(stringr)
library(harmony)

# ------------------------------------------------------------
# 01 LOAD RAW MATRICES
# ------------------------------------------------------------

data_dir <- "/home/lab-user/data/Epi2me_raw_analysis/raw"

subdirs <- list.dirs(
  data_dir,
  recursive = FALSE
)

seurat_objects <- list()

for (subdir in subdirs) {

  matrix_file  <- file.path(subdir, "matrix.mtx.gz")
  feature_file <- file.path(subdir, "features.tsv.gz")
  barcode_file <- file.path(subdir, "barcodes.tsv.gz")

  if (
    file.exists(matrix_file) &
    file.exists(feature_file) &
    file.exists(barcode_file)
  ) {

    dataset_name <- basename(subdir)

    matrix_data <- ReadMtx(
      mtx = matrix_file,
      cells = barcode_file,
      features = feature_file
    )

    # --------------------------------------------------------
    # CREATE SEURAT OBJECT
    # RNA assay = RAW COUNTS
    # --------------------------------------------------------

    seu <- CreateSeuratObject(
      counts = matrix_data,
      project = dataset_name
    )

    # --------------------------------------------------------
    # BASIC METADATA
    # --------------------------------------------------------

    seu$condition <- dataset_name

    seu$sample <- ifelse(
      str_detect(dataset_name, "adeno"),
      "adeno",
      "sham"
    )

    seurat_objects[[dataset_name]] <- seu

  } else {

    message(
      paste("Missing files in:", subdir)
    )
  }
}

# ------------------------------------------------------------
# 02 QC + DECONTX (PER SAMPLE)
# KEEP:
#   RNA assay     = RAW counts
#   decontX assay = CLEAN counts
# ------------------------------------------------------------

seurat_objects_decont <- lapply(
  names(seurat_objects),
  function(nm) {

    seu <- seurat_objects[[nm]]

    # --------------------------------------------------------
    # QC METRICS
    # --------------------------------------------------------

    seu$log10GenesPerUMI <-
      log10(seu$nFeature_RNA) /
      log10(seu$nCount_RNA)

    seu$mitoRatio <-
      PercentageFeatureSet(
        seu,
        pattern = "^mt-"
      ) / 100

    # --------------------------------------------------------
    # OPTIONAL DOWNSAMPLING
    # --------------------------------------------------------
#
#    if (ncol(seu) > 20000) {
#
#      set.seed(1)
#
#      cells_keep <- sample(
#        colnames(seu),
#        20000
#      )
#
#      seu <- subset(
#        seu,
#        cells = cells_keep
#      )
#    }

    # --------------------------------------------------------
    # QC FILTERING
    # --------------------------------------------------------

    seu <- subset(
      seu,
      subset =
        nCount_RNA > 400 &
        nFeature_RNA > 200 &
        nFeature_RNA < 10000 &
        log10GenesPerUMI > 0.80 &
        mitoRatio < 0.15
    )

    # --------------------------------------------------------
    # PRELIMINARY CLUSTERING
    # IMPORTANT FOR GOOD DECONTX
    # --------------------------------------------------------

    DefaultAssay(seu) <- "RNA"

    seu <- NormalizeData(seu)

    seu <- FindVariableFeatures(
      seu,
      selection.method = "vst",
      nfeatures = 3000
    )

    seu <- ScaleData(
      seu,
      vars.to.regress = "mitoRatio"
    )

    seu <- RunPCA(
      seu,
      npcs = 30
    )

    seu <- FindNeighbors(
      seu,
      dims = 1:30
    )

    seu <- FindClusters(
      seu,
      resolution = 0.5
    )

    # --------------------------------------------------------
    # DECONTX
    # --------------------------------------------------------

    sce <- as.SingleCellExperiment(seu)

    sce <- decontX(
      sce,
      z = seu$seurat_clusters
    )

    # --------------------------------------------------------
    # SAVE CONTAMINATION SCORES
    # --------------------------------------------------------

    seu$decontX_contamination <-
      colData(sce)$decontX_contamination

    # --------------------------------------------------------
    # CREATE NEW ASSAY
    # RNA assay remains RAW
    # decontX assay contains corrected counts
    # --------------------------------------------------------

    seu[["decontX"]] <- CreateAssayObject(
      counts = assays(sce)$decontXcounts
    )

    # --------------------------------------------------------
    # RETURN OBJECT
    # --------------------------------------------------------

    return(seu)
  }
)

names(seurat_objects_decont) <- names(seurat_objects)

# ------------------------------------------------------------
# 03 MERGE OBJECTS
# BOTH ASSAYS ARE PRESERVED:
#   RNA
#   decontX
# ------------------------------------------------------------

merged_seurat <- merge(
  x = seurat_objects_decont[[1]],
  y = seurat_objects_decont[-1],
  add.cell.ids = names(seurat_objects_decont),
  project = "IntegratedProject"
)

# ------------------------------------------------------------
# 04 METADATA CLEANUP
# ------------------------------------------------------------

metadata <- merged_seurat@meta.data

metadata$cells <- rownames(metadata)

metadata <- metadata %>%
  dplyr::rename(
    seq_folder = orig.ident,
    nUMI = nCount_RNA,
    nGene = nFeature_RNA
  )

metadata$sample <- ifelse(
  str_detect(metadata$seq_folder, "adeno"),
  "adeno",
  "sham"
)

merged_seurat@meta.data <- metadata

# ------------------------------------------------------------
# 05 CHOOSE ASSAY FOR ANALYSIS
# ------------------------------------------------------------

# RAW COUNTS
# DefaultAssay(merged_seurat) <- "RNA"

# DECONTAMINATED COUNTS
DefaultAssay(merged_seurat) <- "decontX"

# ------------------------------------------------------------
# 06 NORMALIZATION + DIMENSION REDUCTION
# USING CURRENT DEFAULT ASSAY
# ------------------------------------------------------------

merged_seurat <- NormalizeData(
  merged_seurat
)

merged_seurat <- FindVariableFeatures(
  merged_seurat,
  selection.method = "vst",
  nfeatures = 3000
)

merged_seurat <- ScaleData(
  merged_seurat,
  vars.to.regress = "mitoRatio"
)

merged_seurat <- RunPCA(
  merged_seurat,
  npcs = 30
)

# ------------------------------------------------------------
# 07 HARMONY INTEGRATION
# ------------------------------------------------------------

merged_seurat <- RunHarmony(
  object = merged_seurat,
  group.by.vars = "condition",
  dims.use = 1:30
)

# ------------------------------------------------------------
# 08 UMAP + GRAPH + CLUSTERING
# ------------------------------------------------------------

merged_seurat <- RunUMAP(
  merged_seurat,
  reduction = "harmony",
  dims = 1:30
)

merged_seurat <- FindNeighbors(
  merged_seurat,
  reduction = "harmony",
  dims = 1:30
)

merged_seurat <- FindClusters(
  merged_seurat,
  resolution = 0.5
)

# ------------------------------------------------------------
# FINAL OBJECT STRUCTURE
# ------------------------------------------------------------

# Assays(merged_seurat)
#
# [1] "RNA" "decontX"
#
# RNA:
#   RAW counts
#
# decontX:
#   decontaminated counts
#
# SWITCH BETWEEN ASSAYS:
#
# DefaultAssay(merged_seurat) <- "RNA"
#
# or
#
# DefaultAssay(merged_seurat) <- "decontX"
#
# ------------------------------------------------------------