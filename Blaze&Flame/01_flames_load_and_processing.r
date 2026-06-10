# =====================================================
# LIBRARIES
# =====================================================

library(Seurat)
library(Matrix)
library(rtracklayer)
library(SingleCellExperiment)
library(celda)
library(harmony)
library(DoubletFinder)
library(dplyr)
library(tibble)

# =====================================================
# PATHS
# =====================================================

data_dir <- "/media/user/8Tb/blaze_analysis/data_analysis/"
gtf_file <- file.path(data_dir, "genes.gtf")

subdirs <- list.dirs(data_dir, recursive = FALSE)

seurat_objects <- list()

# =====================================================
# LOAD GTF AND CREATE ENSG -> SYMBOL MAP
# =====================================================

gtf <- rtracklayer::import(gtf_file)

gtf_df <- as.data.frame(gtf)

gtf_genes <- gtf_df[gtf_df$type == "gene", ]

gene_map <- data.frame(
  gene_id = gtf_genes$gene_id,
  gene_symbol = gtf_genes$gene_name
)

gene_map <- gene_map[
  !is.na(gene_map$gene_id) &
    !is.na(gene_map$gene_symbol) &
    gene_map$gene_symbol != "",
]

gene_map <- gene_map[
  !duplicated(gene_map$gene_id),
]

map <- gene_map$gene_symbol
names(map) <- gene_map$gene_id

run_doubletfinder_custom <- function(
    seu_sample_subset,
    multiplet_rate = NULL
){

  message("Running DoubletFinder on: ",
          unique(seu_sample_subset$condition))

  # ================================================= #
  # 1. Estimate multiplet rate
  # ================================================= #

  if (is.null(multiplet_rate)) {

    multiplet_rates_10x <- data.frame(
      Multiplet_rate = c(
        0.004, 0.008, 0.0160, 0.023,
        0.031, 0.039, 0.046, 0.054,
        0.061, 0.069, 0.076
      ),
      Recovered_cells = c(
        500, 1000, 2000, 3000,
        4000, 5000, 6000, 7000,
        8000, 9000, 10000
      )
    )

    multiplet_rate <- multiplet_rates_10x %>%
      dplyr::filter(Recovered_cells < ncol(seu_sample_subset)) %>%
      dplyr::slice(which.max(Recovered_cells)) %>%
      dplyr::pull(Multiplet_rate)

    if (length(multiplet_rate) == 0 || is.na(multiplet_rate)) {
      multiplet_rate <- 0.01
    }

    message("Estimated multiplet rate: ", multiplet_rate)
  }

  # ================================================= #
  # 2. Determine PCs safely
  # ================================================= #

  stdv <- seu_sample_subset[["pca"]]@stdev
  percent_stdv <- (stdv / sum(stdv)) * 100
  cumulative <- cumsum(percent_stdv)

  co1 <- which(cumulative > 90 & percent_stdv < 5)[1]

  diff_drop <- c(0, diff(percent_stdv))
  co2 <- which(diff_drop < -0.1)[1]

  min_pc <- min(co1, co2, na.rm = TRUE)

  if (is.na(min_pc) || min_pc < 10) {
    min_pc <- 10
  }

  message("Using PCs: 1:", min_pc)

  # ================================================= #
  # 3. Parameter sweep
  # ================================================= #

  sweep_list <- paramSweep(
    seu_sample_subset,
    PCs = 1:min_pc,
    sct = FALSE
  )

  sweep_stats <- summarizeSweep(sweep_list)
  bcmvn <- find.pK(sweep_stats)

  # ================================================= #
  # 4. FIX pK selection (robust)
  # ================================================= #

  bcmvn <- as.data.frame(bcmvn)

  bcmvn$BCmetric <- suppressWarnings(
    as.numeric(as.character(bcmvn$BCmetric))
  )

  bcmvn <- bcmvn[!is.na(bcmvn$BCmetric), ]

  optimal.pk <- bcmvn$pK[which.max(bcmvn$BCmetric)]
  optimal.pk <- as.numeric(as.character(optimal.pk))

  message("Optimal pK: ", optimal.pk)

  # ================================================= #
  # 5. FIX: safe cluster handling (IMPORTANT)
  # ================================================= #

  annotations <- seu_sample_subset$seurat_clusters
  annotations <- as.character(annotations)

  homotypic.prop <- modelHomotypic(annotations)

  nExp.poi <- round(multiplet_rate * ncol(seu_sample_subset))
  nExp.poi.adj <- round(nExp.poi * (1 - homotypic.prop))

  message("Expected doublets: ", nExp.poi.adj)

  # ================================================= #
  # 6. Clean metadata (prevents xtfrm.data.frame)
  # ================================================= #

  bad_cols <- sapply(seu_sample_subset@meta.data, is.data.frame)
  if (any(bad_cols)) {
    seu_sample_subset@meta.data <- seu_sample_subset@meta.data[, !bad_cols, drop = FALSE]
  }

  # ================================================= #
  # 7. Ensure PCA is valid (critical for DF stability)
  # ================================================= #

  if (any(is.na(seu_sample_subset[["pca"]]@stdev))) {
    seu_sample_subset <- NormalizeData(seu_sample_subset)
    seu_sample_subset <- FindVariableFeatures(seu_sample_subset)
    seu_sample_subset <- ScaleData(seu_sample_subset)
    seu_sample_subset <- RunPCA(seu_sample_subset, verbose = FALSE)
  }

  # ================================================= #
  # 8. Run DoubletFinder
  # ================================================= #

  seu_sample_subset <- doubletFinder(
    seu = seu_sample_subset,
    PCs = 1:min_pc,
    pN = 0.25,
    pK = optimal.pk,
    nExp = nExp.poi.adj,
    reuse.pANN = NULL,
    sct = FALSE
  )

  # ================================================= #
  # 9. Rename output column safely
  # ================================================= #

  df_col <- grep(
    "DF.classifications",
    colnames(seu_sample_subset@meta.data),
    value = TRUE
  )

  if (length(df_col) > 0) {
    colnames(seu_sample_subset@meta.data)[
      colnames(seu_sample_subset@meta.data) == df_col
    ] <- "doublet_finder"
  }

  return(seu_sample_subset)
}

# =====================================================
# PROCESS EACH SAMPLE
# =====================================================

for (subdir in subdirs) {

  matrix_file <- file.path(
    subdir,
    "gene_count.csv"
  )

  if (!file.exists(matrix_file)) {

    message("Missing file in: ", subdir)

    next
  }

  dataset_name <- basename(subdir)

  message("======================================")
  message("Processing: ", dataset_name)
  message("======================================")

  # ------------------------------------------------- #
  # READ MATRIX
  # ------------------------------------------------- #

  matrix_data <- read.table(
    file = matrix_file,
    header = TRUE,
    row.names = 1,
    sep = ",",
    check.names = FALSE
  )

  matrix_data <- as.matrix(matrix_data)

  matrix_data[is.na(matrix_data)] <- 0

  # ------------------------------------------------- #
  # CONVERT ENSG IDs -> GENE SYMBOLS
  # ------------------------------------------------- #

  gene_ids <- sub(
    "\\..*",
    "",
    rownames(matrix_data)
  )

  gene_symbols <- map[gene_ids]

  valid <- !is.na(gene_symbols) &
    gene_symbols != ""

  matrix_data <- matrix_data[valid, ]

  gene_symbols <- gene_symbols[valid]

  rownames(matrix_data) <- gene_symbols

  # collapse duplicated symbols

  matrix_data <- rowsum(
    matrix_data,
    group = rownames(matrix_data)
  )

  matrix_data <- Matrix(
    matrix_data,
    sparse = TRUE
  )

  # ------------------------------------------------- #
  # CREATE SEURAT OBJECT
  # ------------------------------------------------- #

  seu <- CreateSeuratObject(
    counts = matrix_data,
    project = dataset_name
  )

  seu$condition <- dataset_name

  # ------------------------------------------------- #
  # QC
  # ------------------------------------------------- #

  seu$log10GenesPerUMI <- ifelse(
    seu$nCount_RNA == 0,
    0,
    log10(seu$nFeature_RNA) /
      log10(seu$nCount_RNA)
  )

  seu$mitoRatio <- PercentageFeatureSet(
    seu,
    pattern = "^mt-|^Mt-"
  ) / 100

  seu <- subset(
    seu,
    subset =
      nCount_RNA > 1000 &
      nFeature_RNA > 400 &
      nFeature_RNA < 8000 &
      log10GenesPerUMI > 0.80 &
      mitoRatio < 0.35 #25
  )

  message(
    "Cells after QC: ",
    ncol(seu)
  )

  # ------------------------------------------------- #
  # PREPROCESSING
  # ------------------------------------------------- #

  seu <- NormalizeData(seu)

  seu <- FindVariableFeatures(
    seu,
    selection.method = "vst",
    nfeatures = min(2000, nrow(seu))
  )

  seu <- ScaleData(seu)

  seu <- RunPCA(
    seu,
    verbose = FALSE
  )

  seu <- FindNeighbors(
    seu,
    dims = 1:20
  )

  seu <- FindClusters(
    seu,
    resolution = 0.5
  )

  # ------------------------------------------------- #
  # DOUBLETFINDER
  # ------------------------------------------------- #

  seu <- run_doubletfinder_custom(seu)

  message(
    "Doublets detected: ",
    sum(seu$doublet_finder == "Doublet")
  )

  seu <- subset(
    seu,
    subset = doublet_finder == "Singlet"
  )

  message(
    "Cells after DoubletFinder: ",
    ncol(seu)
  )

  # ------------------------------------------------- #
  # DECONTX
  # ------------------------------------------------- #

  sce <- as.SingleCellExperiment(seu)

  sce <- sce[, colnames(seu)]

  sce <- decontX(
    sce,
    z = seu$seurat_clusters
  )

  seu$decontX_contamination <-
    colData(sce)$decontX_contamination

  if (
    "decontX_cluster" %in%
    colnames(colData(sce))
  ) {

    seu$decontX_cluster <-
      colData(sce)$decontX_cluster
  }

  # ------------------------------------------------- #
  # FILTER HIGH CONTAMINATION
  # ------------------------------------------------- #

  seu <- subset(
    seu,
    subset = decontX_contamination < 0.05
  )

  message(
    "Cells after decontX filtering: ",
    ncol(seu)
  )

  # ------------------------------------------------- #
  # FINAL REPROCESSING
  # ------------------------------------------------- #

  seu <- NormalizeData(seu)

  seu <- FindVariableFeatures(
    seu,
    selection.method = "vst",
    nfeatures = 2000)
  

  seu <- ScaleData(seu)

  seu <- RunPCA(
    seu,
    verbose = FALSE
  )

  # ------------------------------------------------- #
  # STORE OBJECT
  # ------------------------------------------------- #

  seurat_objects[[dataset_name]] <- seu

  message(
    dataset_name,
    ": ",
    nrow(seu),
    " genes x ",
    ncol(seu),
    " cells retained"
  )
}

# =====================================================
# MERGE DATASETS
# =====================================================

message("======================================")
message("MERGING DATASETS")
message("======================================")

merged_seurat <- Reduce(
  function(x, y) merge(x, y),
  seurat_objects
)

# =====================================================
# POST-MERGE PROCESSING
# =====================================================

DefaultAssay(merged_seurat) <- "RNA"

# Only needed for Seurat v5
# merged_seurat[["RNA"]] <- JoinLayers(
#   merged_seurat[["RNA"]]
# )

merged_seurat <- NormalizeData(merged_seurat)

merged_seurat <- FindVariableFeatures(
  merged_seurat,
  selection.method = "vst",
  nfeatures = 3000
)
hvgs <- VariableFeatures(merged_seurat)
hvgs_clean <- hvgs[!grepl("^mt-|^Mt-|^Nefh$|^Calca$", hvgs)]

# 4. Sovrascrivi i geni variabili con la lista pulita
VariableFeatures(merged_seurat) <- hvgs_clean
#merged_seurat <- ScaleData(merged_seurat)
merged_seurat <- ScaleData(merged_seurat, features = VariableFeatures(merged_seurat))
merged_seurat <- RunPCA(
  merged_seurat,
  verbose = FALSE
)

# =====================================================
# HARMONY INTEGRATION
# =====================================================

merged_seurat <- RunHarmony(
  merged_seurat,
  group.by.vars = "condition"
)

# =====================================================
# UMAP + CLUSTERING
# =====================================================

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

# =====================================================
# SAVE
# =====================================================

saveRDS(
  merged_seurat,
  file = file.path(
    data_dir,
    "merged_seurat.rds"
  )
)

message("======================================")
message("DONE.")
message("======================================")

