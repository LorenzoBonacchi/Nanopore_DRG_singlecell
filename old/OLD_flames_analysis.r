library(Seurat)
library(Matrix)
library(SingleCellExperiment)
library(celda)
library(biomaRt)

# =====================================================
# Paths
# =====================================================
data_dir <- "/media/user/8Tb/blaze_analysis/data_analysis/"
subdirs <- list.dirs(data_dir, recursive = FALSE)

seurat_objects <- list()

# =====================================================
# biomaRt (load once)
# =====================================================
mart <- useEnsembl(
  biomart = "genes",
  dataset = "mmusculus_gene_ensembl",
  mirror = "usewest"
)

# =====================================================
# Loop over datasets
# =====================================================
for (subdir in subdirs) {

  matrix_file  <- file.path(subdir, "oarfish.count.mtx")
  feature_file <- file.path(subdir, "oarfish.features.txt")
  barcode_file <- file.path(subdir, "oarfish.barcodes.txt")

  if (
    file.exists(matrix_file) &
    file.exists(feature_file) &
    file.exists(barcode_file)
  ) {

    dataset_name <- basename(subdir)

    message("Processing: ", dataset_name)

    # =====================================================
    # Load matrix
    # =====================================================
    m <- readMM(matrix_file)

    barcodes <- readLines(barcode_file)
    features <- readLines(feature_file)

    # =====================================================
    # Orientation fix
    # =====================================================
    if (
      nrow(m) == length(barcodes) &&
      ncol(m) == length(features)
    ) {

      counts <- t(m)

    } else if (
      nrow(m) == length(features) &&
      ncol(m) == length(barcodes)
    ) {

      counts <- m

    } else {

      stop(paste("Dimension mismatch in", dataset_name))
    }

    rownames(counts) <- features
    colnames(counts) <- barcodes

    # =====================================================
    # TRANSCRIPT -> GENE ID
    # =====================================================

    # remove transcript suffix
    gene_ids <- sub("_.*", "", rownames(counts))

    # remove Ensembl version number
    gene_ids <- sub("\\..*", "", gene_ids)

    # collapse transcript counts into genes
    counts_gene <- rowsum(counts, group = gene_ids)

    # =====================================================
    # GENE ID -> GENE SYMBOL
    # =====================================================
    annot <- getBM(
      attributes = c(
        "ensembl_gene_id",
        "external_gene_name"
      ),
      filters = "ensembl_gene_id",
      values = rownames(counts_gene),
      mart = mart
    )

    map <- annot$external_gene_name
    names(map) <- annot$ensembl_gene_id

    gene_symbols <- map[rownames(counts_gene)]

    # =====================================================
    # KEEP ONLY VALID SYMBOLS
    # =====================================================
    valid <- !is.na(gene_symbols) &
             gene_symbols != ""

    counts_gene <- counts_gene[valid, ]
    gene_symbols <- gene_symbols[valid]

    rownames(counts_gene) <- gene_symbols

    # =====================================================
    # COLLAPSE DUPLICATED SYMBOLS
    # =====================================================
    counts_gene <- rowsum(
      counts_gene,
      group = rownames(counts_gene)
    )

    # =====================================================
    # Remove invalid rownames
    # =====================================================
    counts_gene <- counts_gene[
      rownames(counts_gene) != "" &
      !is.na(rownames(counts_gene)),
    ]

    # =====================================================
    # Create Seurat object
    # =====================================================
    seu <- CreateSeuratObject(
      counts = counts_gene,
      project = dataset_name
    )

    seu$condition <- dataset_name

    # =====================================================
    # QC metrics
    # =====================================================
    seu$log10GenesPerUMI <-
      log10(seu$nFeature_RNA) /
      log10(seu$nCount_RNA)

    seu$mitoRatio <-
      PercentageFeatureSet(
        seu,
        pattern = "^mt-|^Mt-"
      ) / 100

    # =====================================================
    # QC filtering
    # =====================================================
    seu <- subset(
      seu,
      subset =
        nCount_RNA > 500 &
        nFeature_RNA > 400 &
        nFeature_RNA < 8000 &
        log10GenesPerUMI > 0.80 &
        mitoRatio < 0.10
    )

    # =====================================================
    # Pre-clustering for decontX
    # =====================================================
    seu <- NormalizeData(seu)

    seu <- ScaleData(seu)

    seu <- RunPCA(seu)

    seu <- FindNeighbors(
      seu,
      dims = 1:20
    )

    seu <- FindClusters(
      seu,
      resolution = 0.5
    )

    # =====================================================
    # decontX
    # =====================================================
    sce <- as.SingleCellExperiment(seu)

    # IMPORTANT alignment fix
    sce <- sce[, colnames(seu)]

    sce <- decontX(
      sce,
      z = seu$seurat_clusters
    )

    # contamination scores
    seu$decontX_contamination <-
      colData(sce)$decontX_contamination

    # decontaminated matrix
    decont_mat <- decontXcounts(sce)

    # =====================================================
    # Replace RNA assay with decontaminated counts
    # =====================================================
    seu[["RNA"]] <- CreateAssayObject(
      counts = decont_mat
    )

    DefaultAssay(seu) <- "RNA"

    # =====================================================
    # Re-run standard workflow AFTER decontamination
    # =====================================================
    seu <- NormalizeData(seu)

    seu <- FindVariableFeatures(seu)

    seu <- ScaleData(seu)

    seu <- RunPCA(seu)

    # =====================================================
    # Store object
    # =====================================================
    seurat_objects[[dataset_name]] <- seu

    message(
      dataset_name,
      ": ",
      nrow(seu),
      " genes retained"
    )

  } else {

    message("Missing files in: ", subdir)
  }
}

# =====================================================
# Merge datasets
# =====================================================
merged_seurat <- Reduce(
  function(x, y) merge(x, y),
  seurat_objects
)

# =====================================================
# Final checks
# =====================================================

# should return character(0)
grep(
  "^ENSMUSG",
  rownames(merged_seurat),
  value = TRUE
)

# check mito genes
head(
  grep(
    "^mt-|^Mt-",
    rownames(merged_seurat),
    value = TRUE
  )
)

# dimensions
merged_seurat