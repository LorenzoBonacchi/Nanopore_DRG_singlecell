library(Seurat)
library(Matrix)
library(SingleCellExperiment)
library(celda)

data_dir <- "/media/user/8Tb/blaze_analysis/data_analysis/"
subdirs <- list.dirs(data_dir, recursive = FALSE)

seurat_objects <- list()

for (subdir in subdirs) {

  matrix_file  <- file.path(subdir, "oarfish.count.mtx")
  feature_file <- file.path(subdir, "oarfish.features.txt")
  barcode_file <- file.path(subdir, "oarfish.barcodes.txt")

  if (file.exists(matrix_file) &
      file.exists(feature_file) &
      file.exists(barcode_file)) {

    dataset_name <- basename(subdir)

    # -------------------------
    # Load matrix (cells x features)
    # -------------------------
    m <- readMM(matrix_file)

    barcodes <- readLines(barcode_file)
    features <- readLines(feature_file)

    # -------------------------
    # VALIDATION
    # -------------------------
    if (nrow(m) == length(barcodes) &&
        ncol(m) == length(features)) {

      # OK: cells x features → transpose for Seurat
      counts <- t(m)

    } else if (nrow(m) == length(features) &&
               ncol(m) == length(barcodes)) {

      # already features x cells
      counts <- m

    } else {
      stop(paste("Dimension mismatch in", dataset_name))
    }

    # -------------------------
    # assign names
    # -------------------------
    rownames(counts) <- features
    colnames(counts) <- barcodes

    # -------------------------
    # Create Seurat object
    # -------------------------
    seu <- CreateSeuratObject(
      counts = counts,
      project = dataset_name
    )

    seu$condition <- dataset_name

    # -------------------------
    # QC metrics
    # -------------------------
    seu$log10GenesPerUMI <-
      log10(seu$nFeature_RNA) / log10(seu$nCount_RNA)

    seu$mitoRatio <-
      PercentageFeatureSet(seu, pattern = "^mt-") / 100

    # -------------------------
    # QC filtering
    # -------------------------
    seu <- subset(
      seu,
      subset =
        nCount_RNA > 500 &
        nFeature_RNA > 400 &
        nFeature_RNA < 8000 &
        log10GenesPerUMI > 0.80 &
        mitoRatio < 0.1
    )

    # -------------------------
    # Pre-clustering
    # -------------------------
    seu <- NormalizeData(seu)
    seu <- FindVariableFeatures(seu)
    seu <- ScaleData(seu)
    seu <- RunPCA(seu)
    seu <- FindNeighbors(seu, dims = 1:20)
    seu <- FindClusters(seu, resolution = 0.5)

    # -------------------------
    # decontX
    # -------------------------
    sce <- as.SingleCellExperiment(seu)

    sce <- decontX(
      sce,
      z = seu$seurat_clusters
    )

    seu$decontX_contamination <-
      colData(sce)$decontX_contamination

    seu[["decontX"]] <- CreateAssayObject(
      counts = decontXcounts(sce)
    )

    DefaultAssay(seu) <- "decontX"

    # -------------------------
    # store
    # -------------------------
    seurat_objects[[dataset_name]] <- seu

  } else {
    message(paste("Missing files in:", subdir))
  }
}

# -------------------------
# merge robusto
# -------------------------
merged_seurat <- Reduce(
  function(x, y) merge(x, y),
  seurat_objects
)


gene_ids <- gsub("-.*", "", features)
library(biomaRt)

mart <- useEnsembl(
  biomart = "genes",
  dataset = "mmusculus_gene_ensembl"
)

annot <- getBM(
  attributes = c("ensembl_gene_id", "external_gene_name"),
  filters = "ensembl_gene_id",
  values = gene_ids,
  mart = mart
)

gene_map <- annot$external_gene_name
names(gene_map) <- annot$ensembl_gene_id

gene_symbols <- gene_map[gene_ids]
gene_symbols[is.na(gene_symbols)] <- gene_ids[is.na(gene_symbols)]
gene_map <- annot$external_gene_name
names(gene_map) <- annot$ensembl_gene_id

gene_symbols <- gene_map[gene_ids]
gene_symbols[is.na(gene_symbols)] <- gene_ids[is.na(gene_symbols)]
library(Matrix)

counts_gene <- rowsum(counts, group = gene_symbols)