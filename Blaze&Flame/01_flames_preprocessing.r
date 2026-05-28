library(Seurat)
library(Matrix)
library(rtracklayer)
library(SingleCellExperiment)
library(celda)
library(harmony)

# ------------------------------------------------- #
# Path set to Blaze & Flame data
# ------------------------------------------------- #
data_dir <- "/media/user/8Tb/blaze_analysis/data_analysis/"
gtf_file <- file.path(data_dir, "genes.gtf")

subdirs <- list.dirs(data_dir, recursive = FALSE)
seurat_objects <- list()

# ------------------------------------------------- #
# load gtf to convrt ENSG to gene symbols
# ------------------------------------------------- #
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

gene_map <- gene_map[!duplicated(gene_map$gene_id), ]
map <- gene_map$gene_symbol
names(map) <- gene_map$gene_id

for (subdir in subdirs) {

  matrix_file <- file.path(subdir, "gene_count.csv")

  if (!file.exists(matrix_file)) {
    message("Missing file in: ", subdir)
    next
  }

  dataset_name <- basename(subdir)
  message("Processing: ", dataset_name)


  matrix_data <- read.table(
    file = matrix_file,
    header = TRUE,
    row.names = 1,
    sep = ",",
    check.names = FALSE
  )

  matrix_data <- as.matrix(matrix_data)
  matrix_data[is.na(matrix_data)] <- 0 # handle NAs if present

  # ------------------------------------------------- #
  # Convert ENSG IDs to gene symbols
  gene_ids <- sub("\\..*", "", rownames(matrix_data))
  gene_symbols <- map[gene_ids]

  valid <- !is.na(gene_symbols) & gene_symbols != ""

  matrix_data <- matrix_data[valid, ]
  gene_symbols <- gene_symbols[valid]
  rownames(matrix_data) <- gene_symbols

  # collapse duplicated symbols
  matrix_data <- rowsum(matrix_data, group = rownames(matrix_data))
  matrix_data <- Matrix(matrix_data, sparse = TRUE)


  seu <- CreateSeuratObject(
    counts = matrix_data,
    project = dataset_name
  )
  seu$condition <- dataset_name

  # ------------------------------------------------- #
  # QC
  seu$log10GenesPerUMI <- ifelse(
    seu$nCount_RNA == 0,
    0,
    log10(seu$nFeature_RNA) / log10(seu$nCount_RNA)
  )

  seu$mitoRatio <- PercentageFeatureSet(
    seu,
    pattern = "^mt-|^Mt-"
  ) / 100

  seu <- subset(
    seu,
    subset =
      nCount_RNA > 500 &
      nFeature_RNA > 500 &
      nFeature_RNA < 8000 &
      log10GenesPerUMI > 0.80 &
      mitoRatio < 0.25
  )

  # ------------------------------------------------- #
  # Setting up for decontX
  seu <- NormalizeData(seu)

  seu <- FindVariableFeatures(
    seu,
    selection.method = "vst",
    nfeatures = min(2000, nrow(seu))
  )

  seu <- tryCatch({
    ScaleData(seu)
  }, error = function(e) {
    message("ScaleData failed: ", dataset_name)
    seu
  })

  seu <- RunPCA(seu, verbose = FALSE)

  seu <- FindNeighbors(seu, dims = 1:20)
  seu <- FindClusters(seu, resolution = 0.5)

  # ------------------------------------------------- #
  # decontX (celda)

  sce <- as.SingleCellExperiment(seu)
  sce <- sce[, colnames(seu)]

  sce <- decontX(
    sce,
    z = seu$seurat_clusters
  )

  seu$decontX_contamination <- colData(sce)$decontX_contamination

  if ("decontX_cluster" %in% colnames(colData(sce))) {
    seu$decontX_cluster <- colData(sce)$decontX_cluster
  }

  # ------------------------------------------------- #
  # Post-decontX processing
  seu <- NormalizeData(seu)

  seu <- FindVariableFeatures(
    seu,
    selection.method = "vst",
    nfeatures = min(2000, nrow(seu))
  )
  seu <- ScaleData(seu)
  seu <- RunPCA(seu, verbose = FALSE)


  seurat_objects[[dataset_name]] <- seu

  message(
    dataset_name, ": ",
    nrow(seu), " genes x ",
    ncol(seu), " cells retained"
  )
}

# =====================================================
# MERGE DATASETS
# =====================================================
merged_seurat <- Reduce(
  function(x, y) merge(x, y),
  seurat_objects
)

# =====================================================
# SAVE
# =====================================================
saveRDS(
  merged_seurat,
  file = file.path(data_dir, "merged_seurat.rds")
)

message("DONE.")

merged_seurat = subset(merged_seurat, subset = decontX_contamination < 0.2)

DefaultAssay(merged_seurat) <- "RNA"

merged_seurat[["RNA"]] <- JoinLayers(merged_seurat[["RNA"]])







