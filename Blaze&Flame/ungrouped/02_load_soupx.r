
for (subdir in subdirs) {
  matrix_file <- file.path(subdir, "gene_count.csv")
  if (!file.exists(matrix_file)) next
  dataset_name <- basename(subdir)
  message("Processing: ", dataset_name)

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
  # GENE MAP
  # ------------------------------------------------- #

  gene_ids <- sub("\\..*", "", rownames(matrix_data))
  gene_symbols <- map[gene_ids]

  valid <- !is.na(gene_symbols) & gene_symbols != ""
  matrix_data <- matrix_data[valid, ]
  gene_symbols <- gene_symbols[valid]

  rownames(matrix_data) <- gene_symbols
  matrix_data <- rowsum(matrix_data, group = rownames(matrix_data))
  matrix_data <- Matrix(matrix_data, sparse = TRUE)

  # =====================================================
  # SOUPX (ambient RNA correction)
  # =====================================================

  seu <- CreateSeuratObject(counts = matrix_data, project = dataset_name)
  seu$condition <- dataset_name

  seu <- NormalizeData(seu)
  seu <- FindVariableFeatures(seu, nfeatures = 2000)
  seu <- ScaleData(seu)
  seu <- RunPCA(seu, verbose = FALSE)

  seu <- FindNeighbors(seu, dims = 1:20)
  seu <- FindClusters(seu, resolution = 1)

  sc <- SoupChannel(
    tod = matrix_data,
    toc = matrix_data
  )

  sc <- setClusters(sc, seu$seurat_clusters)
  sc <- autoEstCont(sc)

  matrix_data_clean <- adjustCounts(sc)

  # rebuild CLEAN Seurat object
  seu <- CreateSeuratObject(counts = matrix_data_clean,
                            project = dataset_name)

  seu$condition <- dataset_name
  # =====================================================
  # STORE
  # =====================================================

  seurat_objects[[dataset_name]] <- seu
}