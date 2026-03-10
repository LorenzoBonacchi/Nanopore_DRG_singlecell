library(Seurat)

data_dir <- "/media/user/8Tb/scRNAseq/seurat_analysis/data_matrices/"
subdirs <- list.dirs(data_dir, recursive = FALSE)
seurat_objects <- list()

for (subdir in subdirs) {

  matrix_file <- file.path(subdir, "matrix.mtx.gz")
  feature_file <- file.path(subdir, "features.tsv.gz")
  barcode_file <- file.path(subdir, "barcodes.tsv.gz")

  if (file.exists(matrix_file) & file.exists(feature_file) & file.exists(barcode_file)) {

    dataset_name <- basename(subdir)

    matrix_data <- ReadMtx(
      mtx = matrix_file,
      cells = barcode_file,
      features = feature_file
    )

    seurat_obj <- CreateSeuratObject(
      counts = matrix_data,
      project = dataset_name
    )
    seurat_obj$condition <- dataset_name #Condition and orig.ident are the same, need to change later for batch reference
    seurat_objects[[dataset_name]] <- seurat_obj

  } else {
    message(paste("Missing files in:", subdir))
  }
}

save(seurat_objects, file="seurat_objects_start.RData")