library(SingleCellExperiment)
library(celda)
library(Seurat)
library(SoupX)

data_dir <- "/media/user/8Tb/raw_analysis/raw"
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




# DecontX
sce_list <- lapply(seurat_objects_filtered, as.SingleCellExperiment)
sce_decont <- lapply(sce_list, decontX)
seurat_decont <- lapply(sce_decont, function(x) {
  counts <- decontXcounts(x)
  
  seu <- CreateSeuratObject(counts = counts)
  return(seu)
})

sce_filtered <- lapply(sce_decont, function(x) {
  x <- x[, x$decontX_contamination < 0.2]
  return(x)
})


# per tornare su Seurat (da cofnermare)
sce_to_seurat <- function(sce) {
  
  counts <- SummarizedExperiment::assay(sce, "decontXcounts")
  
  seu <- Seurat::CreateSeuratObject(
    counts = counts,
    project = "decontX",
    min.cells = 0,
    min.features = 0
  )
  
  seu <- Seurat::AddMetaData(seu, colData(sce))
  
  return(seu)
}

