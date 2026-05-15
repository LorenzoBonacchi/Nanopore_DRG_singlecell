library(SingleCellExperiment)
library(celda)
library(Seurat)
library(SoupX)

data_dir <- "/home/lab-user/data/Epi2me_raw_analysis/raw"
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
seurat_objects_filtered = seurat_objects
sce_list <- lapply(seurat_objects_filtered, as.SingleCellExperiment)
sce_decont <- lapply(sce_list, decontX)
seurat_decont <- lapply(sce_decont, function(x) {
  counts <- decontXcounts(x)
  
  seu <- CreateSeuratObject(counts = counts)
  return(seu)
})

quantile(colData(sce_decont[[1]])$decontX_contamination,
         probs = c(0.5, 0.75, 0.9, 0.95))
quantile(colData(sce_decont[[2]])$decontX_contamination,
         probs = c(0.5, 0.75, 0.9, 0.95))
quantile(colData(sce_decont[[3]])$decontX_contamination,
         probs = c(0.5, 0.75, 0.9, 0.95))
quantile(colData(sce_decont[[4]])$decontX_contamination,
         probs = c(0.5, 0.75, 0.9, 0.95))

thresholds <- c(
  adeno1 = 0.05,
  adeno2 = 0.05,
  sham1 = 0.1,
  sham2 = 0.2
)

sce_filtered <- mapply(function(x, name) {
  thr <- thresholds[[name]]
  x[, x$decontX_contamination < thr]
}, sce_decont, names(sce_decont), SIMPLIFY = FALSE)

# per tornare su Seurat (da cofnermare)
sce_to_seurat <- function(sce) {
    counts <- celda::decontXcounts(sce)
    seu <- Seurat::CreateSeuratObject(
        counts = counts,
        project = "decontX",
        min.cells = 0,
        min.features = 0
    )
   meta <- as.data.frame(SummarizedExperiment::colData(sce))
   seu <- Seurat::AddMetaData(seu, metadata = meta)
   return(seu)
}

seurat_decounted <- lapply(sce_filtered, sce_to_seurat)