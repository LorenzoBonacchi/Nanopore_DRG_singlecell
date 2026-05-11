library(SingleCellExperiment)
library(celda)
library(Seurat)
library(SoupX)

###################################################################
# adeno 1 run only
mat <- read.table(
  "adeno1_matrix_sicelore.txt",
  header = TRUE,
  row.names = 1,
  sep = "\t",
  check.names = FALSE
)
mat <- as.matrix(mat)
obj <- CreateSeuratObject(
  counts = mat,
  project = "adeno1run1",
  min.cells = 3,
  min.features = 200
)
obj

sce <- as.SingleCellExperiment(obj)
sce <- decontX(sce)
decont_counts <- decontXcounts(sce)
obj[["decontX"]] <- CreateAssayObject(counts = decont_counts)
DefaultAssay(obj) <- "decontX"
adeno1 = obj
save(adeno1, file = "adeno1_decontX.RData")
#####################################################################






# DecontX
seurat_objects_filtered = obj
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