library(SingleCellExperiment)
library(celda)
library(Seurat)
library(SoupX)
library(dplyr)
library(ggplot2)
library(celda)
library(stringr)

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
# per tuttee le altre, da confermare
# 01 Loading ------------------------------------------------------ #
data_dir <- "/home/lab-user/data/seurat_sicelore_analysis"
files <- list.files(
  data_dir,
  pattern = "matrix_sicelore\\.txt$",
  full.names = TRUE
)
seurat_objects <- list()
for (file in files) {
  dataset_name <- gsub("_matrix_sicelore\\.txt", "", basename(file))
  matrix_data <- read.table(
    file = file,
    header = TRUE,
    row.names = 1,
    sep = "\t",
    check.names = FALSE
  )
  matrix_data <- as.matrix(matrix_data)
  seurat_obj <- CreateSeuratObject(
    counts = matrix_data,
    project = dataset_name
  )
  seurat_obj$condition <- dataset_name
  seurat_objects[[dataset_name]] <- seurat_obj
}
# 02 QC ----------------------------------------------------------- #

for (i in seq_along(seurat_objects)) {
  seurat_objects[[i]]$log10GenesPerUMI <-
    log10(seurat_objects[[i]]$nFeature_RNA) /
    log10(seurat_objects[[i]]$nCount_RNA)
  seurat_objects[[i]]$mitoRatio <-
    PercentageFeatureSet(
      object = seurat_objects[[i]],
      pattern = "^mt-"
    ) / 100
}
# merge
merged_seurat <- merge(
  x = seurat_objects[[1]],
  y = seurat_objects[-1],
  add.cell.ids = names(seurat_objects),
  project = "IntegratedProject"
)
# metadata
metadata <- merged_seurat@meta.data
metadata$cells <- rownames(metadata)

# Rename columns
metadata <- metadata %>%
        dplyr::rename(seq_folder = orig.ident,
                      nUMI = nCount_RNA,
                      nGene = nFeature_RNA)
metadata$sample <- NA
metadata$sample[which(str_detect(metadata$cells, "^adeno"))] <- "adeno"
merged_seurat@meta.data <- metadata
filtered <- subset(
        merged_seurat,
        subset = 
                 nUMI > 250 & 
                 nGene > 300 &
                 log10GenesPerUMI > 0.80 & 
                 mitoRatio < 0.2
    )
# DecontX
seurat_objects_filtered = seurat_objects
sce_list <- lapply(seurat_objects_filtered, as.SingleCellExperiment)
sce_decont <- lapply(sce_list, decontX)
seurat_decont <- lapply(sce_decont, function(x) {
  counts <- decontXcounts(x)
  seu <- CreateSeuratObject(counts = counts)
  return(seu)
})





