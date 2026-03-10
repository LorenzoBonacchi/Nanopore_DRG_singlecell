# ============================================================== #
# ============================================================== #
# Libraries and data load ====================================== #

library(Seurat)
library(harmony)
library(DoubletFinder)
library(tibble)
library(ggplot2)
setwd("/media/user/8Tb/scRNAseq/seurat_analysis/")
load("seurat_objects_start.RData")


# ============================================================== #
# ============================================================== #
# Scrublet doublet removal ===================================== #

for (object_name in names(seurat_objects)) {
  doublet_score_file <- file.path(
    "scrublet_results",
    paste0(object_name, "_doublet_scores.csv")
  )
  if (!file.exists(doublet_score_file)) {
    warning(paste("File not found:", doublet_score_file, "Skipping..."))
    next
  }
  doublet_scores <- read.csv(doublet_score_file)
  if (nrow(doublet_scores) != ncol(seurat_objects[[object_name]])) {
    stop(paste("Row mismatch for", object_name))
  }
  # Add metadata
  seurat_objects[[object_name]] <- AddMetaData(
    object   = seurat_objects[[object_name]],
    metadata = doublet_scores$DoubletScore,
    col.name = "DoubletScore"
  )
  seurat_objects[[object_name]] <- AddMetaData(
    object   = seurat_objects[[object_name]],
    metadata = as.logical(doublet_scores$PredictedDoublet),
    col.name = "PredictedDoublet"
  )
}


# Loop through each Seurat object and remove doublets
for (object_name in names(seurat_objects)) {
  # Get the current Seurat object
  seu_obj <- seurat_objects[[object_name]]
  seu_obj <- subset(seu_obj, subset = PredictedDoublet == FALSE)
  seurat_objects[[object_name]] <- seu_obj
}

seurat_objects_scrubbed = seurat_objects # as a backup

# ============================================================== #
# ============================================================== #
# Doublet Finder =============================================== #

