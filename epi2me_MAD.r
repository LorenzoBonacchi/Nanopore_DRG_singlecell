
library(Seurat)

# Calcolo percentuale mitocondriale
for (name in names(seurat_objects)) {

  seurat_objects[[name]][["percent.mt"]] <- PercentageFeatureSet(
    seurat_objects[[name]],
    pattern = "^mt-"
  )
}

# Funzione per calcolare mediane e MAD
calculate_mad_thresholds <- function(seurat_obj) {

  # log-transform come nel paper
  log_nCount <- log10(seurat_obj$nCount_RNA + 1)
  log_nFeature <- log10(seurat_obj$nFeature_RNA + 1)
  percent_mt <- seurat_obj$percent.mt

  # mediane
  med_nCount <- median(log_nCount)
  med_nFeature <- median(log_nFeature)
  med_mt <- median(percent_mt)

  # MAD
  mad_nCount <- mad(log_nCount)
  mad_nFeature <- mad(log_nFeature)
  mad_mt <- mad(percent_mt)

  # soglie
  thresholds <- list(

    nCount_lower = med_nCount - 3 * mad_nCount,
    nCount_upper = med_nCount + 3 * mad_nCount,

    nFeature_lower = med_nFeature - 3 * mad_nFeature,
    nFeature_upper = med_nFeature + 3 * mad_nFeature,

    mt_upper = med_mt + 4 * mad_mt
  )

  return(thresholds)
}

# Calcolo per tutti gli oggetti
mad_results <- list()

for (name in names(seurat_objects)) {

  cat("\nProcessing:", name, "\n")

  mad_results[[name]] <- calculate_mad_thresholds(
    seurat_objects[[name]]
  )

  print(mad_results[[name]])
}



for (name in names(seurat_objects)) {

  obj <- seurat_objects[[name]]
  thr <- mad_results[[name]]

  log_nCount <- log10(obj$nCount_RNA + 1)
  log_nFeature <- log10(obj$nFeature_RNA + 1)

  keep_cells <- (
    log_nCount > thr$nCount_lower &
    log_nCount < thr$nCount_upper &
    log_nFeature > thr$nFeature_lower &
    log_nFeature < thr$nFeature_upper &
    obj$percent.mt < thr$mt_upper
  )

  seurat_objects[[name]] <- subset(
    obj,
    cells = colnames(obj)[keep_cells]
  )

  cat(
    name,
    ": kept",
    sum(keep_cells),
    "cells out of",
    length(keep_cells),
    "\n"
  )
}