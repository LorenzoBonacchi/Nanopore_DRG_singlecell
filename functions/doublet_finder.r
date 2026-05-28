

run_doubletfinder_custom <- function(
    seu_sample_subset,
    multiplet_rate = NULL
){

  message("Running DoubletFinder on: ",
          unique(seu_sample_subset$condition))

  # ================================================= #
  # 1. Estimate multiplet rate
  # ================================================= #

  if (is.null(multiplet_rate)) {

    multiplet_rates_10x <- data.frame(
      Multiplet_rate = c(
        0.004, 0.008, 0.0160, 0.023,
        0.031, 0.039, 0.046, 0.054,
        0.061, 0.069, 0.076
      ),
      Recovered_cells = c(
        500, 1000, 2000, 3000,
        4000, 5000, 6000, 7000,
        8000, 9000, 10000
      )
    )

    multiplet_rate <- multiplet_rates_10x %>%
      dplyr::filter(Recovered_cells < ncol(seu_sample_subset)) %>%
      dplyr::slice(which.max(Recovered_cells)) %>%
      dplyr::pull(Multiplet_rate)

    if (length(multiplet_rate) == 0 || is.na(multiplet_rate)) {
      multiplet_rate <- 0.01
    }

    message("Estimated multiplet rate: ", multiplet_rate)
  }

  # ================================================= #
  # 2. Determine PCs safely
  # ================================================= #

  stdv <- seu_sample_subset[["pca"]]@stdev
  percent_stdv <- (stdv / sum(stdv)) * 100
  cumulative <- cumsum(percent_stdv)

  co1 <- which(cumulative > 90 & percent_stdv < 5)[1]

  diff_drop <- c(0, diff(percent_stdv))
  co2 <- which(diff_drop < -0.1)[1]

  min_pc <- min(co1, co2, na.rm = TRUE)

  if (is.na(min_pc) || min_pc < 10) {
    min_pc <- 10
  }

  message("Using PCs: 1:", min_pc)

  # ================================================= #
  # 3. Parameter sweep
  # ================================================= #

  sweep_list <- paramSweep(
    seu_sample_subset,
    PCs = 1:min_pc,
    sct = FALSE
  )

  sweep_stats <- summarizeSweep(sweep_list)
  bcmvn <- find.pK(sweep_stats)

  # ================================================= #
  # 4. FIX pK selection (robust)
  # ================================================= #

  bcmvn <- as.data.frame(bcmvn)

  bcmvn$BCmetric <- suppressWarnings(
    as.numeric(as.character(bcmvn$BCmetric))
  )

  bcmvn <- bcmvn[!is.na(bcmvn$BCmetric), ]

  optimal.pk <- bcmvn$pK[which.max(bcmvn$BCmetric)]
  optimal.pk <- as.numeric(as.character(optimal.pk))

  message("Optimal pK: ", optimal.pk)

  # ================================================= #
  # 5. FIX: safe cluster handling (IMPORTANT)
  # ================================================= #

  annotations <- seu_sample_subset$seurat_clusters
  annotations <- as.character(annotations)

  homotypic.prop <- modelHomotypic(annotations)

  nExp.poi <- round(multiplet_rate * ncol(seu_sample_subset))
  nExp.poi.adj <- round(nExp.poi * (1 - homotypic.prop))

  message("Expected doublets: ", nExp.poi.adj)

  # ================================================= #
  # 6. Clean metadata (prevents xtfrm.data.frame)
  # ================================================= #

  bad_cols <- sapply(seu_sample_subset@meta.data, is.data.frame)
  if (any(bad_cols)) {
    seu_sample_subset@meta.data <- seu_sample_subset@meta.data[, !bad_cols, drop = FALSE]
  }

  # ================================================= #
  # 7. Ensure PCA is valid (critical for DF stability)
  # ================================================= #

  if (any(is.na(seu_sample_subset[["pca"]]@stdev))) {
    seu_sample_subset <- NormalizeData(seu_sample_subset)
    seu_sample_subset <- FindVariableFeatures(seu_sample_subset)
    seu_sample_subset <- ScaleData(seu_sample_subset)
    seu_sample_subset <- RunPCA(seu_sample_subset, verbose = FALSE)
  }

  # ================================================= #
  # 8. Run DoubletFinder
  # ================================================= #

  seu_sample_subset <- doubletFinder(
    seu = seu_sample_subset,
    PCs = 1:min_pc,
    pN = 0.25,
    pK = optimal.pk,
    nExp = nExp.poi.adj,
    reuse.pANN = FALSE,
    sct = FALSE
  )

  # ================================================= #
  # 9. Rename output column safely
  # ================================================= #

  df_col <- grep(
    "DF.classifications",
    colnames(seu_sample_subset@meta.data),
    value = TRUE
  )

  if (length(df_col) > 0) {
    colnames(seu_sample_subset@meta.data)[
      colnames(seu_sample_subset@meta.data) == df_col
    ] <- "doublet_finder"
  }

  return(seu_sample_subset)
}