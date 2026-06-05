# --------------------------------------------------- #
# LIBRARIES
# --------------------------------------------------- #

library(Seurat)
library(Matrix)
library(rtracklayer)
library(SingleCellExperiment)
library(celda)
library(harmony)
library(DoubletFinder)
library(dplyr)
library(tibble)
library(SoupX)
library(harmony)

# --------------------------------------------------- #
# PATHS
# --------------------------------------------------- #

data_dir <- "/media/user/8Tb/blaze_analysis/data_analysis/"
gtf_file <- file.path(data_dir, "genes.gtf")
subdirs <- list.dirs(data_dir, recursive = FALSE)
seurat_objects <- list()

# --------------------------------------------------- #
# LOAD GTF AND CREATE ENSG -> SYMBOL MAP
# --------------------------------------------------- #

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

# --------------------------------------------------- #
# DOUBLETFINDER FUNCTION (unchanged)
# --------------------------------------------------- #

run_doubletfinder_custom <- function(
    seu_sample_subset,
    multiplet_rate = NULL
){

  message("Running DoubletFinder on: ",
          unique(seu_sample_subset$condition))

  if (is.null(multiplet_rate)) {

    multiplet_rates_10x <- data.frame(
      Multiplet_rate = c(0.004,0.008,0.0160,0.023,0.031,0.039,0.046,0.054,0.061,0.069,0.076),
      Recovered_cells = c(500,1000,2000,3000,4000,5000,6000,7000,8000,9000,10000)
    )

    multiplet_rate <- multiplet_rates_10x %>%
      filter(Recovered_cells < ncol(seu_sample_subset)) %>%
      slice(which.max(Recovered_cells)) %>%
      pull(Multiplet_rate)

    if (length(multiplet_rate) == 0 || is.na(multiplet_rate)) {
      multiplet_rate <- 0.01
    }
  }

  stdv <- seu_sample_subset[["pca"]]@stdev
  percent_stdv <- (stdv / sum(stdv)) * 100
  cumulative <- cumsum(percent_stdv)

  co1 <- which(cumulative > 90 & percent_stdv < 5)[1]
  diff_drop <- c(0, diff(percent_stdv))
  co2 <- which(diff_drop < -0.1)[1]

  min_pc <- min(co1, co2, na.rm = TRUE)
  if (is.na(min_pc) || min_pc < 10) min_pc <- 10

  sweep_list <- paramSweep(seu_sample_subset, PCs = 1:min_pc, sct = FALSE)
  sweep_stats <- summarizeSweep(sweep_list)
  bcmvn <- find.pK(sweep_stats)

  bcmvn <- as.data.frame(bcmvn)
  bcmvn$BCmetric <- suppressWarnings(as.numeric(as.character(bcmvn$BCmetric)))
  bcmvn <- bcmvn[!is.na(bcmvn$BCmetric), ]

  optimal.pk <- as.numeric(as.character(bcmvn$pK[which.max(bcmvn$BCmetric)]))

  annotations <- as.character(seu_sample_subset$seurat_clusters)
  homotypic.prop <- modelHomotypic(annotations)

  nExp.poi <- round(multiplet_rate * ncol(seu_sample_subset))
  nExp.poi.adj <- round(nExp.poi * (1 - homotypic.prop))

  seu_sample_subset <- doubletFinder(
    seu = seu_sample_subset,
    PCs = 1:min_pc,
    pN = 0.25,
    pK = optimal.pk,
    nExp = nExp.poi.adj,
    reuse.pANN = NULL,
    sct = FALSE
  )

  df_col <- grep("DF.classifications",
                 colnames(seu_sample_subset@meta.data),
                 value = TRUE)

  colnames(seu_sample_subset@meta.data)[
    colnames(seu_sample_subset@meta.data) == df_col
  ] <- "doublet_finder"

  return(seu_sample_subset)
}

