library(Seurat)
library(SingleCellExperiment)
library(celda)
library(dplyr)
library(stringr)
library(harmony)
library(tibble)
library(presto)
library(celldex)
library(SingleR)


# ---------------------------------------- #

matrix_data <- read.delim(
  "gene_count.csv",
  row.names = 1,
  check.names = FALSE
)

seu <- CreateSeuratObject(
    counts = matrix_data,
    project = "sham1_sup"
    )
seu$condition <- "sham1"
  seu$log10GenesPerUMI <-
    log10(seu$nFeature_RNA) /
    log10(seu$nCount_RNA)
seu$mitoRatio <-
    PercentageFeatureSet(
      seu,
      pattern = "^mt-"
    ) / 100
  # QC filtering
  seu <- subset(
    seu,
    subset =
      nCount_RNA > 500 &
      nFeature_RNA > 200 &
      log10GenesPerUMI > 0.80 &
      mitoRatio < 0.1
  )