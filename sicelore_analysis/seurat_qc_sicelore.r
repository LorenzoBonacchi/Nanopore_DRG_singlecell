# Adeno1 run only
library(Seurat)
library(SoupX)
library(SingleCellExperiment)
library(dplyr)
library(ggplot2)
library(celda)
library(stringr)

load("adeno1_decontX.RData")
adeno1$log10GenesPerUMI <- log10(adeno1$nFeature_decontX) / log10(adeno1$nCount_decontX)
adeno1$mitoRatio <- PercentageFeatureSet(object = adeno1, pattern = "^mt-")
adeno1$mitoRatio <- adeno1$mitoRatio / 100

metadata <- adeno1@meta.data
metadata$cells <- rownames(metadata)

# Rename columns
metadata <- metadata %>%
        dplyr::rename(seq_folder = orig.ident,
                      nUMI = nCount_decontX,
                      nGene = nFeature_decontX)
metadata$sample <- NA
metadata$sample[which(str_detect(metadata$cells, "^adeno"))] <- "adeno"
adeno1@meta.data <- metadata
adeno1$sample = "adeno1" #fixed

filtered <- subset(
        adeno1,
        subset = 
                 nUMI > 250 & 
                 nGene > 300 &
                 log10GenesPerUMI > 0.80 & 
                 mitoRatio < 0.2
    )
 
# around 1300 cells pass the filters. Let's see 


save(filtered, adeno1, file = "adeno1_decontX_filtered.RData")