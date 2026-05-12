# Exploratory scripts for assessing the quality of the Seurat objects generated from the SiceLore pipeline.
# This is a work in progress, and will be updated as we finish the SiceLore pipeline and generate the Seurat objects for all the samples.

library(Seurat)
library(SoupX)
library(SingleCellExperiment)
library(dplyr)
library(ggplot2)
library(celda)
library(stringr)


# 01 Loading ------------------------------------------------------ #
# ----------------------------------------------------------------- #

load("exploratory_prefilters_qc.RData")
load("exploratory_postfilters_qc.RData")

# 02 format and extract metadata ---------------------------------- #
# ----------------------------------------------------------------- #

# 2a post --------------------------------------------------------- #
post = merged_seurat

metadata <- post@meta.data
metadata$cells <- rownames(metadata)
metadata <- metadata %>%
        dplyr::rename(seq_folder = orig.ident,
                      nUMI = nCount_RNA,
                      nGene = nFeature_RNA)

metadata$sample <- NA
metadata$sample[which(str_detect(metadata$cells, "^adeno"))] <- "adeno"
metadata$sample[which(str_detect(metadata$cells, "^sham"))] <- "sham"
post@meta.data <- metadata


# 02b pre ---------- ---------------------------------------------- #
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
pre = merged_seurat


# 03 plots -------------------------------------------------------- #
# ----------------------------------------------------------------- #


pre_meta  <- pre@meta.data
post_meta <- post@meta.data
pre_meta$filter_status  <- "pre"
post_meta$filter_status <- "post"

common_cols <- intersect(colnames(pre_meta), colnames(post_meta))
metadata_compare <- bind_rows(
  pre_meta[, common_cols],
  post_meta[, common_cols]
)

metadata_compare$filter_status <- factor(
  metadata_compare$filter_status,
  levels = c("pre", "post")
)
# Cell counts
metadata_compare %>%
  ggplot(aes(x = sample, fill = filter_status)) +
  geom_bar(position = "dodge") +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
    plot.title = element_text(hjust = 0.5, face = "bold")
  ) +
  ggtitle("NCells: pre vs post filtering")


# UMI counts
metadata_compare %>%
  ggplot(aes(
    x = nUMI,
    color = filter_status,
    fill = filter_status
  )) +
  geom_density(alpha = 0.2) +
  scale_x_log10() +
  theme_classic() +
  ylab("Cell density") +
  geom_vline(xintercept = 500)

# genes per cell
metadata_compare %>%
  ggplot(aes(
    x = nGene,
    color = filter_status,
    fill = filter_status
  )) +
  geom_density(alpha = 0.2) +
  scale_x_log10() +
  theme_classic() +
  geom_vline(xintercept = 300)


# UMI x Gene per cell
metadata_compare %>%
  ggplot(aes(
    x = mitoRatio,
    color = filter_status,
    fill = filter_status
  )) +
  geom_density(alpha = 0.2) +
  scale_x_log10() +
  theme_classic() +
  geom_vline(xintercept = 0.2)


# Mito ratio
metadata_compare %>%
  ggplot(aes(
    x = nUMI,
    y = nGene,
    color = filter_status
  )) +
  geom_point(alpha = 0.4, size = 0.5) +
  scale_x_log10() +
  scale_y_log10() +
  theme_classic() +
  geom_vline(xintercept = 500) +
  geom_hline(yintercept = 250)