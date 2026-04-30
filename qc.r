library(Seurat)
library(dplyr)
library(ggplot2)
library(celda)
library(stringr)
merged_seurat = seurat_decounted
View(merged_seurat@meta.data)

merged_seurat$adeno1$log10GenesPerUMI <- log10(merged_seurat$adeno1$nFeature_RNA) / log10(merged_seurat$adeno1$nCount_RNA)
merged_seurat$adeno1$mitoRatio <- PercentageFeatureSet(object = merged_seurat$adeno1, pattern = "^mt-")
merged_seurat$adeno1$mitoRatio <- merged_seurat$adeno1$mitoRatio / 100

merged_seurat$adeno2$log10GenesPerUMI <- log10(merged_seurat$adeno2$nFeature_RNA) / log10(merged_seurat$adeno2$nCount_RNA)
merged_seurat$adeno2$mitoRatio <- PercentageFeatureSet(object = merged_seurat$adeno2, pattern = "^mt-")
merged_seurat$adeno2$mitoRatio <- merged_seurat$adeno2$mitoRatio / 100

merged_seurat$sham1$log10GenesPerUMI <- log10(merged_seurat$sham1$nFeature_RNA) / log10(merged_seurat$sham1$nCount_RNA)
merged_seurat$sham1$mitoRatio <- PercentageFeatureSet(object = merged_seurat$sham1, pattern = "^mt-")
merged_seurat$sham1$mitoRatio <- merged_seurat$sham1$mitoRatio / 100

merged_seurat$sham2$log10GenesPerUMI <- log10(merged_seurat$sham2$nFeature_RNA) / log10(merged_seurat$sham2$nCount_RNA)
merged_seurat$sham2$mitoRatio <- PercentageFeatureSet(object = merged_seurat$sham2, pattern = "^mt-")
merged_seurat$sham2$mitoRatio <- merged_seurat$sham2$mitoRatio / 100

merged_seurat <- merge(merged_seurat$adeno1, 
                       y = merged_seurat[c("adeno2", "sham1", "sham2")], 
                       add.cell.ids = names(merged_seurat[c("adeno1","adeno2", "sham1", "sham2")]), 
                       project = "IntegratedProject")

metadata <- merged_seurat@meta.data
metadata$cells <- rownames(metadata)

# Rename columns
metadata <- metadata %>%
        dplyr::rename(seq_folder = orig.ident,
                      nUMI = nCount_RNA,
                      nGene = nFeature_RNA)


metadata$sample <- NA
metadata$sample[which(str_detect(metadata$cells, "^adeno"))] <- "adeno"
metadata$sample[which(str_detect(metadata$cells, "^sham"))] <- "sham"
merged_seurat@meta.data <- metadata

# Cell counts
metadata %>% 
  	ggplot(aes(x=sample, fill=sample)) + 
  	geom_bar() +
  	theme_classic() +
  	theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  	theme(plot.title = element_text(hjust=0.5, face="bold")) +
  	ggtitle("NCells")

# UMI counts
metadata %>% 
  	ggplot(aes(color=sample, x=nUMI, fill= sample)) + 
  	geom_density(alpha = 0.2) + 
  	scale_x_log10() + 
  	theme_classic() +
  	ylab("Cell density") +
  	geom_vline(xintercept = 500)

metadata %>% 
  	ggplot(aes(color=sample, x=nUMI, fill= ident)) + 
  	geom_density(alpha = 0.2) + 
  	scale_x_log10() + 
  	theme_classic() +
  	ylab("Cell density") +
  	geom_vline(xintercept = 500)

# genes per cell
metadata %>% 
  	ggplot(aes(color=sample, x=nGene, fill= sample)) + 
  	geom_density(alpha = 0.2) + 
  	theme_classic() +
  	scale_x_log10() + 
  	geom_vline(xintercept = 300)

metadata %>% 
  	ggplot(aes(color=sample, x=nGene, fill= ident)) + 
  	geom_density(alpha = 0.2) + 
  	theme_classic() +
  	scale_x_log10() + 
  	geom_vline(xintercept = 300)

# UMI x Gene per cell
metadata %>% 
  	ggplot(aes(x=nUMI, y=nGene, color=mitoRatio)) + 
  	geom_point() + 
	scale_colour_gradient(low = "gray90", high = "black") +
  	stat_smooth(method=lm) +
  	scale_x_log10() + 
  	scale_y_log10() + 
  	theme_classic() +
  	geom_vline(xintercept = 500) +
  	geom_hline(yintercept = 250) +
  	facet_wrap(~sample)

metadata %>% 
  	ggplot(aes(x=nUMI, y=nGene, color=mitoRatio)) + 
  	geom_point() + 
	scale_colour_gradient(low = "gray90", high = "black") +
  	stat_smooth(method=lm) +
  	scale_x_log10() + 
  	scale_y_log10() + 
  	theme_classic() +
  	geom_vline(xintercept = 500) +
  	geom_hline(yintercept = 250) +
  	facet_wrap(~ident)

# Mito ratio
metadata %>% 
  	ggplot(aes(color=sample, x=mitoRatio, fill=sample)) + 
  	geom_density(alpha = 0.2) + 
  	scale_x_log10() + 
  	theme_classic() +
  	geom_vline(xintercept = 0.2)

metadata %>% 
  	ggplot(aes(color=sample, x=mitoRatio, fill=ident)) + 
  	geom_density(alpha = 0.2) + 
  	scale_x_log10() + 
  	theme_classic() +
  	geom_vline(xintercept = 0.2)


filtered <- subset(
        merged_seurat,
        subset = 
                 nUMI > 500 & 
                 nGene > 250 &
                 log10GenesPerUMI > 0.80 & 
                 mitoRatio < 0.2
    )