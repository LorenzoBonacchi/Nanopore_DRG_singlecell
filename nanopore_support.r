library(Seurat)
library(dplyr)
library(ggplot2)
library(celda)
library(stringr)

# LOAD DATA
# subdirs are the output of Epi2me sc workflow 
# adeno.gene_raw_feature_bc_matrix
# sham.gene_raw_feature_bc_matrix

data_dir <- "/home/lab-user/data/nanopore_support/"
subdirs <- list.dirs(data_dir, recursive = FALSE)
seurat_objects <- list()

for (subdir in subdirs) {

  matrix_file <- file.path(subdir, "matrix.mtx.gz")
  feature_file <- file.path(subdir, "features.tsv.gz")
  barcode_file <- file.path(subdir, "barcodes.tsv.gz")

  if (file.exists(matrix_file) & file.exists(feature_file) & file.exists(barcode_file)) {

    dataset_name <- basename(subdir)

    matrix_data <- ReadMtx(
      mtx = matrix_file,
      cells = barcode_file,
      features = feature_file
    )

    seurat_obj <- CreateSeuratObject(
      counts = matrix_data,
      project = dataset_name
    )
    seurat_obj$condition <- dataset_name #Condition and orig.ident are the same, need to change later for batch reference
    seurat_objects[[dataset_name]] <- seurat_obj

  } else {
    message(paste("Missing files in:", subdir))
  }
}

# QC METRICS AND FILTERING 
seurat_objects$adeno.gene_raw_feature_bc_matrix$log10GenesPerUMI <- log10(seurat_objects$adeno.gene_raw_feature_bc_matrix$nFeature_RNA) / log10(seurat_objects$adeno.gene_raw_feature_bc_matrix$nCount_RNA)
seurat_objects$adeno.gene_raw_feature_bc_matrix$mitoRatio <- PercentageFeatureSet(object = seurat_objects$adeno.gene_raw_feature_bc_matrix, pattern = "^mt-")
seurat_objects$adeno.gene_raw_feature_bc_matrix$mitoRatio <- seurat_objects$adeno.gene_raw_feature_bc_matrix$mitoRatio / 100

seurat_objects$sham.gene_raw_feature_bc_matrix$log10GenesPerUMI <- log10(seurat_objects$sham.gene_raw_feature_bc_matrix$nFeature_RNA) / log10(seurat_objects$sham.gene_raw_feature_bc_matrix$nCount_RNA)
seurat_objects$sham.gene_raw_feature_bc_matrix$mitoRatio <- PercentageFeatureSet(object = seurat_objects$sham.gene_raw_feature_bc_matrix, pattern = "^mt-")
seurat_objects$sham.gene_raw_feature_bc_matrix$mitoRatio <- seurat_objects$sham.gene_raw_feature_bc_matrix$mitoRatio / 100

# MERGE AND FORMAT METADATA
merged_seurat <- merge(
  x = seurat_objects[[1]],
  y = list(
    seurat_objects[[2]]
  ),
  add.cell.ids = c("adeno", "sham"),
  project = "IntegratedProject"
)

metadata <- merged_seurat@meta.data
metadata$cells <- rownames(metadata)
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

# genes per cell
metadata %>% 
  	ggplot(aes(color=sample, x=nGene, fill= sample)) + 
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


# Mito ratio
metadata %>% 
  	ggplot(aes(color=sample, x=mitoRatio, fill=sample)) + 
  	geom_density(alpha = 0.2) + 
  	scale_x_log10() + 
  	theme_classic() +
  	geom_vline(xintercept = 0.2)


# FINAL FILTER ACCORDING TO "GENERAL THRESHOLDS"

filtered <- subset(
        merged_seurat,
        subset = 
                 nUMI > 500 & 
                 nGene > 250 &
                 log10GenesPerUMI > 0.80 & 
                 mitoRatio < 0.2
    )

# merged seurat 112162 | adeno 57129  sham 55033 
# filtered 38495 |  adeno 611 sham 37884 
