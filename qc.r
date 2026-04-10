library(Seurat)

data_dir <- "/home/lab-user/data/scRNAseq_epi2me/seurat_analysis/data_matrices_raw"
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

save(seurat_objects, file="seurat_objects_start.RData")


# Aggiungi metriche QC a ciascun oggetto
for (obj_name in names(seurat_objects)) {
  obj <- seurat_objects[[obj_name]]
  
  # Percentuale geni mitocondriali
  obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^mt-")
  
  # Percentuale geni eritrocitari
  obj[["percent.redcell"]] <- PercentageFeatureSet(
    obj,
    features = c("Hba-a1","Hba-a2","Hbb-bs","Hbb-bt")
  )
  
  # Salva numero di cellule
  obj$sample_name <- obj_name
  
  seurat_objects[[obj_name]] <- obj
}
qc_summary <- data.frame()

for (obj_name in names(seurat_objects)) {
  obj <- seurat_objects[[obj_name]]
  meta <- obj@meta.data
  
  temp <- data.frame(
    sample = obj_name,
    n_cells = ncol(obj),
    median_nFeature = median(meta$nFeature_RNA),
    median_nCount = median(meta$nCount_RNA),
    median_percent_mt = median(meta$percent.mt),
    median_percent_redcell = median(meta$percent.redcell),
    mean_percent_mt = mean(meta$percent.mt),
    mean_percent_redcell = mean(meta$percent.redcell)
  )
  
  qc_summary <- rbind(qc_summary, temp)
}

qc_summary

for (obj_name in names(seurat_objects)) {
  obj <- seurat_objects[[obj_name]]
  
  obj <- subset(
    obj,
    subset = nFeature_RNA > 200 &
             nFeature_RNA < 10000 &
             percent.mt < 10 &
             percent.redcell < 10
  )
  
  seurat_objects[[obj_name]] <- obj
}

png("qc_summary_n_cells_filtered.png", width = 800, height = 600)
ggplot(qc_summary, aes(x = sample, y = n_cells, fill = sample)) +
  geom_col() +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  ) +
  labs(
    title = "Numero di cellule per campione",
    x = "Campione",
    y = "Numero cellule"
  )
dev.off()

qc_summary_2 <- data.frame()

for (obj_name in names(seurat_objects_filtered)) {
  obj <- seurat_objects_filtered[[obj_name]]
  meta <- obj@meta.data
  
  temp <- data.frame(
    sample = obj_name,
    n_cells = ncol(obj),
    median_nFeature = median(meta$nFeature_RNA),
    median_nCount = median(meta$nCount_RNA),
    median_percent_mt = median(meta$percent.mt),
    median_percent_redcell = median(meta$percent.redcell),
    mean_percent_mt = mean(meta$percent.mt),
    mean_percent_redcell = mean(meta$percent.redcell)
  )
  
  qc_summary_2 <- rbind(qc_summary_2, temp)
}

qc_summary_2



png("qc_violin.png", width = 800, height = 600)
for (obj_name in names(seurat_objects)) {
  obj <- seurat_objects[[obj_name]]
  
  p <- VlnPlot(
    obj,
    features = c(
      "nFeature_RNA",
      "nCount_RNA",
      "percent.mt",
      "percent.redcell"
    ),
    ncol = 4,
    pt.size = 0
  ) 
  
  print(p)
}
dev.off()