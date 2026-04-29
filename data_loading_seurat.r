library(Seurat)
library(SoupX)
data_dir <- "/media/user/8Tb/raw_analysis/raw"
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
# ===================================================== #
# Load processed objects and markers ========================= #
# ===================================================== #
data_dir_proc <- "/media/user/8Tb/raw_analysis/proc"
subdirs <- list.dirs(data_dir_proc, recursive = FALSE)
seurat_objects_proc <- list()
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
    seurat_obj_proc <- CreateSeuratObject(
      counts = matrix_data,
      project = dataset_name
    )
    seurat_obj_proc$condition <- dataset_name #Condition and orig.ident are the same, need to change later for batch reference
    seurat_objects_proc[[dataset_name]] <- seurat_obj_proc

  } else {
    message(paste("Missing files in:", subdir))
  }
}


#save(seurat_objects, file="seurat_objects_start.RData")
clean_objects <- list()

for (name in names(seurat_objects)) {

  raw_obj <- seurat_objects[[name]]
  proc_obj <- seurat_objects_proc[[name]]

  # preprocessing sui processed (non sui raw!)
  proc_obj <- NormalizeData(proc_obj)
  proc_obj <- FindVariableFeatures(proc_obj)
  proc_obj <- ScaleData(proc_obj)
  proc_obj <- RunPCA(proc_obj)
  proc_obj <- FindNeighbors(proc_obj, dims = 1:20)
  proc_obj <- FindClusters(proc_obj, resolution = 0.2)

  # estrai matrici
  raw_counts  <- GetAssayData(raw_obj, assay = "RNA", layer = "counts")
  filt_counts <- GetAssayData(proc_obj, assay = "RNA", layer = "counts")

  # SoupX corretto
  
  sc <- SoupChannel(
    tod = raw_counts,   # tutte le droplets
    toc = filt_counts   # cellule filtrate
  )

  sc <- setClusters(sc, proc_obj$seurat_clusters)
  sc <- autoEstCont(sc)

  clean_counts <- adjustCounts(sc)

  clean_obj <- CreateSeuratObject(counts = clean_counts)
  clean_obj$condition <- name

  clean_objects[[name]] <- clean_obj

  print(paste("SoupX done for", name))
}



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


adeno <- merge(
  x = seurat_objects[[1]],
  y = list(
    seurat_objects[[2]]
  ),
  add.cell.ids = c("adeno1", "adeno2"),
  project = "adeno"
)
sham <- merge(
  x = seurat_objects[[3]],
  y = list(
    seurat_objects[[4]]
  ),
  add.cell.ids = c("sham1", "sham2"),
  project = "sham"
)