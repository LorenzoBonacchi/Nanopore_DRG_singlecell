library(Seurat)
library(SoupX)


data_dir_raw <- "/media/user/8Tb/raw_analysis/raw"
subdirs_raw <- list.dirs(data_dir_raw, recursive = FALSE)

raw_objects <- list()
for (subdir in subdirs_raw) {

  matrix_file  <- file.path(subdir, "matrix.mtx.gz")
  features     <- file.path(subdir, "features.tsv.gz")
  barcodes     <- file.path(subdir, "barcodes.tsv.gz")

  if (file.exists(matrix_file)) {
    sample <- basename(subdir)
    counts <- ReadMtx(
      mtx = matrix_file,
      cells = barcodes,
      features = features
    )
    # RAW NON filtrato → SoupX TOD
    raw_objects[[sample]] <- counts
  }
}

data_dir_proc <- "/media/user/8Tb/raw_analysis/proc"
subdirs_proc <- list.dirs(data_dir_proc, recursive = FALSE)

proc_objects <- list()

for (subdir in subdirs_proc) {

  matrix_file  <- file.path(subdir, "matrix.mtx.gz")
  features     <- file.path(subdir, "features.tsv.gz")
  barcodes     <- file.path(subdir, "barcodes.tsv.gz")

  if (file.exists(matrix_file)) {
    sample <- basename(subdir)
    counts <- ReadMtx(
      mtx = matrix_file,
      cells = barcodes,
      features = features
    )
    proc_objects[[sample]] <- counts
  }
}

stopifnot(identical(names(raw_objects), names(proc_objects)))

clean_objects <- list()

for (name in names(raw_objects)) {

  cat("\nProcessing:", name, "\n")

  raw_counts  <- raw_objects[[name]]
  filt_counts <- proc_objects[[name]]

  # -------------------------
  # cell intersection NON serve qui
  # SoupX lo gestisce meglio così
  # -------------------------
  common_genes <- intersect(rownames(raw_counts), rownames(filt_counts))
  raw_counts  <- raw_counts[common_genes, ]
  filt_counts <- filt_counts[common_genes, ]
  # -------------------------
  # SoupX object
  # -------------------------
  sc <- SoupChannel(
    tod = raw_counts,
    toc = filt_counts
  )
  # clustering NON obbligatorio in pipeline base
  sc <- setClusters(sc, rep(1, ncol(filt_counts)))
  # -------------------------
  # robust estimation (IMPORTANT FIX)
  # -------------------------
  sc <- setContaminationFraction(sc, 0.05)
  clean_counts <- adjustCounts(sc)
  # -------------------------
  # Seurat output
  # -------------------------
  obj <- CreateSeuratObject(clean_counts)
  obj$condition <- name
  clean_objects[[name]] <- obj
  cat("SoupX done:", name, "\n")
}


library(Seurat)
library(SoupX)
library(dplyr)

# =====================================================
# PATHS
# =====================================================
data_dir_raw  <- "/home/lab-user/data/Epi2me_raw_analysis/raw"
data_dir_proc <- "/home/lab-user/data/Epi2me_raw_analysis/proc"

raw_dirs  <- list.dirs(data_dir_raw,  recursive = FALSE)
proc_dirs <- list.dirs(data_dir_proc, recursive = FALSE)

seurat_objects <- list()
seurat_objects_proc <- list()
clean_objects <- list()

# =====================================================
# 1. LOAD RAW + PROCESSED (NO NORMALIZATION)
# =====================================================

for (dir in raw_dirs) {
  matrix_file  <- file.path(dir, "matrix.mtx.gz")
  feature_file <- file.path(dir, "features.tsv.gz")
  barcode_file <- file.path(dir, "barcodes.tsv.gz")
  if (file.exists(matrix_file) & file.exists(feature_file) & file.exists(barcode_file)) {
    name <- basename(dir)
    counts <- ReadMtx(
      mtx = matrix_file,
      cells = barcode_file,
      features = feature_file
    )
    seurat_objects[[name]] <- CreateSeuratObject(counts = counts)
  }
}

for (dir in proc_dirs) {

  matrix_file  <- file.path(dir, "matrix.mtx.gz")
  feature_file <- file.path(dir, "features.tsv.gz")
  barcode_file <- file.path(dir, "barcodes.tsv.gz")

  if (file.exists(matrix_file) & file.exists(feature_file) & file.exists(barcode_file)) {
    name <- basename(dir)
    counts <- ReadMtx(
      mtx = matrix_file,
      cells = barcode_file,
      features = feature_file
    )

    seurat_objects_proc[[name]] <- CreateSeuratObject(counts = counts)
  }
}

# =====================================================
# 2. SOUPX LOOP (CORRECTED)
# =====================================================

for (name in names(seurat_objects)) {
  message("Processing: ", name)
  raw_obj  <- seurat_objects[[name]]
  proc_obj <- seurat_objects_proc[[name]]

  # -----------------------------
  # RAW COUNTS (TOD)
  # -----------------------------
  raw_counts <- GetAssayData(raw_obj, layer = "counts")

  # -----------------------------
  # FILTERED COUNTS (TOC)
  # -----------------------------
  filt_counts <- GetAssayData(proc_obj, layer = "counts")

  # -----------------------------
  # ALIGN GENES
  # -----------------------------
  common_genes <- intersect(rownames(raw_counts), rownames(filt_counts))

  raw_counts  <- raw_counts[common_genes, ]
  filt_counts <- filt_counts[common_genes, ]

  filt_counts <- filt_counts[rownames(raw_counts), ]

  # -----------------------------
  # CREATE SOUP CHANNEL
  # -----------------------------
  sc <- SoupChannel(
    tod = raw_counts,
    toc = filt_counts
  )

  # =====================================================
  # IMPORTANT: CLUSTERING TEMP OBJECT (FIX)
  # =====================================================
  tmp <- CreateSeuratObject(counts = filt_counts)
  tmp <- NormalizeData(tmp)
  tmp <- FindVariableFeatures(tmp)
  tmp <- ScaleData(tmp)
  tmp <- RunPCA(tmp)
  tmp <- FindNeighbors(tmp, dims = 1:20)
  tmp <- FindClusters(tmp, resolution = 0.8)

  sc <- setClusters(sc, tmp$seurat_clusters)

  # -----------------------------
  # ESTIMATE SOUP
  # -----------------------------
sc <- autoEstCont(sc, 
                  forceAccept = TRUE,
                  tfidfMin = 0.5,
                  soupQuantile = 0.25)

  rho <- sc$rho
  message("Estimated contamination (rho): ", rho)

  # -----------------------------
  # CORRECT COUNTS
  # -----------------------------
  clean_counts <- adjustCounts(sc)

  # -----------------------------
  # FINAL SEURAT OBJECT
  # -----------------------------
  clean_obj <- CreateSeuratObject(counts = clean_counts)
  clean_obj$condition <- name

  # ONLY NOW FULL PIPELINE
  clean_obj <- NormalizeData(clean_obj)
  clean_obj <- FindVariableFeatures(clean_obj)
  clean_obj <- ScaleData(clean_obj)
  clean_obj <- RunPCA(clean_obj)
  clean_obj <- FindNeighbors(clean_obj, dims = 1:20)
  clean_obj <- FindClusters(clean_obj, resolution = 0.2)

  clean_objects[[name]] <- clean_obj

  message("DONE: ", name)
}


