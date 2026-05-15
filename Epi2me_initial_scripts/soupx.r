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

