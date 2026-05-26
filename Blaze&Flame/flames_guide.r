library(rtracklayer)
library(Seurat)
library(DropletUtils)
library(gridExtra)
library(data.table)
library(BiocParallel)
library(celda)
library(SingleCellExperiment)
library(DoubletFinder)
library(stringr)
library(cowplot)
library(grid)
library(patchwork)
library(tidyverse)
library(ORFik)
library(GenomicFeatures)
library(Gviz)
library(BSgenome.Mmusculus.UCSC.mm39)


# Set working directory and create folders for output files
setwd(".")  # Set this to correct location
dir.create("./output_files/ref_files", recursive = TRUE, showWarnings = FALSE)
dir.create("./output_files/counts", recursive = TRUE, showWarnings = FALSE)
dir.create("./output_files/seu_objects", recursive = TRUE, showWarnings = FALSE)
dir.create("./output_files/empty_drops", recursive = TRUE, showWarnings = FALSE)
dir.create("./output_files/decontx", recursive = TRUE, showWarnings = FALSE)
dir.create("./output_files/QC", recursive = TRUE, showWarnings = FALSE)
dir.create("./output_files/DE", recursive = TRUE, showWarnings = FALSE)
dir.create("./output_files/multi_sample", recursive = TRUE, showWarnings = FALSE)

make_isoform_gene_symbol_dict <- function(FLAMES_gtf, 
                                          reference_gtf, 
                                          output_file) {
  # Import the first GTF file (transcripts GTF)
  gtf1 <- import(FLAMES_gtf)
  gtf1_df <- as.data.frame(gtf1)
  
  # Select relevant columns from the first GTF
  selected_columns1 <- gtf1_df[, c("transcript_id", "gene_id")]
  unique_selected_cols <- unique(selected_columns1)
  
  # Import the second GTF file (reference GTF with gene symbols)
  gtf2 <- import(reference_gtf)
  gtf2_df <- as.data.frame(gtf2)
  
  # Select relevant columns from the second GTF
  selected_columns2 <- gtf2_df[, c("gene_name", "gene_id")]
  unique_gene_symbol <- unique(selected_columns2)
  
  # Merge the two data frames on 'gene_id'
  combined_data <- merge(unique_selected_cols, 
                         unique_gene_symbol, 
                         by = "gene_id", 
                         all.x = TRUE)
  
  # If 'gene_name' is missing, replace it with 'gene_id'
  combined_data$gene_symbol <- ifelse(is.na(combined_data$gene_name), 
                                      combined_data$gene_id, 
                                      combined_data$gene_name)
  
  # Select relevant columns
  final_combined_data <- combined_data[, c("transcript_id", "gene_id", "gene_symbol")]
  
  # Write to a CSV file
    write.csv(final_combined_data, file = file.path("output_files/ref_files", output_file), row.names = FALSE)

  
  return(final_combined_data)
}

# The FLAMES ref can be found in your selected output folder after running the Flames pipeline. 
FLAMES_gtf_file <- "/media/user/8Tb/blaze_analysis/blaze_adeno1/ourdir/flames_out/isoform_annotated.gff3" #ensure file is unzipped
reference_gtf_file <- "/media/user/8Tb/blaze_analysis/blaze_adeno1/ourdir/genes.gtf" # ensure file is unzipped
output_file <- "isoform_gene_dict.csv"

# Call the helper function defined in code block above to create a dictionary containing corresponding gene information for each isoform
# This may take a few minutes 
isoform_gene_dict <- make_isoform_gene_symbol_dict(FLAMES_gtf_file,
                                                   reference_gtf_file,
                                                   output_file)



# ----------------------------------------- #
# ----------------------------------------- #
# ----------------------------------------- #


convert_ENSGID_to_geneSymbol <- function(gene_count_matrix_path, 
                                         id_symbol_df = isoform_gene_dict, 
                                         output_file,
                                         return_df = FALSE) {
  
  # Load the reference dictionary we made earlier - select gene-level cols
  id_symbol_df <- as_tibble(id_symbol_df) %>%
    dplyr::select(gene_id, gene_symbol)
  
  # Load the data object with ENSGID row names
  gene_count_matrix <- fread(gene_count_matrix_path, header = TRUE)
  colnames(gene_count_matrix)[1] <- "gene_id"
  
  # Replace ENSGIDs with gene symbols in original flames gene-level count matrix
  formatted_gene_count_matrix <- gene_count_matrix %>%
    merge(id_symbol_df, by.x = 'gene_id', by.y = 'gene_id') %>%   # Add gene symbol information
    distinct(gene_symbol, .keep_all = TRUE) %>%   # Remove duplicates based on gene symbol
    dplyr::select(-gene_id) %>%   # Remove the ENSGID column
    column_to_rownames(var = "gene_symbol")   # use the gene symbols we added as rownames
  
  # Write out the processed data frame
  fwrite(formatted_gene_count_matrix, 
            output_file, 
            row.names = TRUE)
  
  # Return the processed count matrix for further use if needed
  if(return_df){
    return(formatted_gene_count_matrix)
  }
}


# convert Gene_id to gene symbol for gene counts
convert_ENSGID_to_geneSymbol(
  gene_count_matrix_path = "/media/user/8Tb/blaze_analysis/blaze_adeno1/ourdir/flames_out/gene_count.csv",
  output_file = "./output_files/counts/geneSymbol_gene_count.csv"
)

# convert Gene_id to gene symbol for background counts
convert_ENSGID_to_geneSymbol(
  gene_count_matrix_path = "./data/background/gene_count.csv",
  output_file = "./output_files/counts/background_geneSymbol_gene_count.csv"
)





library(Matrix)
library(data.table)
library(Seurat)
library(dplyr)

counts <- readMM("/media/user/8Tb/blaze_analysis/blaze_adeno1/ourdir/flames_out/oarfish.count.mtx")
counts <- t(counts)   # 🔥 FIX CRUCIALE
barcodes <- fread("/media/user/8Tb/blaze_analysis/blaze_adeno1/ourdir/flames_out/barcodes.txt", header = FALSE)
features <- fread("/media/user/8Tb/blaze_analysis/blaze_adeno1/ourdir/flames_out/features.txt", header = FALSE)

colnames(counts) <- barcodes$V1
rownames(counts) <- features$V1   # spesso gene_id o transcript_id

features$gene_id <- sub("_.*", "", features$V1)
rownames(counts) <- features$gene_id

seu <- CreateSeuratObject(counts = counts)

id_symbol_df <- isoform_gene_dict %>%
  dplyr::select(gene_id, gene_symbol) %>%
  distinct()

gene_map <- setNames(id_symbol_df$gene_symbol, id_symbol_df$gene_id)

new_names <- gene_map[rownames(seu)]

# fallback per NA
new_names[is.na(new_names)] <- rownames(seu)[is.na(new_names)]

rownames(seu) <- new_names

seu <- AggregateExpression(seu, group.by = "features", assays = "RNA", slot = "counts")$RNA
seu <- CreateSeuratObject(seu)

head(rownames(seu))




library(Seurat)
library(Matrix)
library(SingleCellExperiment)
library(celda)
library(biomaRt)

data_dir <- "/media/user/8Tb/blaze_analysis/data_analysis/"
subdirs <- list.dirs(data_dir, recursive = FALSE)

seurat_objects <- list()

# -------------------------
# biomaRt ONCE (IMPORTANT optimization)
# -------------------------
mart <- useEnsembl("genes", dataset = "mmusculus_gene_ensembl")

for (subdir in subdirs) {

  matrix_file  <- file.path(subdir, "oarfish.count.mtx")
  feature_file <- file.path(subdir, "oarfish.features.txt")
  barcode_file <- file.path(subdir, "oarfish.barcodes.txt")

  if (file.exists(matrix_file) &
      file.exists(feature_file) &
      file.exists(barcode_file)) {

    dataset_name <- basename(subdir)

    # -------------------------
    # Load matrix
    # -------------------------
    m <- readMM(matrix_file)

    barcodes <- readLines(barcode_file)
    features <- readLines(feature_file)

    # -------------------------
    # orientation fix
    # -------------------------
    if (nrow(m) == length(barcodes) &&
        ncol(m) == length(features)) {
      counts <- t(m)
    } else if (nrow(m) == length(features) &&
               ncol(m) == length(barcodes)) {
      counts <- m
    } else {
      stop(paste("Dimension mismatch in", dataset_name))
    }

    rownames(counts) <- features
    colnames(counts) <- barcodes

    # =====================================================
    # 🔥 1. REGION → GENE ID collapse
    # =====================================================
    gene_ids <- sub("_.*", "", rownames(counts))

    counts_gene <- rowsum(counts, group = gene_ids)

    # =====================================================
    # 🔥 2. GENE ID → SYMBOL mapping (SAFE)
    # =====================================================
    annot <- getBM(
      attributes = c("ensembl_gene_id", "external_gene_name"),
      filters = "ensembl_gene_id",
      values = rownames(counts_gene),
      mart = mart
    )

    map <- annot$external_gene_name
    names(map) <- annot$ensembl_gene_id

    gene_symbols <- map[rownames(counts_gene)]

    # fallback safe
    gene_symbols[is.na(gene_symbols) | gene_symbols == ""] <- rownames(counts_gene)

    rownames(counts_gene) <- gene_symbols

    # =====================================================
    # 🔥 3. FINAL COLLAPSE (handle duplicated symbols)
    # =====================================================
    counts_gene <- rowsum(counts_gene, group = rownames(counts_gene))

    # =====================================================
    # 🔥 4. REMOVE INVALID ROWNAMES (Seurat safety)
    # =====================================================
    counts_gene <- counts_gene[
      rownames(counts_gene) != "" &
      !is.na(rownames(counts_gene)),
    ]

    # =====================================================
    # Create Seurat object
    # =====================================================
    seu <- CreateSeuratObject(
      counts = counts_gene,
      project = dataset_name
    )

    seu$condition <- dataset_name

    # -------------------------
    # QC metrics
    # -------------------------
    seu$log10GenesPerUMI <-
      log10(seu$nFeature_RNA) / log10(seu$nCount_RNA)

    seu$mitoRatio <-
      PercentageFeatureSet(seu, pattern = "^mt-") / 100

    # -------------------------
    # QC filtering
    # (your original logic preserved)
    # -------------------------
    seu <- subset(
      seu,
      subset =
        nCount_RNA > 500 &
        nFeature_RNA > 400 &
        nFeature_RNA < 8000 &
        log10GenesPerUMI > 0.80 &
        mitoRatio < 0.1
    )

    # -------------------------
    # Pre-clustering
    # -------------------------
    seu <- NormalizeData(seu)
    seu <- FindVariableFeatures(seu)
    seu <- ScaleData(seu)
    seu <- RunPCA(seu)
    seu <- FindNeighbors(seu, dims = 1:20)
    seu <- FindClusters(seu, resolution = 0.5)

    # -------------------------
    # decontX
    # -------------------------
    sce <- as.SingleCellExperiment(seu)
    sce <- sce[, colnames(seu)]  # 🔥 IMPORTANT alignment fix
    sce <- decontX(sce,
        z = seu$seurat_clusters
    )

    seu$decontX_contamination <- colData(sce)$decontX_contamination

    decont_mat <- decontXcounts(sce)

    seu[["decontX"]] <- CreateAssayObject(counts = decont_mat)
    # -------------------------
    # store
    # -------------------------
    seurat_objects[[dataset_name]] <- seu

  } else {
    message(paste("Missing files in:", subdir))
  }
}

# -------------------------
# merge robusto
# -------------------------
merged_seurat <- Reduce(
  function(x, y) merge(x, y),
  seurat_objects
)