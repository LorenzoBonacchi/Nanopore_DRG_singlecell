# ============================================================== #
# ============================================================== #
# Libraries and data load ====================================== #

library(Seurat)
library(harmony)
library(DoubletFinder)
library(tibble)
library(ggplot2)
setwd("/media/user/8Tb/scRNAseq/seurat_analysis/")
load("seurat_objects_start.RData")


# ============================================================== #
# ============================================================== #
# Scrublet doublet removal ===================================== #

for (object_name in names(seurat_objects)) {
  doublet_score_file <- file.path(
    "scrublet_results",
    paste0(object_name, "_doublet_scores.csv")
  )
  if (!file.exists(doublet_score_file)) {
    warning(paste("File not found:", doublet_score_file, "Skipping..."))
    next
  }
  doublet_scores <- read.csv(doublet_score_file)
  if (nrow(doublet_scores) != ncol(seurat_objects[[object_name]])) {
    stop(paste("Row mismatch for", object_name))
  }
  # Add metadata
  seurat_objects[[object_name]] <- AddMetaData(
    object   = seurat_objects[[object_name]],
    metadata = doublet_scores$DoubletScore,
    col.name = "DoubletScore"
  )
  seurat_objects[[object_name]] <- AddMetaData(
    object   = seurat_objects[[object_name]],
    metadata = as.logical(doublet_scores$PredictedDoublet),
    col.name = "PredictedDoublet"
  )
}


# Loop through each Seurat object and remove doublets
for (object_name in names(seurat_objects)) {
  # Get the current Seurat object
  seu_obj <- seurat_objects[[object_name]]
  seu_obj <- subset(seu_obj, subset = PredictedDoublet == FALSE)
  seurat_objects[[object_name]] <- seu_obj
}

seurat_objects_scrubbed = seurat_objects # as a backup

# ============================================================== #
# ============================================================== #
# Doublet Finder =============================================== #

# DoubletFinder
# Functions 
run_doubletfinder_custom <- function(seu_sample_subset, multiplet_rate = NULL){
  # for debug
  #seu_sample_subset <- samp_split[[1]]
  # Print sample number
  print(paste0("Sample ", unique(seu_sample_subset[['orig.ident']]), '...........')) 
  
  if(is.null(multiplet_rate)){
    print('multiplet_rate not provided....... estimating multiplet rate from cells in dataset')
    
    # 10X multiplet rates table
    #https://rpubs.com/kenneditodd/doublet_finder_example
    multiplet_rates_10x <- data.frame('Multiplet_rate'= c(0.004, 0.008, 0.0160, 0.023, 0.031, 0.039, 0.046, 0.054, 0.061, 0.069, 0.076),
                                      'Loaded_cells' = c(800, 1600, 3200, 4800, 6400, 8000, 9600, 11200, 12800, 14400, 16000),
                                      'Recovered_cells' = c(500, 1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000, 10000))
    
    print(multiplet_rates_10x)
    
    multiplet_rate <- multiplet_rates_10x %>% dplyr::filter(Recovered_cells < nrow(seu_sample_subset@meta.data)) %>% 
      dplyr::slice(which.max(Recovered_cells)) %>% # select the min threshold depending on your number of samples
      dplyr::select(Multiplet_rate) %>% as.numeric(as.character()) # get the expected multiplet rate for that number of recovered cells
    
    print(paste('Setting multiplet rate to', multiplet_rate))
  }
  
  # Pre-process seurat object with standard seurat workflow --- 
  sample <- NormalizeData(seu_sample_subset)
  sample <- FindVariableFeatures(sample)
  sample <- ScaleData(sample)
  sample <- RunPCA(sample, nfeatures.print = 10)
  
  # Find significant PCs
  stdv <- sample[["pca"]]@stdev
  percent_stdv <- (stdv/sum(stdv)) * 100
  cumulative <- cumsum(percent_stdv)
  co1 <- which(cumulative > 90 & percent_stdv < 5)[1] 
  co2 <- sort(which((percent_stdv[1:length(percent_stdv) - 1] - 
                       percent_stdv[2:length(percent_stdv)]) > 0.1), 
              decreasing = T)[1] + 1
  min_pc <- min(co1, co2)
  
  # Finish pre-processing with min_pc
  sample <- RunUMAP(sample, dims = 1:min_pc)
  sample <- FindNeighbors(object = sample, dims = 1:min_pc)              
  sample <- FindClusters(object = sample, resolution = 0.1)
  
  # pK identification (no ground-truth) 
  #introduces artificial doublets in varying props, merges with real data set and 
  # preprocesses the data + calculates the prop of artficial neighrest neighbours, 
  # provides a list of the proportion of artificial nearest neighbours for varying
  # combinations of the pN and pK
  sweep_list <- paramSweep(sample, PCs = 1:min_pc, sct = FALSE)   
  sweep_stats <- summarizeSweep(sweep_list)
  bcmvn <- find.pK(sweep_stats) # computes a metric to find the optimal pK value (max mean variance normalised by modality coefficient)
  # Optimal pK is the max of the bimodality coefficient (BCmvn) distribution
  optimal.pk <- bcmvn %>% 
    dplyr::filter(BCmetric == max(BCmetric)) %>%
    dplyr::select(pK)
  optimal.pk <- as.numeric(as.character(optimal.pk[[1]]))
  
  ## Homotypic doublet proportion estimate
  annotations <- sample@meta.data$seurat_clusters # use the clusters as the user-defined cell types
  homotypic.prop <- modelHomotypic(annotations) # get proportions of homotypic doublets
  
  nExp.poi <- round(multiplet_rate * nrow(sample@meta.data)) # multiply by number of cells to get the number of expected multiplets
  nExp.poi.adj <- round(nExp.poi * (1 - homotypic.prop)) # expected number of doublets
  
  # run DoubletFinder
  sample <- doubletFinder(seu = sample, 
                          PCs = 1:min_pc, 
                          pK = optimal.pk, # the neighborhood size used to compute the number of artificial nearest neighbours
                          nExp = nExp.poi.adj) # number of expected real doublets
  # change name of metadata column with Singlet/Doublet information
  colnames(sample@meta.data)[grepl('DF.classifications.*', colnames(sample@meta.data))] <- "doublet_finder"
  
  # Subset and save
  # head(sample@meta.data['doublet_finder'])
  # singlets <- subset(sample, doublet_finder == "Singlet") # extract only singlets
  # singlets$ident
  double_finder_res <- sample@meta.data['doublet_finder'] # get the metadata column with singlet, doublet info
  double_finder_res <- rownames_to_column(double_finder_res, "row_names") # add the cell IDs as new column to be able to merge correctly
  return(double_finder_res)
}

#### Apply DoubletFinder
doublet_results <- list()
for (name in names(seurat_objects)) {
  seurat_obj <- seurat_objects[[name]]
  doublet_results[[name]] <- run_doubletfinder_custom(seu_sample_subset = seurat_obj)
  print(paste("DoubletFinder completed for", name))
}

#### Remove douplets
seurat_objects_filtered <- list()
for (name in names(seurat_objects)) {
  seurat_obj <- seurat_objects[[name]]
  doublet_res <- doublet_results[[name]]
  
  # Ensure doublet_res is a data frame with proper columns
  if (!is.null(doublet_res)) {
    # Filter the Seurat object to retain only Singlets
    singlet_ids <- doublet_res$row_names[doublet_res$doublet_finder == "Singlet"]
    seurat_objects_filtered[[name]] <- subset(seurat_obj, cells = singlet_ids)
    
    # Print progress
    print(paste("Filtered Seurat object for", name, "- Singlets retained:", length(singlet_ids)))
  } else {
    # If no doublet results, retain the original Seurat object
    seurat_objects_filtered[[name]] <- seurat_obj
    print(paste("No doublet results for", name, "- Original object retained."))
  }
}

# QC filtering after doublets removal
for (obj_name in names(seurat_objects_filtered)) {
    data <- seurat_objects_filtered[[obj_name]]
    # Calculate the percentage of mitochondrial genes
    # Adjust the pattern "^MT-" to match your dataset's mitochondrial gene prefixes if needed
    data[["percent.mt"]] <- PercentageFeatureSet(data, pattern = "^MT-")
    # Calculate the percentage of red cell genes
    data[["percent.redcell"]] <- PercentageFeatureSet(data, features = c("HBA1", "HBA2", "HBB"))
    # Filter based on percent.mt, percent.redcell and standard criteria
    data <- subset(
        data,
        subset = 
                 nFeature_RNA > 200 & 
                 nFeature_RNA < 10000 & 
                 percent.mt < 10 & 
                 percent.redcell < 10
    )
    # Normalize and preprocess the data
    data <- NormalizeData(data)
    data <- FindVariableFeatures(data, selection.method = "vst", nfeatures = 2000)
    
    # Update the filtered object back into the list
    seurat_objects_filtered[[obj_name]] <- data
}

save(seurat_objects_filtered,seurat_objects,file="objects_post_filtering.RData")