library(Seurat)
library(harmony)
library(DoubletFinder)
library(tibble)
library(celldex)
library(SingleR)
library(ggplot2)
library(tidyverse)
library(dplyr)
library(DESeq2)
library(limma)
library(Matrix)
library(purrr)

head(annotated@meta.data)

# counts matrix
dim(GetAssayData(annotated, layer = "counts"))
GetAssayData(annotated, layer = "counts")[1:6,1:6]

kids <- purrr::set_names(levels(annotated$celltype))
kids
# Total number of clusters
nk <- length(kids)
nk
# Named vector of sample names
annotated$condition <- factor(annotated$condition)
sids <- purrr::set_names(levels(annotated$condition))
sids
# Total number of samples 
ns <- length(sids)
ns

## Determine the number of cells per sample
table(annotated$orig.ident)

## Numeric vector of cells per sample
n_cells <- as.numeric(table(annotated$orig.ident))

## Sample IDs
sids <- names(table(annotated$orig.ident))

## Match sample order
m <- match(sids, unique(annotated@meta.data$orig.ident))

annotated$group <- ifelse(
  grepl("adeno", annotated$orig.ident),
  "Adeno",
  "Sham"
)

ei <- annotated@meta.data %>%
  distinct(orig.ident, group) %>%
  mutate(
    n_cells = as.numeric(table(annotated$orig.ident))
  )

ei

# ------------------------------------------------ #
# Count aggregation
# Aggregate the counts per sample_id and cluster_id

cts <- LayerData(annotated, assay = "RNA", layer = "counts")
groups <- annotated@meta.data[, c("celltype", "orig.ident")]
group_index <- interaction(groups$celltype, groups$orig.ident)
pb <- rowsum(as.matrix(t(cts)), group_index)


pb <- t(pb)


meta_columns <- c("orig.ident", "condition")
meta <- annotated@meta.data %>%
            select(meta_columns) %>%
            unique() %>%
            remove_rownames()

meta
bulk <- AggregateExpression(
            annotated,
            return.seurat = T,
            assays = "RNA",
            group.by = c("seurat_clusters", "orig.ident", "condition")
)

n_cells <- annotated@meta.data %>% 
              dplyr::count(orig.ident, seurat_clusters) 
            
n_cells$orig.ident <- str_replace(n_cells$orig.ident, "_", "-")

meta_bulk <- left_join(bulk@meta.data, n_cells)
rownames(meta_bulk) <- meta_bulk$orig.ident

bulk@meta.data$orig.ident <- stringr::str_extract(bulk@meta.data$orig.ident, "BP\\d+-skin|HC\\d+-skin")
meta_bulk <- left_join(bulk@meta.data, n_cells, by = c("orig.ident", "seurat_clusters"))
bulk@meta.data <- meta_bulk

# Turn condition into a factor
bulk$condition <- factor(bulk$condition, levels = c("BP", "HC"))

bulk@meta.data %>% head()
bulk_kera <- subset(bulk, subset= (celltype == "keratinocytes") & (condition %in% c("BP", "HC")))

png("counts_deg_kera.png",width=720,height=720)
ggplot(bulk_kera@meta.data, aes(x=sample, y=n.y, fill=condition)) +
    geom_bar(stat="identity", color="black") +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
    labs(x="Sample name", y="Number of cells") +
    geom_text(aes(label=n.y), vjust=-0.5)
dev.off()


cluster_counts <- FetchData(bulk_kera, layer="counts", vars=rownames(bulk_kera))
dds <- DESeqDataSetFromMatrix(t(cluster_counts),
                                colData = bulk_kera@meta.data,
                                design = ~ condition)


dds
dds$condition <- relevel(dds$condition, ref = "control")
dds <- DESeq(dds)  # Run the analysis
results <- results(dds)
results = na.omit(results)
filtered_df <- results[!is.na(results$padj) & results$padj < 0.1, ]
filtered_df <- filtered_df[filtered_df$log2FoldChange > 0.58 | filtered_df$log2FoldChange < -0.58,]