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