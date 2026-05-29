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
library(scater)


head(annotated@meta.data)

annotated$condition <- ifelse(
  grepl("adeno", merged_seurat$condition),
  "adeno",
  "sham"
)

# Using 'label' and 'sample' as our two factors; each column of the output
# corresponds to one unique combination of these two factors.
summed <- aggregateAcrossCells(
    as.SingleCellExperiment(annotated),
    id = colData(as.SingleCellExperiment(annotated))[, c("celltype", "orig.ident")]
)
label <- "Immune"
current <- summed[,label==summed$celltype]

colnames(current) <- paste(
    current$celltype,
    current$orig.ident,
    sep="_"
)

# Creating up a DGEList object for use in edgeR:
library(edgeR)
y <- DGEList(counts(current), samples=colData(current))
y
discarded <- current$ncells < 10
y <- y[,!discarded]
summary(discarded)

keep <- filterByExpr(y, group=current$condition)
y <- y[keep,]
summary(keep)

y <- calcNormFactors(y)
y$samples

par(mfrow=c(2,3))
for (i in seq_len(ncol(y))) {
    plotMD(y, column=i)
}

plotMDS(cpm(y, log=TRUE), 
    col = ifelse(y$samples$condition == "adeno", "red", "blue"))
y$samples$condition <- factor(y$samples$condition)
y$samples$condition <- relevel(y$samples$condition, ref = "sham")
design <- model.matrix(~ condition, data = y$samples)
design


fit <- glmQLFit(y, design, robust=TRUE)
summary(fit$var.prior)
summary(fit$df.prior)
plotQLDisp(fit)

res <- glmQLFTest(fit, coef=ncol(design))
summary(decideTests(res))

# ------------------------------------------------- #
# Limma-voom approach
library(edgeR)
library(limma)


all.ct <- unique(summed$celltype)
results <- list()
for (ct in all.ct) {
    # subset celltype
    current <- summed[, summed$celltype == ct]
    current <- current[, current$ncells >= 10] # remove low-cell pseudobulks
    if (ncol(current) < 2) next # skip if too few samples
    y <- DGEList(
        counts(current),
        samples = colData(current)
    )
    keep <- filterByExpr(
        y,
        group = y$samples$condition
    )
    y <- y[keep, , keep.lib.sizes=FALSE]
    if (nrow(y) == 0) next
    y <- calcNormFactors(y)
    y$samples$condition <- factor(y$samples$condition)
    y$samples$condition <- relevel(
        y$samples$condition,
        ref = "sham"
    )
    design <- model.matrix(~ condition,data = y$samples)
    v <- voom(y, design)
    fit <- lmFit(v, design)
    fit <- eBayes(fit)
    res <- topTable(fit, coef = "conditionadeno", number = Inf, sort.by = "P")
    results[[ct]] <- res
}


# ------------------------------------------------- #


# INPUT: results (list di topTable per celltype)
de.results <- results

# 1. Unione geni globali
all.genes <- unique(unlist(lapply(de.results, rownames)))


# 2. MATRICE logFC (gene × celltype)
lfc.mat <- sapply(de.results, function(df) {
    res <- rep(NA, length(all.genes))
    names(res) <- all.genes
    res[rownames(df)] <- df$logFC
    return(res)
})

lfc.mat <- as.data.frame(lfc.mat)

# 3. MATRICE P.VALUE (gene × celltype)
p.mat <- sapply(de.results, function(df) {
    res <- rep(NA, length(all.genes))
    names(res) <- all.genes
    res[rownames(df)] <- df$P.Value
    return(res)
})

p.mat <- as.data.frame(p.mat)

# 4. BINARIZZAZIONE DE (p < 0.05)
is.de <- p.mat < 0.05
is.de[is.na(is.de)] <- FALSE


# 5. UP / DOWN CONSISTENCY ACROSS CELL TYPES
# upregulated (logFC > 0 = adeno up)
up.mat <- lfc.mat > 0
up.mat[is.na(up.mat)] <- FALSE

# downregulated (logFC < 0 = adeno down)
down.mat <- lfc.mat < 0
down.mat[is.na(down.mat)] <- FALSE


# 6. TOP CONSISTENT GENES ACROSS CELL TYPES
# genes most consistently up in adeno
top.up <- sort(rowMeans(up.mat), decreasing = TRUE)
head(top.up, 10)

# genes most consistently down in adeno
top.down <- sort(rowMeans(down.mat), decreasing = TRUE)
head(top.down, 10)

# =========================
# 7. CONSISTENTLY SIGNIFICANT GENES
# =========================

sig.consistency <- rowMeans(is.de)
head(sort(sig.consistency, decreasing = TRUE), 10)

# =========================
# 8. CELLTYPE-SPECIFICITY (ESEMPIO: Immune)
# =========================

cell <- "Immune"

# genes not DE in other cell types
not.de <- !is.de
not.de.other <- rowMeans(not.de[, colnames(not.de) != cell, drop = FALSE]) == 1

# DE only in this cell type
unique.degs <- is.de[, cell] & not.de.other
unique.degs <- names(which(unique.degs))

# subset results for that cell type
de.cell <- de.results[[cell]]
de.cell <- de.cell[unique.degs, , drop = FALSE]
de.cell <- de.cell[order(de.cell$P.Value), ]

# show top cell-type-specific genes
head(de.cell)