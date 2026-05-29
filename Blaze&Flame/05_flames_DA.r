
library(edgeR)

# ------------------------------------------------- #
# Abundances matrix

abundances <- table(annotated$celltype, annotated$orig.ident)
abundances <- unclass(abundances)

meta <- data.frame(orig.ident = annotated$orig.ident, condition = annotated$condition)
meta <- meta[!duplicated(meta$orig.ident), ]
meta <- meta[match(colnames(abundances), meta$orig.ident), ]
rownames(meta) <- meta$orig.ident

# ------------------------------------------------- #
# DGEList and filtering
y.ab <- DGEList(abundances, samples = meta)
keep <- filterByExpr(y.ab, group = y.ab$samples$condition)
y.ab <- y.ab[keep, , keep.lib.sizes = FALSE]
y.ab <- calcNormFactors(y.ab)

# ------------------------------------------------- #
# Design (adeno vs sham)
y.ab$samples$condition <- factor(y.ab$samples$condition)
y.ab$samples$condition <- relevel(y.ab$samples$condition, ref = "sham")
design <- model.matrix(~ condition, data = y.ab$samples)
y.ab <- estimateDisp(y.ab, design)
fit.ab <- glmQLFit(y.ab, design)

# ------------------------------------------------- #
# Test
res <- glmQLFTest(fit.ab, coef = "conditionadeno")
summary(decideTests(res))
topTags(res)