library(Seurat)
library(harmony)
library(dplyr)



DefaultAssay(merged_seurat) <- "RNA"
merged_seurat <- NormalizeData(merged_seurat)
merged_seurat <- FindVariableFeatures(
  merged_seurat,
  selection.method = "vst",
  nfeatures = 2000
)

var_genes <- VariableFeatures(merged_seurat)
var_genes <- var_genes[
  !grepl("^mt-|^Mt-", var_genes)
]
VariableFeatures(merged_seurat) <- var_genes

merged_seurat <- ScaleData(
  merged_seurat,
  features = VariableFeatures(merged_seurat)
)
merged_seurat <- RunPCA(
  merged_seurat,
  features = VariableFeatures(merged_seurat)
)
merged_seurat <- RunHarmony(
  merged_seurat,
  group.by.vars = "condition",
)
merged_seurat <- RunUMAP(
  merged_seurat,
  reduction = "harmony",
  dims = 1:30
)
merged_seurat <- FindNeighbors(
  merged_seurat,
  reduction = "harmony",
  dims = 1:30
)
merged_seurat <- FindClusters(
  merged_seurat,
  resolution = 0.5
)

# Find All Markers -------------------- #
merged_seurat = JoinLayers(merged_seurat)
markers <- FindAllMarkers(
  merged_seurat,
  only.pos = TRUE,
  test.use = "wilcox",
  logfc.threshold = 0.25,
  min.pct = 0.25
)
markers_top <- markers %>%
  filter(p_val_adj < 0.05,
         pct.1 > 0.25,
         avg_log2FC > 0.25)
markers_top <- markers_top %>%
  group_by(cluster) %>%
  slice_max(avg_log2FC, n = 20, with_ties = FALSE) %>%
  ungroup()



# Dotplots ---------------------------- #
neurons = c("Snap25", "Rbfox3", "Pvalb")
sgc = c("Fabp7", "Ednrb")
myelinating_schwann = c("Mpz", "Mbp")
nonmyelinating_schwann = c("Scn7a")
immune = c("Ptprc", "Ccr2")
endothelial = c("Pecam1", "Flt1")
pericytes = c("Notch3", "Kcnj8")
fibroblasts = c("Pdgfra", "Tbx18")
all = c(neurons, sgc, myelinating_schwann, nonmyelinating_schwann, immune, endothelial, pericytes, fibroblasts)




# Wang annotation
wang_neurons = c("Snap25","Syt1","Avil","Gap43","Nefl","Nefm")
wang_schwann = c("Mpz","Plp1")
wang_fibroblasts = c("Dcn","Apod")
wang_endothelial = c("Cldn5","Ly6c1")
wang_smooth_muscle = c("Acta2","Tagln")
wang_macrophages = c("Mrc1","Cd68")
wang_capillary = c("S100a8")
wang_immune = c("Cd74")
wang_rbc = c("Hba-a1")
wang_rubbish = c("Malat1")
wang = c(wang_schwann,wang_fibroblasts,wang_endothelial,wang_smooth_muscle,wang_macrophages,wang_capillary,wang_immune,wang_rbc,wang_neurons,wang_rubbish)


# SingleR annotation ---------------------------- #
library(celldex)
library(SingleR)

merged_seurat = JoinLayers(merged_seurat)
combined = merged_seurat # in order to have a backup
ref = celldex::MouseRNAseqData()
sc_counts = GetAssayData(combined,layers="data") 

# SingleR annotation
pred = SingleR(test= sc_counts, ref= ref, label= ref$label.fine) 
combined$singleR.fine = pred$labels[match(rownames(combined@meta.data),rownames(pred))]
pred = SingleR(test= sc_counts, ref= ref, label= ref$label.main) 
combined$singleR.main = pred$labels[match(rownames(combined@meta.data),rownames(pred))]

# First annotation attempt with major cell type markers
# Rename idents by celltype
new.cluster.ids = c("Neurons","Neurons","Schwann","SGC","Schwann","Endothelial","Immune","Neurons","cluster8","cluster9","cluster10")
annotated = combined
names(new.cluster.ids) <- levels(annotated)
annotated <- RenameIdents(annotated, new.cluster.ids)
annotated$celltype = Idents(annotated)
save(annotated,combined,markers,file="start_pseudobulk.RData")


neurons <- subset(annotated, celltype == "Neurons")
non_neurons <- subset(annotated, celltype != "Neurons")



