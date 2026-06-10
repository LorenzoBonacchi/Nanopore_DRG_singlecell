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
markers_top = as.data.frame(markers_top)


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
wang_neurons = c("Snap25","Syt1","Avil","Gap43","Nefl","Nefm","Nefh","Calca")
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


#score

axon_genes <- c("Nefl", "Tubb2a", "Tubb2b", "Stmn3", "Gphn", "Cadm2")
myelin_genes <- c("Mbp", "Mpz", "Pmp2", "Plp1", "Mal", "Mag")
glia_support_genes <- c("Apoe", "Fabp7", "Ndrg1", "S100b", "Cnp", "Fxyd1")
neuron_genes <- c("Snap25","Syp","Rbfox3","Tubb3","Map2","Nefl","Syn1","Syt1")

merged_seurat <- AddModuleScore(merged_seurat,features = list(axon_genes),name = "axon_score")
merged_seurat <- AddModuleScore(merged_seurat, features = list(myelin_genes), name = "myelin_score")
merged_seurat <- AddModuleScore(merged_seurat, features = list(glia_support_genes), name = "glia_score")
merged_seurat <- AddModuleScore(merged_seurat, features = list(neuron_genes), name = "neuron_score")

FeaturePlot(merged_seurat, features=c("axon_score1","myelin_score1","glia_score1","neuron_score1"))

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
new.cluster.ids = c("Neurons_0","Neurons_1","Schwann_2","SGC_3","Schwann_4","Endothelial","Neurons_6","Neurons_7","SGC_8","Immune_9","Immune_10","Neurons_11")
annotated = combined
names(new.cluster.ids) <- levels(annotated)
annotated <- RenameIdents(annotated, new.cluster.ids)
annotated$celltype = Idents(annotated)
save(annotated,combined,markers,file="start_pseudobulk.RData")


#neurons <- subset(annotated, celltype == "Neurons")
#non_neurons <- subset(annotated, celltype != "Neurons")

# fine annotation
annotated <- FindClusters(annotated, resolution = 1)
annotated = JoinLayers(annotated)
markers_fine <- FindAllMarkers(
  annotated,
  only.pos = TRUE,
  test.use = "wilcox",
  logfc.threshold = 0.25,
  min.pct = 0.25
)
markers_top <- markers_fine %>%
  filter(p_val_adj < 0.05,
         pct.1 > 0.25,
         avg_log2FC > 0.25)
markers_top <- markers_top %>%
  group_by(cluster) %>%
  slice_max(avg_log2FC, n = 20, with_ties = FALSE) %>%
  ungroup()
markers_top = as.data.frame(markers_top)


neuronal = c("Tac1", "Calca", "Ntrk1", "Ntrk2", "Ntrk3", "Ret", "Piezo2", "Mrgprd", "Mrgpra3", "Pvalb", "Runx3", "Th", "Sst", "Nefh", "Scn10a", "Scn11a")
neuronal_specific = c("Ntrk3", "Nefh","Gphn","Pvalb", "Fxyd2","Pcp4", "Calca", "Prph", "Pcsk1n")

VlnPlot(
  neurons,
  features = c(
    "mitoRatio",
    "nFeature_RNA",
    "nCount_RNA"
  ),
  group.by = "seurat_clusters"
)


