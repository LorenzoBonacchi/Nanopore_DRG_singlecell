

library(dplyr)
library(harmony)

merged_seurat <- NormalizeData(merged_seurat)
merged_seurat <- FindVariableFeatures(merged_seurat, nfeatures = 3000)
merged_seurat <- ScaleData(merged_seurat)

merged_seurat <- RunPCA(merged_seurat, npcs = 30)

merged_seurat <- RunHarmony(
  merged_seurat,
  group.by = "orig.ident"
)

#resolutions <- c(0.1, 0.2, 0.3, 0.4, 0.5, 1.0) 
resolutions <- c(0.3) 
merged_seurat <- RunUMAP(merged_seurat, reduction = "harmony", dims = 1:30)
merged_seurat <- FindNeighbors(merged_seurat, reduction = "harmony", dims = 1:30)
merged_seurat <- FindClusters(merged_seurat, resolution = resolutions)

neurons = c("Snap25", "Rbfox3", "Pvalb")
sgc = c("Fabp7", "Ednrb")
myelinating_schwann = c("Mpz", "Mbp")
nonmyelinating_schwann = c("Scn7a")
immune = c("Ptprc", "Ccr2")
endothelial = c("Pecam1", "Flt1")
pericytes = c("Notch3", "Kcnj8")
fibroblasts = c("Pdgfra", "Tbx18")
all = c(neurons, sgc, myelinating_schwann, nonmyelinating_schwann, immune, endothelial, pericytes, fibroblasts)

library(celldex)
library(SingleR)


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
new.cluster.ids = c("Neurons","Neurons","Neurons","Neurons","Satellite","Mschwann","Endothelial","Immune","Neurons","Neurons")
annotated = combined
names(new.cluster.ids) <- levels(annotated)
annotated <- RenameIdents(annotated, new.cluster.ids)
annotated$celltype = Idents(annotated)
save(annotated,combined,markers,file="start_pseudobulk.RData")
