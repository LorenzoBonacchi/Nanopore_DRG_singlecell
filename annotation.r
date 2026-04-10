library(Seurat)
library(harmony)
library(DoubletFinder)
library(tibble)
library(celldex)
library(SingleR)

load("data_markers.RData")
combined = merged_seurat # in order to have a backup
ref = celldex::HumanPrimaryCellAtlasData()
sc_counts = GetAssayData(combined,layers="data") 

# SingleR annotation
pred = SingleR(test= sc_counts, ref= ref, label= ref$label.fine) 
combined$singleR.fine = pred$labels[match(rownames(combined@meta.data),rownames(pred))]
pred = SingleR(test= sc_counts, ref= ref, label= ref$label.main) 
combined$singleR.main = pred$labels[match(rownames(combined@meta.data),rownames(pred))]

# Set markers used by Liu et al. 
kt = c("KRT1","KRT5", "KRT10", "KRT14") #keratinocytes 
fb = c("DCN", "COL1A1", "COL1A2") #fibroblasts
dmac = c("PTPRC", "CD68", "CD1C") #dendritic_macrophages
t_nk = c("PTPRC", "CD3D", "GNLY", "NKG7") 
endo = c("CD93", "ACKR1", "AQP1") #endothelial 
lym = c("CCL21", "LYVE1", "TFF3") #lymphatic
mel = c("TYRP1", "PMEL", "DCT") #melanocytes
sweat = c("DCD", "KRT19", "AQP5") #sweatglandcells
smooth = c("TAGLN", "ACTA2", "TPM2") #smooth_muscle


all_markers = c("KRT1","KRT5", "KRT10", "KRT14","DCN", "COL1A1", "COL1A2","PTPRC", "CD68", 
    "CD1C","CD3D", "GNLY","NKG7","CD93", "ACKR1", "AQP1","CCL21", "LYVE1", "TFF3",
    "TYRP1", "PMEL", "DCT","DCD", "KRT19", "AQP5", "TAGLN", "ACTA2", "TPM2") 
# I used all_markers in order to check all of them in a dotplot
png("dotplot_sc_bp.png",width=1800,height=1800)
DotPlot(annotated,features = all_markers)
dev.off()

# Rename idents by celltype
new.cluster.ids = c("keratinocytes","keratinocytes","keratinocytes","keratinocytes",
"smoothmuscle","fibroblasts","endothelial","keratinocytes","keratinocytes","dmac",
"melanocytes","keratinocytes","lymphatic","T_cells","keratinocytes","keratinocytes","keratinocytes",
"sweet_glands_cells","endothelial","fibroblasts","keratinocytes")

annotated = combined
names(new.cluster.ids) <- levels(annotated)
annotated <- RenameIdents(annotated, new.cluster.ids)
annotated$celltype = Idents(annotated)

save(annotated,combined,markers,file="start_pseudobulk.RData")