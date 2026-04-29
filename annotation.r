library(Seurat)
library(harmony)
library(DoubletFinder)
library(tibble)
library(celldex)
library(SingleR)

load("data_markers.RData")
combined = merged_seurat_filtered # in order to have a backup
ref = celldex::MouseRNAseqData()
sc_counts = GetAssayData(combined,layers="data") 

# SingleR annotation
pred = SingleR(test= sc_counts, ref= ref, label= ref$label.fine) 
combined$singleR.fine = pred$labels[match(rownames(combined@meta.data),rownames(pred))]
pred = SingleR(test= sc_counts, ref= ref, label= ref$label.main) 
combined$singleR.main = pred$labels[match(rownames(combined@meta.data),rownames(pred))]

# Set markers used by Liu et al. 

neurons = c("Snap25","Syt1")
glial = c("Sox10","Plp1")
fibroblasts = c("Col6a1", "Col1a1")
endothelial = c("Pecam1","Flt1")
immune = c("Aif1", "Ptprc")
mural = c("Rgs5","Acta2")
erythrocytes = c("Hba-a1")
    

all_markers = c("Snap25","Syt1","Sox10","Plp1","Col6a1", "Col1a1","Pecam1","Flt1","Aif1", "Ptprc","Rgs5","Acta2","Hba-a1","Mpz")
png("dotplot_drg.png",width=1800,height=1800)
DotPlot(combined,features = all_markers)
dev.off()


snap <- FetchData(combined, "Snap25")[,1]
mpz  <- FetchData(combined, "Mpz")[,1]
pec  <- FetchData(combined, "Pecam1")[,1]
ptprc <- FetchData(combined, "Ptprc")[,1]
col1a1 <- FetchData(combined, "Col1a1")[,1]
combined$celltype_simple <- "unknown"
combined$celltype_simple[snap > 1] <- "Neuron"
combined$celltype_simple[mpz > 1] <- "Schwann"
combined$celltype_simple[pec > 1] <- "Endothelial"
combined$celltype_simple[ptprc > 1] <- "Immune"
combined$celltype_simple[col1a1 > 0.3] <- "Fibroblast"

## Neurons
neurons <- subset(combined, subset = Snap25 > 1)
neurons <- DietSeurat(
  neurons,
  counts = TRUE,
  data = TRUE,
  scale.data = TRUE
)
neurons <- NormalizeData(neurons)
neurons <- FindVariableFeatures(neurons)
neurons <- ScaleData(neurons)
neurons <- RunPCA(neurons)
neurons <- FindNeighbors(neurons, dims = 1:30)  
neurons <- FindClusters(neurons, resolution = 0.8)
neurons <- RunUMAP(neurons, dims = 1:30)

## Non neurons
non_neurons <- subset(combined, Snap25 <= 1)
non_neurons@graphs <- list()
non_neurons@neighbors <- list()
non_neurons@reductions <- list()
non_neurons <- NormalizeData(non_neurons)
non_neurons <- FindVariableFeatures(non_neurons)
non_neurons <- ScaleData(non_neurons)
non_neurons <- RunPCA(non_neurons)
non_neurons <- FindNeighbors(non_neurons)
non_neurons <- FindClusters(non_neurons, resolution = 0.8)
non_neurons <- RunUMAP(non_neurons)

markers <- c("Snap25","Syt1","Sox10","Plp1","Col6a1","Col1a1",
             "Pecam1","Flt1","Aif1","Ptprc","Rgs5","Acta2","Hba-a1","Mpz")
# =================================================== #
# Nuova annotazione con markers
# =================================================== #                                            
neurons = c("Snap25","Syt1","Avil","Gap43","Nefl","Nefm")
glial = c("Sox10","Plp1","Fabp7","Tyrp1","Hmgcs2","Slc1a3")
fibroblasts = c("Col6a1", "Col1a1","Pdgfra","Dpp4","Lum","Egfr")
endothelial = c("Pecam1","Flt1","Cldn5","Emcn","Prom1")
macro = c("Mrc1","Ptprc","Cd68","Cd86","Csf1r")
mural = c("Rgs5","Acta2","Tagln","Des")
erythrocytes = c("Hba-a1")
schwann = c("Mpz","Mag","Scn7a","Pou3f1","Ncam1")
avg_exp_neurons <- AverageExpression(combined, features = neurons, group.by = "seurat_clusters")
avg_exp_glial <- AverageExpression(combined, features = glial, group.by = "seurat_clusters")
avg_exp_fibro <- AverageExpression(combined, features = fibroblasts, group.by = "seurat_clusters")
avg_exp_endo <- AverageExpression(combined, features = endothelial, group.by = "seurat_clusters")
avg_exp_macro <- AverageExpression(combined, features = macro, group.by = "seurat_clusters")
avg_exp_mural <- AverageExpression(combined, features = mural, group.by = "seurat_clusters")
avg_exp_erythro <- AverageExpression(combined, features = erythrocytes, group.by = "seurat_clusters")
avg_exp_schwann <- AverageExpression(combined, features = schwann, group.by = "seurat_clusters")

avg_exp_neurons$RNA
avg_exp_glial$RNA
avg_exp_fibro$RNA
avg_exp_endo$RNA
avg_exp_macro$RNA
avg_exp_mural$RNA
avg_exp_erythro$RNA
avg_exp_schwann$RNA

neurons = c("Snap25","Syt1","Avil","Gap43","Nefl","Nefm")
glial = c("Sox10","Plp1","Fabp7","Tyrp1","Hmgcs2","Slc1a3")
fibroblasts = c("Col6a1", "Col1a1","Pdgfra","Dpp4","Lum","Egfr")
endothelial = c("Pecam1","Flt1","Cldn5","Emcn","Prom1")
macro = c("Mrc1","Ptprc","Cd68","Cd86","Csf1r")
mural = c("Rgs5","Acta2","Tagln","Des")
erythrocytes = c("Hba-a1")
schwann = c("Mpz","Mag","Scn7a","Pou3f1","Ncam1")
avg_exp_neurons <- AverageExpression(non_neurons, features = neurons, group.by = "seurat_clusters")
avg_exp_glial <- AverageExpression(non_neurons, features = glial, group.by = "seurat_clusters")
avg_exp_fibro <- AverageExpression(non_neurons, features = fibroblasts, group.by = "seurat_clusters")
avg_exp_endo <- AverageExpression(non_neurons, features = endothelial, group.by = "seurat_clusters")
avg_exp_macro <- AverageExpression(non_neurons, features = macro, group.by = "seurat_clusters")
avg_exp_mural <- AverageExpression(non_neurons, features = mural, group.by = "seurat_clusters")
avg_exp_erythro <- AverageExpression(non_neurons, features = erythrocytes, group.by = "seurat_clusters")
avg_exp_schwann <- AverageExpression(non_neurons, features = schwann, group.by = "seurat_clusters")

avg_exp_neurons$RNA
avg_exp_glial$RNA
avg_exp_fibro$RNA
avg_exp_endo$RNA
avg_exp_macro$RNA
avg_exp_mural$RNA
avg_exp_erythro$RNA
avg_exp_schwann$RNA
all = c(neurons,glial,fibroblasts,endothelial,macro,mural,erythrocytes,schwann)
DotPlot(non_neurons,features = all) 




# =================================================== #
# Neuron subtypes annotation
# =================================================== #

n1 = c("Tac1","Sst","Stmn3","Scg2","Ntrk1")
n2 = c("Kcnip4","Tmem35a","Alg2","Rab4a","Serp2")
n3 = c("Strn3","Notch1","Lgi4","Cdh19","Itgb4")
n4 = c("Bcl2","Cdh15","Tap1","Col12a1","Pcdh9")
n5 = c("Acot1","Maff","Phlda3","Cnp","Gsta4")
n6 = c("Ralgps1","Prdm8","Ccser1","Ptpn5","Col6a5")


# =================================================== #
# Gabri celltype annotation
# =================================================== #

neurons = c("Snap25","Syt1","Avil","Gap43","Nefl","Nefm")
glial = c("Sox10","Plp1","Fabp7","Tyrp1","Hmgcs2","Slc1a3")
fibroblasts = c("Col6a1", "Col1a1","Pdgfra","Dpp4","Lum","Egfr")
endothelial = c("Pecam1","Flt1","Cldn5","Emcn","Prom1")
macro = c("Mrc1","Ptprc","Cd68","Cd86","Csf1r")
mural = c("Rgs5","Acta2","Tagln","Des")
erythrocytes = c("Hba-a1")
schwann = c("Mpz","Mag","Scn7a","Pou3f1","Ncam1","S100b","Ngfr")

neurons = c("Snap25","Syt1","Avil","Gap43","Nefl","Nefm","Tubb3","Rbfox3","Elavl4","Prph")
glial = c("Sox10","Plp1","Fabp7","Slc1a3","Glul","Kcnj10","S100b","Gja1")
fibroblasts = c("Col6a1","Col1a1","Col1a2","Col3a1","Pdgfra","Dpp4","Lum","Pi16","Dpt","Sfrp4","Prrx1","Aebp1")
endothelial = c("Pecam1","Flt1","Cldn5","Emcn","Prom1","Cdh5","Kdr","Tek")
macro = c("Mrc1","Ptprc","Cd68","Cd86","Csf1r","Adgre1","Tyrobp","Lyz2","Fcgr1","C1qa")
mural = c("Rgs5","Acta2","Tagln","Des","Pdgfrb","Cspg4","Mcam","Notch3")
erythrocytes = c("Hba-a1","Hbb-bs","Slc4a1","Gypa","Alas2")
schwann = c("Mpz","Mag","Scn7a","Pou3f1","Ncam1","S100b","Ncmap","Cdh19","Ngfr","Sox10","Plp1")

##########################################################

non_peptidergic_nociceptors_NP1 = c("Nppb")
non_peptidergic_nociceptors_NP2 = c("Mrgprd","Cd55")
non_peptidergic_nociceptors_NP3 = c("Mrgpra3","Cd55")
peptidergic_nociceptors_PEP1_PEP2 = c("Tac1")
peptidergic_nociceptors_PEP3_PEP4 = c("Tac1","Sstr2")
peptidergic_nociceptors_PEP5_cold_associated = c("Trpm8")
peptidergic_nociceptors_PEP6_heat_associated = c("Trpv1")
myelinated_Abeta_low_threshold_mechanoreceptors_NF1_NF2 = c("Nefh","Scn1b")
C_fiber_low_threshold_mechanoreceptors_cLTMR = c("Fam19a4","Th")
injury_associated_neurons = c("Atf3","Sprr1a")
########################################################

# ==================================================== #
# ==================================================== #
# Wang annotation
neurons = c("Snap25","Syt1","Avil","Gap43","Nefl","Nefm")
schwann = c("Mpz","Plp1")
fibroblasts = c("Dcn","Apod")
endothelial = c("Cldn5","Ly6c1")
smooth_muscle = c("Acta2","Tagln")
macrophages = c("Mrc1","Cd68")
capillary = c("S100a8")
immune = c("Cd74")
rbc = c("Hba-a1")
all = c(schwann,fibroblasts,endothelial,smooth_muscle,macrophages,capillary,immune,rbc,neurons)

# ==================================================== #
# Elife Schwann annotation


# ==================================================== #
# ==================================================== #
# ==================================================== #
# ==================================================== #
# ==================================================== #
# ==================================================== #
# ==================================================== #
# ==================================================== #
### TO DO DOPO ANNOTAZIONE finale:

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