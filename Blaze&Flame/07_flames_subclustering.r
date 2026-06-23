
neurons <- RunPCA(neurons)
neurons <- FindNeighbors(neurons, dims = 1:30)
neurons <- FindClusters(neurons, resolution = 1)
neurons <- RunUMAP(neurons, dims = 1:30)

non_neurons <- RunPCA(non_neurons)
non_neurons <- FindNeighbors(non_neurons, dims = 1:30)
non_neurons <- FindClusters(non_neurons, resolution = 1)
non_neurons <- RunUMAP(non_neurons, dims = 1:30)

library(celldex)
library(SingleR)

non_neurons = JoinLayers(non_neurons)
combined_non_neurons = non_neurons # in order to have a backup
ref = celldex::MouseRNAseqData()
sc_counts = GetAssayData(combined_non_neurons,layers="data") 

# SingleR annotation
pred = SingleR(test= sc_counts, ref= ref, label= ref$label.fine) 
combined_non_neurons$singleR.fine = pred$labels[match(rownames(combined_non_neurons@meta.data),rownames(pred))]
pred = SingleR(test= sc_counts, ref= ref, label= ref$label.main) 
combined_non_neurons$singleR.main = pred$labels[match(rownames(combined_non_neurons@meta.data),rownames(pred))]




qNSCs = c("Hes5","Sox2","Notch1","Gfap","Id1","Hey1","Sox9","Bmp2","Hes1")
NPCs = c("Mki67","Top2a","Pcna","Ascl1","Neurod1","Dcx","Stmn2","Egfr","Sox4")
transition = c("Ascl1","Egfr","Sox4","Hes1")

neurons_markers = c("Snap25","Rbfox3","Pvalb","Gfra3","Tac1","Calca","Gal","Cldn9","Zcchc12","Sstr2","Dcn","Trpm8",
  "Rxfb1","Nppb","Th","Fam19a4","Mrgprb4","Mrgpra3","Mrgprd","Lpar3","Gm7271","S100b","Nefh","Wnt7a","Trappc3l","Ntrk3","Gfra1","Prokr2","smr2","Baiap2l1","Atf3")

## Notebook markers
neurons_gen = c("Snap25", "Rbfox3", "Pvalb", "Nefh", "Tubb3")
sgc_gen = c("Fabp7", "Ednrb", "Apoe", "Slc1a3", "Gja1")
mysc = c("Mpz", "Mbp", "Pllp")
nomysc = c("Scn7a", "Ngfr", "Ncam1", "L1cam") 
fibro = c("Pdgfra", "Tbx18", "Dcn", "Apod")
endo = c("Pecam1", "Flt1", "Cldn5")
immune = c("Ptprc", "Ccr2", "Cd74", "Aif1", "Cd68")
pericytes = c("Notch3", "Kcnj8", "Pdgfrb")
proprioceptors = c("Pvalb", "Ntrk3", "Etv1", "Runx3")
lowthr_alphabeta = c("Ntrk3", "Ntrk2")
lowthr_alphaomega = c("Trappc3l", "Ntrk2", "Gfra2")
nociceptors = c("Tac1", "Calca", "Adcyap1")
termo = c("Trpv1", "Trpm8")
pep_alphaomega = c("Ntrk1", "Nefh")
np1 = c("Mrgprd", "Lpar3", "Gfra2")
np2 = c("Mrgpra3", "Mrgprx1")
np3 = c("Sst", "Nppb", "Il31ra", "Osmr", "Jak1")
cltmr = c("Th", "Fam19a4")
scprecursor = c("Sox10", "Mki67", "Top2a")
immature_sc = c("Ngfr", "Ncam1", "L1cam")
pro_myelin_sc = c("Pou3f1", "Cdkn1c")


all_neu = c(neurons_gen, sgc_gen, mysc, nomysc, fibro, endo, immune, pericytes, proprioceptors, lowthr_alphabeta, lowthr_alphaomega, nociceptors, termo, pep_alphaomega, np1, np2, np3, cltmr)
all_sc = c(sgc_gen, mysc, nomysc, scprecursor, immature_sc, pro_myelin_sc)




















# 1. Isola i cluster identificati come neuroni nel primo round
neurons_subset <- subset(annotated, idents = c("Neurons_0", "Neurons_1", "Neurons_6","Neurons_7", "Neurons_11"))

# 2. Ripristina l'assay e ripulisci i layer (se usi Seurat v5)
DefaultAssay(neurons_subset) <- "RNA"

# 3. Ricalcola le feature variabili SOLO sui neuroni 
# (ora i geni come Scn10a, Tac1, Th guideranno la PCA!)
neurons_subset <- NormalizeData(neurons_subset)
neurons_subset <- FindVariableFeatures(neurons_subset, nfeatures = 2000)

# Rimuovi di nuovo i mitocondriali dalle HVG per sicurezza
var_genes <- VariableFeatures(neurons_subset)
var_genes <- var_genes[!grepl("^mt-|^Mt-", var_genes)]
VariableFeatures(neurons_subset) <- var_genes

# 4. Riallinea con Harmony per correggere i batch tra i neuroni
neurons_subset <- ScaleData(neurons_subset, features = VariableFeatures(neurons_subset))
neurons_subset <- RunPCA(neurons_subset, verbose = FALSE)
neurons_subset <- RunHarmony(neurons_subset, group.by.vars = "condition")

# 5. Calcola UMAP e Neighbors sullo spazio Harmony dei soli neuroni
neurons_subset <- RunUMAP(neurons_subset, reduction = "harmony", dims = 1:15) # Spesso servono meno PC per i soli neuroni
neurons_subset <- FindNeighbors(neurons_subset, reduction = "harmony", dims = 1:15)

# 6. Clustering a risoluzione medio-alta per catturare i sottotipi
neurons_subset <- FindClusters(neurons_subset, resolution = 0.8)

pep_genes <- list(c("Tac1", "Calca", "Ntrk1", "Sst"))
np_genes <- list(c("Mrgprd", "Mrgpra3", "Chrna3"))
nf_genes <- list(c("Nefh", "Pvalb", "Ntrk3", "Runx3"))

neurons_subset <- AddModuleScore(neurons_subset, features = pep_genes, name = "PEP_Score")
neurons_subset <- AddModuleScore(neurons_subset, features = np_genes, name = "NP_Score")
neurons_subset <- AddModuleScore(neurons_subset, features = nf_genes, name = "NF_Score")

# Visualizzali graficamente sulla UMAP dei soli neuroni
FeaturePlot(neurons_subset, features = c("PEP_Score1", "NP_Score1", "NF_Score1"), ncol = 3)