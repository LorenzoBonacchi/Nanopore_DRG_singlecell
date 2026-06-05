
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