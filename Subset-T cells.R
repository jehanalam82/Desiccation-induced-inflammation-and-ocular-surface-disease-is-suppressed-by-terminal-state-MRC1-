library(Seurat)
library(dplyr)
wd <- "H:/Manuscripts/DS-Model_ScRNA_Seq_Analysis/Subset-Tcells" # MODIFY
setwd(wd)
filename <- "H:/Manuscripts/DS-Model_ScRNA_Seq_Analysis/DS_merged.rds" # MODIFY
x <-readRDS(filename)
DefaultAssay(x) <- "RNA"
x <- FindClusters(object = x, resolution = 0.4)# Find clusters and visualize
DimPlot(object = x, group.by = "RNA_snn_res.0.4", label = TRUE)   ### MODIFY based on resolution you used
######################Sub-Setting-Tcells###########################################################################
saveRDS(subset(x, ident = c("6")), "H:/Manuscripts/DS-Model_ScRNA_Seq_Analysis/Subset-Tcells/Cluster_6.rds")
Tcells <- "H:/Manuscripts/DS-Model_ScRNA_Seq_Analysis/Subset-Tcells/Cluster_6.rds" # MODIFY
x <-readRDS(Tcells)
DefaultAssay(x) <- "RNA"
x<- NormalizeData(x, normalization.method = "LogNormalize", scale.factor = 10000)
x<- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
x<- ScaleData(x, verbose = FALSE)
x<- RunPCA(x, npcs = 50, verbose = FALSE)
x<- RunUMAP(x, reduction = "pca", dims = 1:30)
x<- FindNeighbors(x, reduction = "pca", dims = 1:30)
x<- FindClusters(x, resolution = 0.4)
p1 <- DimPlot(x, reduction = "umap" , group.by = "sample")
p2 <- DimPlot(x, reduction = "umap", label = TRUE, repel = TRUE)
p1
p2
current.cluster.ids <- c(0, 1, 2, 3, 4, 5)
new.cluster.ids <- c("CD4-EFF", "Treg", "CD8-Tn", "CD4-Tn", "NKT", "CD8-EFF")
x@active.ident <- plyr::mapvalues(x = x@active.ident, from = current.cluster.ids, to = new.cluster.ids)
DimPlot(x, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend() 


#######Cells-count################################################################################################################
x[["RNA_snn_res.0.4"]] <- Idents(object = x)
table(x@meta.data$RNA_snn_res.0.4, x@meta.data$days)
####################Cluster-Specific-markers##########################
Clusters.markers <- FindAllMarkers(x, only.pos = FALSE, min.pct = 0.25, logfc.threshold = 0.25)
#Selecting Top 1000 Genes 
top1000 <- Clusters.markers %>%
  filter(p_val_adj < 0.05) %>% #significant filter
  slice_max(abs(avg_log2FC), n = 1000) %>% #Using absolute for high or low genes
  arrange(avg_log2FC) #Order by value
top1000genes <- rownames(top1000)
write.csv(top1000, file = "H:/Manuscripts/DS-Model_ScRNA_Seq_Analysis/Subset-Tcells/TCells_refine" )

#############################################################Dot-Plot#############################################################
DotPlot(x, features = c("Ptprc","Cd3d","Cd3g", "Cd3e","Cd4","Cd8a",
                        "Il17a", "Ifng","Prf1",
                        "Gzma", "Gzmb","Gzmk","Nkg7",
                        "Ccl5","Xcl1","Cd7","Ctsw",
                        "Cd28",
                        "Ccr7","Lef1","Sell", 
                        "Ctla4","Foxp3","Ikzf2"
                        
                        ),
        cols = c("blue", "red"), dot.scale = 8) + theme(axis.text.x = element_text(angle = 90, hjust=1))

