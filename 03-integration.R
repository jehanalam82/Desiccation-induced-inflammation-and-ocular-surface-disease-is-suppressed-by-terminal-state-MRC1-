#!/usr/bin/env Rscript
# vim: set noexpandtab tabstop=2:

library(dplyr)
library(Seurat)

src <- "/storage/singlecell/kxiong/DSmodel/clustering/"

# uncomment below to perform CCA integration
#prefix <- paste0("/storage/singlecell/kxiong/DSmodel/integration/DS_integrated")
#mergeOnly = FALSE
#all.features = FALSE

# uncomment below to merge datasets
prefix <- paste0("/storage/singlecell/kxiong/DSmodel/integration/DS_merged")
mergeOnly = TRUE
all.features = FALSE

if (mergeOnly == TRUE) {
  clusters <- "RNA_snn_res.0.4"
} else {
  clusters <- "integrated_snn_res.0.5"
}

###############################################################################

# Helper function for Seurat data integration
integrate <- function(data.list,
                      normalization.method = "LogNormalize",
                      nfeatures = 2000,
                      dims = 1:30, k.weight = 100,
                      all.features = FALSE) {

  # Log normalization method
  if (normalization.method == "LogNormalize") {
    data.list <- lapply(data.list, function(x) {
      x <- NormalizeData(x)
      x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = nfeatures)
      return(x)
    })
    features <- SelectIntegrationFeatures(object.list = data.list, nfeatures = nfeatures)
  }

  # SCTransform method
  if (normalization.method == "SCT") {
    data.list <- lapply(X = data.list, FUN = SCTransform)
    features <- SelectIntegrationFeatures(object.list = data.list, nfeatures = nfeatures)
    data.list <- PrepSCTIntegration(object.list = data.list, anchor.features = features)
  }

  # Perform integration
  data.anchors <- FindIntegrationAnchors(object.list = data.list, normalization.method = normalization.method, anchor.features = features, dims = dims)
  if (all.features == TRUE) {
    all.genes <- Reduce(intersect, lapply(data.list, rownames))
    x <- IntegrateData(anchorset = data.anchors, normalization.method = normalization.method, dims = dims, k.weight = k.weight, features.to.integrate = all.genes)
  } else {
    x <- IntegrateData(anchorset = data.anchors, normalization.method = normalization.method, dims = dims, k.weight = k.weight)
  }

  DefaultAssay(x) <- "integrated"

  # Run the standard workflow for visualization and clustering
  x <- ScaleData(x, verbose = FALSE)
  x <- RunPCA(x, npcs = 30, verbose = FALSE)
  x <- RunTSNE(x, dims = 1:30)
  x <- RunUMAP(x, reduction = "pca", dims = 1:30)
  x <- FindNeighbors(x, reduction = "pca", dims = 1:30)
  x <- FindClusters(x, resolution = 0.9)
  x <- FindClusters(x, resolution = 0.7)
  x <- FindClusters(x, resolution = 0.5)
  x <- FindClusters(x, resolution = 0.3)
  x <- FindClusters(x, resolution = 0.1)

  return(x)
}

# Helper function to merge Seurat objects
mergeSeurat <- function(data.list) {
  x <- merge(data.list[[1]], y = data.list[-1])
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
  x <- ScaleData(x, verbose = FALSE)
  x <- RunPCA(x, npcs = 30, verbose = FALSE)
  x <- RunTSNE(x, dims = 1:30)
  x <- RunUMAP(x, reduction = "pca", dims = 1:30)
  x <- FindNeighbors(x, reduction = "pca", dims = 1:30)
  x <- FindClusters(x, resolution = 0.9)
  x <- FindClusters(x, resolution = 0.7)
  x <- FindClusters(x, resolution = 0.5)
  x <- FindClusters(x, resolution = 0.3)
  x <- FindClusters(x, resolution = 0.1)
  return(x)
}


# Load individual datasets 
sample.list <- list.dirs(path = src, full.names = FALSE, recursive = FALSE)
data.list <- lapply(sample.list, function(sample) {
  file <- paste0(src, "/", sample, "/", sample, "_cluster.rds")
  x <- readRDS(file)
  x <- DietSeurat(x)
  x$sample <- sample
	x$cell.id <- paste0(x$sample,"_",names(x$orig.ident))
  return(x)
})

if (mergeOnly == TRUE) {
	# Merge datasets to check for batch effects
	x <- mergeSeurat(data.list)
} else {
	# Integrate datasets
	x <- integrate(data.list, all.features = all.features)
	#x <- integrate(data.list, normalization.method = "SCT", nfeatures = 3000, all.features = all.features)
}
saveRDS(x, paste0(prefix, ".rds"))


#################################################################################

# Modify metadata
x <- readRDS(paste0(prefix, ".rds"))
x <- FindClusters(x, resolution = 0.4)
x$days <- plyr::mapvalues(x$sample, from = c("NS","DS1","DS3","DS5","DS10"), to = c(0, 1, 3, 5, 10))
x$days <- as.numeric(x$days)

# Cell type annotation by clusters
x <- SetIdent(x, value=clusters)
Idents(x) <- factor(x = Idents(x), levels = as.character(0:20))
x$seurat_clusters <- Idents(x)
new.cluster.ids <- c("MF","cDC2/MoDC","MF/cDC2/MoDC","Monocyte","gd-T cell", "ILC2",
                     "T cell","NK cell","Mast cell","cDC1","MF",
                     "MF","MF","Neutrophil", "Pre-B/Pre-T cell","T cell/NK cell",
                     "B cell","mDC?", "B cell/DC","Stromal","pDC")
names(new.cluster.ids) <- levels(x)
x <- RenameIdents(x, new.cluster.ids)
x$celltype <- Idents(x)

pdf(paste0(prefix,"_celltype.pdf"), width=8, height=8)

# Plot UMAP
p2 <- DimPlot(x, reduction = "umap", group.by = "celltype", label=TRUE) + NoLegend()
print(p2)

# Plot TSNE
p2 <- DimPlot(x, reduction = "tsne", group.by = "celltype", label=TRUE) + NoLegend()
print(p2)

dev.off()

saveRDS(x, paste0(prefix, ".rds"))

###############################################################################


# Visualization
x <- readRDS(paste0(prefix, ".rds"))

pdf(paste0(prefix,".pdf"), width=8, height=8)

# Plot UMAP
p1 <- DimPlot(x, reduction = "umap", group.by = "sample", shuffle = TRUE)
p2 <- DimPlot(x, reduction = "umap", group.by = clusters, label=TRUE, label.size = 8) + NoLegend()
print(p1)
print(p2)

# Plot TSNE
p1 <- DimPlot(x, reduction = "tsne", group.by = "sample", shuffle = TRUE)
p2 <- DimPlot(x, reduction = "tsne", group.by = clusters, label=TRUE, label.size = 8) + NoLegend()
print(p1)
print(p2)

# Plot markers
DefaultAssay(x) <- "RNA"

features <- c("Cd79a","Cd28","S100a9","Gata3","Gzma","Mcpt4","Cd209a","Trdc","Fscn1","Xcr1","Cd4","Foxp3","Cd8a")

p1 <- DotPlot(x, features = features, group.by = clusters) + RotatedAxis()
print(p1)
dev.off()

png(paste0(prefix, "_features.png"), width=26, height=20, units="in", res=300)
print(FeaturePlot(x, features = features))
dev.off()

# Plot split by sample
png(paste0(prefix,"_split.png"), width=12, height=4, units = "in", res=300)
p1 <- DimPlot(x, reduction = "umap", group.by = clusters, split.by = "sample", label = TRUE, repel = TRUE, label.size = 5)
print(p1)
dev.off()

