#!/usr/bin/env Rscript
# vim: set noexpandtab tabstop=2:

library("optparse")
 
option_list = list(
  make_option(c("-f", "--file"), type="character", default=NULL, 
              help="Seurat Object .rds file name", metavar="character"),
  make_option(c("-n", "--name"), type="character", default=NULL, 
              help="sample name", metavar="character"),
  make_option(c("-o", "--out"), type="character", default=NULL, 
              help="output directory (ie. /dir/) [default= %default]", metavar="character")
); 
 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);


# Load libraries
library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(tools)

f=opt$file
prefix=file.path(opt$out, opt$name)
dir.create(opt$out, recursive = TRUE)

# Load data
x <- readRDS(f)
x <- x[['RNA']]@counts
x <- CreateSeuratObject(counts = x, min.cells = 3, min.features = 200)
x[["percent.mt"]] <- PercentageFeatureSet(x, pattern = "^MT-|^mt-")
x

# Normalize the data
x <- NormalizeData(x)

# HVGs
x <- FindVariableFeatures(object = x, selection.method = 'vst', nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(x), 10)
top10

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(x)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
png(paste0(prefix, "_hvg.png"), width=8, height=8, units="in", res=300)
print(plot2)
dev.off()

# Scale the data
all.genes <- rownames(x)
x <- ScaleData(x, features = all.genes)

# Dimensional reduction
x <- RunPCA(x, features = VariableFeatures(object=x))
png(paste0(prefix, "_elbow.png"), width=8, height=8, units="in", res=300)
print(ElbowPlot(x))
dev.off()

x <- RunTSNE(x, dims=1:30)
x <- RunUMAP(x, dims=1:30)

# Clustering
x=FindNeighbors(x, dims=1:30)
resolutions=c(0.1, 0.2, 0.3,0.4,0.5,0.6,0.7,0.8,0.9)
for (res in resolutions) {
  
  x=FindClusters(x, resolution=res)
  resa=as.character(res)
  
  # Find markers by cluster
  x.markers=FindAllMarkers(x, only.pos=T, min.pct=0.25, logfc.threshold=0.25)
  x.markers %>%
    group_by(cluster) %>%
    top_n(n=10, wt=avg_log2FC) -> top10
  png(paste(prefix, resa, "top10.png", sep="_"), width=10, height=10, units="in", res=300)
  print(DoHeatmap(x, features=top10$gene) + NoLegend())
  dev.off()
  
  # Plot clusters
  png(paste(prefix, resa, "umap.png", sep="_"), width=8, height=8, units="in", res=300)
  print(DimPlot(object = x, reduction = "umap",label=T,label.size = 8))
  dev.off()
  png(paste(prefix, resa, "percent.mt.png", sep="_"), width=4, height=4, units="in", res=300)
  print(VlnPlot(x, features="percent.mt") + NoLegend())
  dev.off()
  png(paste(prefix, resa, "nCount_RNA.png", sep="_"), width=4, height=4, units="in", res=300)
  print(VlnPlot(x, features="nCount_RNA") + NoLegend())
  dev.off()
  png(paste(prefix, resa, "nFeature_RNA.png", sep="_"), width=4, height=4, units="in", res=300)
  print(VlnPlot(x, features="nFeature_RNA") + NoLegend())
  dev.off()
  
  # Plot known retina cell type markers
  features <- c("Cd79a","Cd28","S100a9","Gata3","Gzma","Mcpt4","Cd209a","Trdc","Fscn1","Xcr1","Cd4","Foxp3","Cd8a")
  png(paste(prefix, resa, "known_dot.png", sep="_"), width=10, height=8, units="in", res=300)
  print(DotPlot(x,features = features) + RotatedAxis())
  dev.off()

  # Plot known RPE markers
  features <- c("Cd79a","Cd28","S100a9","Gata3","Gzma","Mcpt4","Cd209a","Trdc","Fscn1","Xcr1","Cd4","Foxp3","Cd8a")
  png(paste(prefix, resa, "known_feature.png", sep="_"), width=10, height=8, units="in", res=300)
  print(FeaturePlot(x, features=features))
	dev.off()
}

saveRDS(x, paste(prefix, "cluster.rds", sep="_"))

