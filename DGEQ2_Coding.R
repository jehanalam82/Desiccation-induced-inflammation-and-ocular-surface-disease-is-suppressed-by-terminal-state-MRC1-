# Load libraries
library(Seurat)
library(DESeq2)
library(EnhancedVolcano)
library(DEGreport)
library(VennDiagram)
library(RColorBrewer)
# Set working directory
wd <- "H:/Manuscripts/DS-Model_ScRNA_Seq_Analysis/DESEQ2/" # MODIFY
setwd(wd)
#Prepare subsets¶
# Load full object with all clusters and all days

filename <- "H:/Manuscripts/DS-Model_ScRNA_Seq_Analysis/DS_merged.rds" # MODIFY
x <-readRDS(filename)

# Edit metadata if needed
x$sample <- factor(x$sample, levels = c("NS","DS1","DS3","DS5","DS10"))
x$days <- as.numeric(x$days)
## Subset cluster of interest
x <- SetIdent(x, value = "RNA_snn_res.0.4") #Set identity as cluster number labels
x.subset <- subset(x, ident = "21") # MODIFY, Subset neutrophil cluster
saveRDS(x.subset, "C21.rds")
# Alternatively run this code to subset all macrophage clusters
# Subset multiple clusters
x <- SetIdent(x, value = "RNA_snn_res.0.4") #Set identity as cluster number labels
x.subset <- subset(x, ident = c("0","2","3","10","11","12")) #Subset all macrophage clusters
saveRDS(x.subset, "MF.rds")
#DESeq2 - DEGs across 10-day period, fold change is calculated as the average of 10-day period 
#The LRT is comparing the full model to the reduced model to identify significant genes. 
#The p-values are determined solely by the difference in deviance between the 'full' and 'reduced' 
#model formula (not log2 fold changes). Essentially the LRT test is testing whether the term(s) 
#removed in the 'reduced' model explains a significant amount of variation in the data.
#If you use these results, I would focus on the p-values and not the log2 fold changes. 
#Genes with high p-values indicate they vary significantly with the variable "days". 
#I included DEG pattern analysis here since I think this is the most effective way to analyze 
#these results without the need of Venn diagrams. I would also recommend plotting violin plots 
#to see the gene expression patterns of the most significant genes.
prefix <- "C20" # MODIFY, filename of subsetted object without suffix
#Run DESeq2¶
# Load subsetted object
x <- readRDS(paste0(prefix,".rds"))

# Prepare metadata for DESeq2 - scale days to 0 to 1
x$days <- as.numeric(x$days)
range01 <- function(x){(x-min(x))/(max(x)-min(x))}
x$days.range <- range01(x$days) # Same as dividing by 10
# Prepare DESeq2 inputs
counts <- x@assays$RNA@counts
metadata <- x@meta.data
design <- ~ days.range
reduced <- ~ 1
# Optional, remove lowly expressed genes before dispersion estimates to improve dispersion fit and reduce computation time/load
# Optional, remove lowly expressed genes before dispersion estimates to improve dispersion fit and reduce computation time/load
dds <- DESeqDataSetFromMatrix(countData = round(counts),
                              colData = metadata,
                              design = design)


min.pct <- 0.1 # MODIFY Remove genes that are expressed in less than 10% of cells in all samples

opts <- lapply(unique(colData(dds)$sample), function(sample) {
  subset <- counts(dds)[,which(colData(dds)$sample == sample)] #Subset by sample
  keep <- rowSums(subset > 0) >= min.pct * ncol(subset) #Keep genes expressed in more than min.pct of that sample
})
keep <- Reduce("|", opts) #Merge passing genes for all samples
print(paste0("Number of genes before filtering: ", dim(counts)[1]))
counts <- counts[which(keep),] #Keep filtered genes
print(paste0("Number of genes after filtering: ", dim(counts)[1]))

# Run DESeq2
dds <- scran::computeSumFactors(dds)
dds <- DESeq(dds, test="LRT", reduced = reduced, useT=TRUE, minmu=1e-6, minReplicatesForReplace=Inf, fitType="glmGamPoi")
saveRDS(dds, paste0(prefix, "_days.rds"))
print(paste0("DESeq2 results saved to ", prefix, "_days.rds"))
#Load and save DESeq2 results
# Load results
dds <- readRDS(paste0(prefix, "_days.rds"))
res <- results(dds, alpha = 0.05, cooksCutoff = FALSE, independentFiltering = FALSE)
res$symbol <- mcols(dds)$symbol
summary(res)
# Add extra gene expression info to result table
res <- as.data.frame(res)

# Add average Seurat normalized expression to result table
data <- x@assays$RNA@data
data <- data[rownames(res),]
opts <- lapply(levels(colData(dds)$sample), function(sample) {
  subset <- data[,which(colData(dds)$sample == sample)] #Subset by sample
  avg <- rowMeans(subset) #Calculate average expression of each gene from Seurat normalized counts
})
names(opts) <- paste0(levels(colData(dds)$sample),".avg.Seurat")
avg.Seurat <- as.data.frame(do.call(cbind, opts))

# Add average DESeq2 normalized expression to result table
opts <- lapply(levels(colData(dds)$sample), function(sample) {
  subset <- counts(dds)[,which(colData(dds)$sample == sample)] #Subset by sample
  avg <- rowMeans(subset) #Calculate average expression of each gene from DESeq2 normalized counts
})
names(opts) <- paste0(levels(colData(dds)$sample),".avg.DESeq2")
avg.DESeq2 <- as.data.frame(do.call(cbind, opts))

# Add percent expression to result table
opts <- lapply(levels(colData(dds)$sample), function(sample) {
  subset <- counts(dds, normalized = TRUE)[,which(colData(dds)$sample == sample)] #Subset normalized counts by sample
  pct <- rowSums(subset > 0) / ncol(subset) #Calculate percent of cells expressing each gene
})
names(opts) <- paste0(levels(colData(dds)$sample),".pct")
pct <- as.data.frame(do.call(cbind, opts))

# Merge with DESeq2 results
res <- cbind(res, avg.Seurat, avg.DESeq2, pct)
head(res)

# Save results to csv file
resOrdered <- res[order(res$padj),]
write.csv(resOrdered, file=paste0(prefix, "_resOrdered.csv")) #Save all genes, ordered by adjusted p-value
print(paste0("Full gene list from DESeq2 saved to ", prefix, "_resOrdered.csv"))
# MODIFY, Set filter cutoffs
padj.cutoff <- 0.0001 # adjusted p-value cutoff
fc.cutoff <- 0        # log2 fold change cutoff
min.pct <- 0.2        # minimum percentage of cells expressed, in at least one sample
# Filter by min.pct
opts <- lapply(unique(colData(dds)$sample), function(sample) {
  subset <- counts(dds)[,which(colData(dds)$sample == sample)]
  keep <- rowSums(subset > 0) >= min.pct * ncol(subset)
})
keep <- Reduce("|", opts)
resSig <- res[which(keep),]

# Filter by log2 fold change and adjusted p-value
resSig <- resSig[which(abs(resSig$log2FoldChange) > fc.cutoff & resSig$padj < padj.cutoff), ]

head(resSig)

# Save results to csv file
resSig <- resSig[order(resSig$padj),]
write.csv(resSig, file=paste0(prefix, "_resSig.csv")) #Save significant genes, ordered by adjusted p-value
print(paste0("Significant gene list from DESeq2 saved to ", prefix, "_resSig.csv"))
#Plot volcano plot¶
# MODIFY, Set filter cutoffs
padj.cutoff <- 0.0001 # adjusted p-value cutoff
fc.cutoff <- 0        # log2 fold change cutoff
min.pct <- 0.2        # minimum percentage of cells expressed, in at least one sample
# Load gene lists
res <- read.csv(file=paste0(prefix, "_resOrdered.csv"), row.names = 1)

# Apply filters again
opts <- lapply(unique(colData(dds)$sample), function(sample) {
  subset <- counts(dds)[,which(colData(dds)$sample == sample)]
  keep <- rowSums(subset > 0) >= min.pct * ncol(subset)
})
keep <- Reduce("|", opts)
resSig <- res[which(keep),]
resSig <- resSig[which(abs(resSig$log2FoldChange) > fc.cutoff & resSig$padj < padj.cutoff), ]

# Plot volcano plot
#png(paste0(prefix, "_volcano.png")), width = 10, height = 10, units = "in", res = 300)
p <- EnhancedVolcano(res,
                     lab = rownames(res),
                     x = 'log2FoldChange',
                     y = 'padj',
                     pCutoff = padj.cutoff,
                     FCcutoff = fc.cutoff,
                     title = paste0(prefix, " volcano plot"),
                     subtitle = paste0(nrow(resSig)," significant genes"))

print(p)
#dev.off()
########################Find DEG patterns¶#################
# Load significant gene list and corresponding subsetted Seurat object
time_key <- "sample"
resSig <- read.csv(file=paste0(prefix, "_resSig.csv"), row.names = 1)
x <- readRDS(paste0(prefix,".rds"))
data <- x@assays$RNA@data
cluster_rlog <- data[rownames(resSig),]

# Find DEG patterns
#png(filename = paste0(prefix, "_clusters.png"), width = 800, height = 600)
clusters <- degPatterns(as.matrix(cluster_rlog), metadata = dds@colData, time = time_key, col = NULL, plot = T, minc = 1)
#dev.off()

# Check cluster patterns with other cutoffs
# png(filename = paste0(prefix, "_clusters_benchmarking.png"), width = 800, height = 600)
p <- clusters$benchmarking
print(p)
# dev.off()

# Save degPatterns object
saveRDS(clusters, paste0(prefix,"_clusters.rds"))

# Add cluster information to resSig table and save to csv file
resSig <- cbind(resSig, clusters$df$cluster)
write.csv(resSig, file=paste0(prefix, "_resSig_clusters.csv")) #Save significant genes
###############################################################################################
# Optional, plot other cutoffs from benchmarking and save clusters
cutoff <- "cutoff0.64" # MODIFY, choose cutoff from benchmarking plot

# Plot clusters corresponding to chosen cutoff
degPlotCluster(clusters$normalized, time_key, min_genes = 1, cluster_column = cutoff)

# Add cluster information to resSig table and save to csv file
cluster <- clusters$normalized[!duplicated(clusters$normalized$genes),cutoff]
resSig <- cbind(resSig, cluster)
write.csv(resSig, file=paste0(prefix, "_resSig_", cutoff,".csv")) #Save significant genes
##################################################################################################
## Misc plotting with Seurat
####################################################################################################
# Various plots of most significant genes
resSig <- read.csv(file=paste0(prefix, "_resSig.csv"), row.names = 1)
resSig <- resSig[order(resSig$padj),]
resSig$gene <- rownames(resSig)

x <- readRDS(paste0(prefix,".rds"))

# Violin plot with or without points
VlnPlot(x, features = resSig$gene[1:10], group.by = "sample")
VlnPlot(x, features = resSig$gene[1:10], group.by = "sample", pt.size = 0)

# Dot plot
DotPlot(x, features = resSig$gene[1:10], group.by = "sample") + coord_flip()

# Average expression heatmap
x.avg <- AverageExpression(x, group.by = "sample", return.seurat = TRUE)
x.avg$sample <- factor(rownames(x.avg@meta.data), levels = c("NS","DS1","DS3","DS5","DS10"))
DoHeatmap(x.avg, features = resSig$gene[1:100], group.by = "sample", draw.lines = FALSE)
#################################################################################################
#DEGs between each sample, for venn diagram
# (If you want a venn diagram to compare between samples, I can think of a few different methods/scenarios.
#Option 1, you use "sample" (categorical variable) instead of "days" (continuous variable) as the design formula 
#for the LRT test. The test will give you a global p-value across all samples, which basically tells you how significant 
#the gene varies among all samples. The log fold change given in the results will be pairwise log fold change.
#Option 2, you can do pairwise comparisons using the LRT test, where you subset the two samples you want to compare, 
#and then repeat option 1. This way gives a separate p-value and fold change value for each comparison. The fold change 
#values are expected to be the same as in option 1.
#Option 3, you can do pairwise comparisons using the Wald test. This is built into DESeq2 so you do not need to subset 
#the samples for each comparison. This gives you a separate p-value and fold change value for each comparison. The fold 
#change values should be the same or similar to option 1 and 2. https://support.bioconductor.org/p/100333/
#Option 4, because you are ignoring the time information in the variable "days" with this approach, I think the simpler 
#and more tried-and-true approach would be to just use the FindMarkers function in Seurat using the wilcoxon rank sum test 
#(which is known to work well with single-cell data, vs. DESeq2 is originally designed for bulk RNA-seq data).
#I think option 4 makes the most sense, so that is the code I include below. If you want code for one of the other options, 
#I can also write it.
# MODIFY, Set filter cutoffs

padj.cutoff <- 0.05 # adjusted p-value cutoff
fc.cutoff <- 1        # log2 fold change cutoff
min.pct <- 0.1        # minimum percentage of cells expressed, in at least one sample

# Specify comparisons - comment/uncomment one below or change accordingly
comparisons <- list(c("NS", "DS1"),
                    c("DS1", "DS3"),
                    c("DS3", "DS5"),
                    c("DS5", "DS10"))
comparisons <- list(c("NS", "DS1"),
                    c("NS", "DS3"),
                    c("NS", "DS5"),
                    c("NS", "DS10"))

names(comparisons) <- lapply(comparisons, function(comparison) {paste0(comparison[1],".vs.",comparison[2])})
# Load subsetted object
x <- readRDS(paste0(prefix,".rds"))
# Run FindMarkers
deg.list <- lapply(comparisons, function(comparison) {
  deg <- FindMarkers(x, ident.1 = comparison[1], ident.2 = comparison[2], group.by = "sample",
                     min.pct = min.pct, logfc.threshold = fc.cutoff)
  deg <- deg[which(deg$p_val_adj < padj.cutoff),] #Filter by adjusted p-value
  
  # Save pairwise comparison results
  write.csv(deg, paste0(prefix,"_",comparison[1],".vs.",comparison[2],".csv"))
  print(paste0("FindMarkers results saved to ", prefix,"_",comparison[1],".vs.",comparison[2],".csv"))
  
  return(deg)
})
# Extract genes for venn diagram
deg.genes <- lapply(deg.list, function(comparison) {
  return(rownames(comparison))
})
deg.genes
# Plot venn diagram
futile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger")
grid.newpage() 
myCol <- brewer.pal(4, "Pastel2")
venn <- venn.diagram(
  x = deg.genes,
  #filename = paste0(prefix, "_venn.png"), #Uncomment to save
  filename = NULL, #Comment to save
  output=TRUE,
  fill = myCol
)

grid.draw(venn)
# To obtain intersecting genes, MOIDFY code below
intersect(deg.genes[['NS.vs.DS5']], deg.genes[['NS.vs.DS10']])

#############################DEGs for macrophages across pseudotime or stage#############################################
prefix <- "MF" # MODIFY, filename of subsetted object without suffix
cellrank_filename <- "H:/Manuscripts/DS-Model_ScRNA_Seq_Analysis/DESEQ2/MF_adata_processed_deterministic_terminal_states.csv"
###########################Run DESeq2 by pseudotime#############################
# Load subsetted object
x <- readRDS(paste0(prefix,".rds"))

# Add cellrank metadata
cellrank <- read.csv(cellrank_filename, row.names = 1)
for (col in c("palantir_pseudotime", "absorption_prob_3", "absorption_prob_2", "absorption_prob_0", "stage")) {
  x[[col]] <- cellrank[[col]]
}

# Subset macrophage lineage
x <- subset(x, subset = absorption_prob_0 > 0.6)

# Prepare metadata for DESeq2
x$palantir_pseudotime <- as.numeric(x$palantir_pseudotime)
range01 <- function(x){(x-min(x))/(max(x)-min(x))}
x$pseudotime <- range01(x$palantir_pseudotime) # Rescale pseudotime after subsetting
x$stage <- plyr::mapvalues(x$stage, from = c("stage 1", "stage 2", "stage 3"), to = c("1", "2", "3"))
x$stage <- factor(x$stage, levels = c("1", "2", "3"))
# Prepare DESeq2 inputs
counts <- x@assays$RNA@counts
metadata <- x@meta.data
design <- ~ pseudotime
reduced <- ~ 1
dds <- DESeqDataSetFromMatrix(countData = round(counts),
                              colData = metadata,
                              design = design)
# Optional, remove lowly expressed genes before dispersion estimates to improve dispersion fit and reduce computation time/load
min.pct <- 0.1 # MODIFY Remove genes that are expressed in less than 10% of cells in all samples

opts <- lapply(unique(metadata$sample), function(sample) {
  subset <- counts(dds)[,which(metadata$sample == sample)] #Subset by sample
  keep <- rowSums(subset > 0) >= min.pct * ncol(subset) #Keep genes expressed in more than min.pct of that sample
})
keep <- Reduce("|", opts) #Merge passing genes for all samples
print(paste0("Number of genes before filtering: ", dim(counts)[1]))
counts <- counts[which(keep),] #Keep filtered genes
print(paste0("Number of genes after filtering: ", dim(counts)[1]))
# Run DESeq2
dds <- scran::computeSumFactors(dds)
dds <- DESeq(dds, test="LRT", reduced = reduced, useT=TRUE, minmu=1e-6, minReplicatesForReplace=Inf, fitType="glmGamPoi")
saveRDS(dds, paste0(prefix, "_pseudotime.rds"))
print(paste0("DESeq2 results saved to ", prefix, "_pseudotime.rds"))
################################################################################################################################################
#Load and save DESeq2 results
# Load results
dds <- readRDS(paste0(prefix, "_pseudotime.rds"))
res <- results(dds, alpha = 0.05, cooksCutoff = FALSE, independentFiltering = FALSE)
res$symbol <- mcols(dds)$symbol
summary(res)
# Add extra gene expression info to result table
res <- as.data.frame(res)

# Add average Seurat normalized expression to result table
data <- x@assays$RNA@data
data <- data[rownames(res),]
opts <- lapply(levels(colData(dds)$sample), function(sample) {
  subset <- data[,which(colData(dds)$sample == sample)] #Subset by sample
  avg <- rowMeans(subset) #Calculate average expression of each gene from Seurat normalized counts
})
names(opts) <- paste0(levels(colData(dds)$sample),".avg.Seurat")
avg.Seurat <- as.data.frame(do.call(cbind, opts))

# Add average DESeq2 normalized expression to result table
opts <- lapply(levels(colData(dds)$sample), function(sample) {
  subset <- counts(dds)[,which(colData(dds)$sample == sample)] #Subset by sample
  avg <- rowMeans(subset) #Calculate average expression of each gene from DESeq2 normalized counts
})
names(opts) <- paste0(levels(colData(dds)$sample),".avg.DESeq2")
avg.DESeq2 <- as.data.frame(do.call(cbind, opts))

# Add percent expression to result table
opts <- lapply(levels(colData(dds)$sample), function(sample) {
  subset <- counts(dds, normalized = TRUE)[,which(colData(dds)$sample == sample)] #Subset normalized counts by sample
  pct <- rowSums(subset > 0) / ncol(subset) #Calculate percent of cells expressing each gene
})
names(opts) <- paste0(levels(colData(dds)$sample),".pct")
pct <- as.data.frame(do.call(cbind, opts))

# Merge with DESeq2 results
res <- cbind(res, avg.Seurat, avg.DESeq2, pct)
head(res)

# Save results to csv file
resOrdered <- res[order(res$padj),]
write.csv(resOrdered, file=paste0(prefix, "_resOrdered.csv")) #Save all genes, ordered by adjusted p-value
print(paste0("Full gene list from DESeq2 saved to ", prefix, "_resOrdered.csv"))
# MODIFY, Set filter cutoffs
padj.cutoff <- 0.0001 # adjusted p-value cutoff
fc.cutoff <- 0        # log2 fold change cutoff
min.pct <- 0.2        # minimum percentage of cells expressed, in at least one sample
# Filter by min.pct
opts <- lapply(unique(colData(dds)$sample), function(sample) {
  subset <- counts(dds)[,which(colData(dds)$sample == sample)]
  keep <- rowSums(subset > 0) >= min.pct * ncol(subset)
})
keep <- Reduce("|", opts)
resSig <- res[which(keep),]

# Filter by log2 fold change and adjusted p-value
resSig <- resSig[which(abs(resSig$log2FoldChange) > fc.cutoff & resSig$padj < padj.cutoff), ]

head(resSig)

# Save results to csv file
resSig <- resSig[order(resSig$padj),]
write.csv(resSig, file=paste0(prefix, "_resSig.csv")) #Save significant genes, ordered by adjusted p-value
print(paste0("Significant gene list from DESeq2 saved to ", prefix, "_resSig.csv"))
#Plot volcano plot¶
# MODIFY, Set filter cutoffs
padj.cutoff <- 0.0001 # adjusted p-value cutoff
fc.cutoff <- 0        # log2 fold change cutoff
min.pct <- 0.2        # minimum percentage of cells expressed, in at least one sample
# Load gene lists
res <- read.csv(file=paste0(prefix, "_resOrdered.csv"), row.names = 1)

# Apply filters again
opts <- lapply(unique(colData(dds)$sample), function(sample) {
  subset <- counts(dds)[,which(colData(dds)$sample == sample)]
  keep <- rowSums(subset > 0) >= min.pct * ncol(subset)
})
keep <- Reduce("|", opts)
resSig <- res[which(keep),]
resSig <- resSig[which(abs(resSig$log2FoldChange) > fc.cutoff & resSig$padj < padj.cutoff), ]

# Plot volcano plot
#png(paste0(prefix, "_volcano.png")), width = 10, height = 10, units = "in", res = 300)
p <- EnhancedVolcano(res,
                     lab = rownames(res),
                     x = 'log2FoldChange',
                     y = 'padj',
                     pCutoff = padj.cutoff,
                     FCcutoff = fc.cutoff,
                     title = paste0(prefix, " volcano plot"),
                     subtitle = paste0(nrow(resSig)," significant genes"))

print(p)
#dev.off()
########################Find DEG patterns¶#################
# Load significant gene list and corresponding subsetted Seurat object
time_key <- "sample"
resSig <- read.csv(file=paste0(prefix, "_resSig.csv"), row.names = 1)
x <- readRDS(paste0(prefix,".rds"))
data <- x@assays$RNA@data
cluster_rlog <- data[rownames(resSig),]

# Find DEG patterns
#png(filename = paste0(prefix, "_clusters.png"), width = 800, height = 600)
clusters <- degPatterns(as.matrix(cluster_rlog), metadata = dds@colData, time = time_key, col = NULL, plot = T, minc = 1)
#dev.off()

# Check cluster patterns with other cutoffs
# png(filename = paste0(prefix, "_clusters_benchmarking.png"), width = 800, height = 600)
p <- clusters$benchmarking
print(p)
# dev.off()

# Save degPatterns object
saveRDS(clusters, paste0(prefix,"_clusters.rds"))

# Add cluster information to resSig table and save to csv file
resSig <- cbind(resSig, clusters$df$cluster)
write.csv(resSig, file=paste0(prefix, "_resSig_clusters.csv")) #Save significant genes
###############################################################################################
# Optional, plot other cutoffs from benchmarking and save clusters
cutoff <- "cutoff0.64" # MODIFY, choose cutoff from benchmarking plot

# Plot clusters corresponding to chosen cutoff
degPlotCluster(clusters$normalized, time_key, min_genes = 1, cluster_column = cutoff)

# Add cluster information to resSig table and save to csv file
cluster <- clusters$normalized[!duplicated(clusters$normalized$genes),cutoff]
resSig <- cbind(resSig, cluster)
write.csv(resSig, file=paste0(prefix, "_resSig_", cutoff,".csv")) #Save significant genes
##################################################################################################
## Misc plotting with Seurat
####################################################################################################
# Various plots of most significant genes
resSig <- read.csv(file=paste0(prefix, "_resSig.csv"), row.names = 1)
resSig <- resSig[order(resSig$padj),]
resSig$gene <- rownames(resSig)

x <- readRDS(paste0(prefix,".rds"))

# Violin plot with or without points
VlnPlot(x, features = resSig$gene[1:10], group.by = "sample")
VlnPlot(x, features = resSig$gene[1:10], group.by = "sample", pt.size = 0)

# Dot plot
DotPlot(x, features = resSig$gene[1:10], group.by = "sample") + coord_flip()

# Average expression heatmap
x.avg <- AverageExpression(x, group.by = "sample", return.seurat = TRUE)
x.avg$sample <- factor(rownames(x.avg@meta.data), levels = c("NS","DS1","DS3","DS5","DS10"))
DoHeatmap(x.avg, features = resSig$gene[1:100], group.by = "sample", draw.lines = FALSE)
#################################################################################################
#DEGs between each sample, for venn diagram
# (If you want a venn diagram to compare between samples, I can think of a few different methods/scenarios.
#Option 1, you use "sample" (categorical variable) instead of "days" (continuous variable) as the design formula 
#for the LRT test. The test will give you a global p-value across all samples, which basically tells you how significant 
#the gene varies among all samples. The log fold change given in the results will be pairwise log fold change.
#Option 2, you can do pairwise comparisons using the LRT test, where you subset the two samples you want to compare, 
#and then repeat option 1. This way gives a separate p-value and fold change value for each comparison. The fold change 
#values are expected to be the same as in option 1.
#Option 3, you can do pairwise comparisons using the Wald test. This is built into DESeq2 so you do not need to subset 
#the samples for each comparison. This gives you a separate p-value and fold change value for each comparison. The fold 
#change values should be the same or similar to option 1 and 2. https://support.bioconductor.org/p/100333/
#Option 4, because you are ignoring the time information in the variable "days" with this approach, I think the simpler 
#and more tried-and-true approach would be to just use the FindMarkers function in Seurat using the wilcoxon rank sum test 
#(which is known to work well with single-cell data, vs. DESeq2 is originally designed for bulk RNA-seq data).
#I think option 4 makes the most sense, so that is the code I include below. If you want code for one of the other options, 
#I can also write it.
# MODIFY, Set filter cutoffs

padj.cutoff <- 0.05 # adjusted p-value cutoff
fc.cutoff <- 1        # log2 fold change cutoff
min.pct <- 0.1        # minimum percentage of cells expressed, in at least one sample

# Specify comparisons - comment/uncomment one below or change accordingly
comparisons <- list(c("NS", "DS1"),
                    c("DS1", "DS3"),
                    c("DS3", "DS5"),
                    c("DS5", "DS10"))
comparisons <- list(c("NS", "DS1"),
                    c("NS", "DS3"),
                    c("NS", "DS5"),
                    c("NS", "DS10"))

names(comparisons) <- lapply(comparisons, function(comparison) {paste0(comparison[1],".vs.",comparison[2])})
# Load subsetted object
x <- readRDS(paste0(prefix,".rds"))
# Run FindMarkers
deg.list <- lapply(comparisons, function(comparison) {
  deg <- FindMarkers(x, ident.1 = comparison[1], ident.2 = comparison[2], group.by = "sample",
                     min.pct = min.pct, logfc.threshold = fc.cutoff)
  deg <- deg[which(deg$p_val_adj < padj.cutoff),] #Filter by adjusted p-value
  
  # Save pairwise comparison results
  write.csv(deg, paste0(prefix,"_",comparison[1],".vs.",comparison[2],".csv"))
  print(paste0("FindMarkers results saved to ", prefix,"_",comparison[1],".vs.",comparison[2],".csv"))
  
  return(deg)
})
# Extract genes for venn diagram
deg.genes <- lapply(deg.list, function(comparison) {
  return(rownames(comparison))
})
deg.genes
# Plot venn diagram
futile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger")
grid.newpage() 
myCol <- brewer.pal(4, "Pastel2")
venn <- venn.diagram(
  x = deg.genes,
  #filename = paste0(prefix, "_venn.png"), #Uncomment to save
  filename = NULL, #Comment to save
  output=TRUE,
  fill = myCol
)

grid.draw(venn)
# To obtain intersecting genes, MOIDFY code below
intersect(deg.genes[['NS.vs.DS5']], deg.genes[['NS.vs.DS10']])

