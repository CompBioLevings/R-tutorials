# Load packages
library(doParallel)
library(BiocParallel)
library(grid)
library(gridExtra)
library(labeling)
library(scales)
library(reshape2)
library(ggplot2)
library(splines)
library(knitr)
library(kableExtra)
library(colorspace)
library(RColorBrewer)
library(openxlsx)
library(dplyr)
library(magrittr)
library(pheatmap)
library(ClusterR)

# Read in 'example' data
gene.expression.df <- read.xlsx(xlsxFile = "Example_clustering_data.xlsx", rowNames = T, colNames = T)

# Normalize by mean and standard error using Z-score normalization - may not be necessary for some/most
# Scale scales along columns so we have to 'transpose'/switch columns and rows to scale by genes 
# ACROSS experiments (log2FC = log2(fold-change expression between treatment and control))
log2FC.DGE.df <- t(scale(t(gene.expression.df)))

# Make 'backup' to restore this original df if needed
log2FC.DGE.df.bkup <- log2FC.DGE.df

# Generate gene.matrix for clustering
gene.matrix <- as.matrix(x=as.data.frame(log2FC.DGE.df))
colnames(gene.matrix) <- colnames(log2FC.DGE.df.bkup)

# Make sure there are no 'empty' cells
is.na(gene.matrix) %>% sum()

# Specify the number of clusters you want to 'separate' into using hclust 'cuttree'
# You will have to rerun and do this many times to find the right 'point' potentially
num.centers <- 9

# Get an idea for the scale for the data
min(log2FC.DGE.df)
max(log2FC.DGE.df)

# Set up the heatmap colors - 6.7 is ~max absolute value for Z-norm data
hmcol <- colorRampPalette(rev(brewer.pal(11,"RdBu")))(75)
custombreaks <- seq(from=-2.7, to=2.7, length.out=75)

# Print heatmap made from hierarchical clustering
# pdf("Hclust-heatmap_Z-norm.pdf", width = 3, height = 17, onefile = T)
pheatmap(mat=gene.matrix, cellwidth=7, cellheight=1, color=hmcol, breaks=custombreaks, 
         clustering_distance_rows = "euclidean", clustering_method = "complete", fontsize_col = 7,
         cluster_cols = F, fontsize_row = 1, main = "   Expr Heatmap - Hclust", cutree_rows = num.centers)
# dev.off()

# We may want to show the 'raw' log2FC values instead of normalized values in the heatmap
# To do so, do the clustering first manually, then use the raw values AND clustering object with Pheatmap
rowclust <- hclust(dist(gene.matrix), method = "complete")

# Reset 'clustering matrix' to raw values
gene.matrix <- gene.expression.df

# Redo the scaling for colors to match raw log2FC
# Get an idea for the scale for the data
min(gene.matrix)
max(gene.matrix)

# Set up the heatmap colors - 6.7 is ~max absolute value for Z-norm data
hmcol <- colorRampPalette(rev(brewer.pal(11,"RdBu")))(75)
custombreaks <- seq(from=-5, to=5, length.out=75)

# now plot
# pdf("Hclust-heatmap_Z-norm_rawvals.pdf", width = 3, height = 17, onefile = T)
pheatmap(mat=gene.matrix, cellwidth=7, cellheight=1, color=hmcol, breaks=custombreaks, 
         cluster_rows = rowclust, cluster_cols = F, fontsize_col = 7,fontsize_row = 1, 
         main = "   Expr Heatmap - Hclust", cutree_rows = num.centers)
# dev.off()


## Now an alternate clustering method - Kmeans
# Check for optimal number of clusters
set.seed(123)
optimal.K <- Optimal_Clusters_KMeans(log2FC.DGE.df, max_clusters = 25, plot_clusters = T,
                                     criterion = "dissimilarity")

# Use a Kmeans clustering method, splitting the data into 9 groups (visually try 
# many different numbers, and check for convergence of Optimal_Clusters and visual check on
# a good number of clusters - 9 in this example)
num.centers <- 9
set.seed(123)
Kclust <- KMeans_rcpp(data=log2FC.DGE.df, clusters=num.centers, max_iters = 100, initializer="kmeans++")
ordered.log2FC.clust.df <- gene.expression.df
ordered.log2FC.clust.df$Kmeans.cluster <- Kclust$cluster
order_df <- data.frame(Kmeans.cluster = Kclust$cluster,
                       supp_ordering = rowMeans(ordered.log2FC.clust.df))
ordered.log2FC.clust.df <- ordered.log2FC.clust.df[with(order_df, order(Kmeans.cluster, supp_ordering)), ]

# Get the expression data and sample names for graphing in a heatmp
gene.matrix <- as.matrix(x=as.data.frame(ordered.log2FC.clust.df[,1:9]))
colnames(gene.matrix) <- colnames(ordered.log2FC.clust.df[,1:9])
rownames(gene.matrix) <- rownames(ordered.log2FC.clust.df)

# Get an idea for the scale for the data
min(gene.matrix)
max(gene.matrix)

# Set up the heatmap colors
hmcol <- colorRampPalette(rev(brewer.pal(11,"RdBu")))(75)
custombreaks <- seq(from=-5, to=5, length.out=75)

# Set up 'annotation' column for showing the clusters/boundaries
annotation.df <- data.frame(Cluster = ordered.log2FC.clust.df$Kmeans.cluster)
rownames(annotation.df) <- rownames(gene.matrix)
annotation.colors <- list(Cluster = colorRampPalette(brewer.pal(n = 9, name = "Set1"))(num.centers) %>% setNames(1:num.centers))

# Plot Kmeans with pheatmap - can also directly plot with pheatmap
pheatmap(mat=gene.matrix, cellwidth=7, cellheight=1, color=hmcol, breaks=custombreaks,
         cluster_cols = F, cluster_rows = F, fontsize_col = 6, fontsize_row = 1,
         show_rownames = F, show_colnames = T, main = "Expr Heatmap - Kmeans Clusters",
         annotation_row = annotation.df, annotation_colors = annotation.colors,
         annotation_legend = F, filename = "   Kclust-heatmap_Z-norm_rawvals.pdf",
         width = 3, height = 17
)
