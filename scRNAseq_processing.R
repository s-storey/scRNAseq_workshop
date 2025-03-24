###############################################################
# Title: Comprehensive Seurat Workflow Script
# Author: [Your Name]
# Date: [Workshop Date]
# 
# Description:
#   This script provides a hands-on workflow for single-cell RNA-seq
#   data analysis using the Seurat R package. We will go through:
#     1. Data import (using an example dataset)
#     2. Quality control (QC) and filtering
#     3. Normalization
#     4. Identifying highly variable features
#     5. Scaling and dimensional reduction
#     6. Clustering
#     7. Non-linear dimension reduction (UMAP)
#     8. Differential expression (DE)
#     9. Visualization and result exporting
#
#   Detailed comments are provided at each step explaining:
#     - Why each step is performed
#     - How to adjust parameters based on your data
#     - Potential pitfalls and best practices
###############################################################

##########################
# 0) Environment Setup
##########################

# Install (if needed) and load the required packages.
# In practice, you might need to install them beforehand.
# install.packages("Seurat")         # For single-cell analysis
# install.packages("patchwork")      # For combining plots
# install.packages("dplyr")          # For data manipulation
# install.packages("ggplot2")        # For plotting

library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)

# For demonstration, we will load the pbmc3k dataset (3k Peripheral Blood Mononuclear Cells).
# If you have your own data, you can skip this part or adapt it.


##########################
# 1) Data Import
##########################

# Move your downloaded "filtered_gene_bc_matrices" file to your working directory
# In RStudio, find your directory, click 'More' and select 'Set As Working Directory'
# Click into 'filtered_gene_bc_matrices' folder, then click the hg19 folder
# After verifying your files are there, click 'More' then 'Copy Folder Path to Clipboard'

# paste the path here
data_dir <- "your/path/to/hg9"
# Load the PBMC dataset
# if we had a .h5 file, we would use Read10X_h5 instead
pbmc.data <- Read10X(data.dir = data_dir)

# Initialize the Seurat object with the raw (non-normalized data).
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", min.cells = 3, min.features = 200)
# Explanation:
# - project: a name for the dataset/project.
# - min.cells, min.features: filter out genes/cells that have very 
#   low coverage to remove possible technical artifacts or empty droplets.
# - min.features 

# Lets take a peak at our data
# Below are some options to visualize this data
pbmc # gives information about your object. 
# An object of class Seurat 
# 13714 features across 2700 samples within 1 assay 
# Active assay: RNA (13714 features, 0 variable features)
# 1 layer present: counts

#### Slots
# Common slots for a seurat object include: 
# Assays
# meta.data
# reductions

# Assays hold expression data and processed versions of this data 
# Layers are representations of expression data within an assay
# counts: raw count data
# data: the normalized data
# scale.data: the scaled data usually prepared for downstream analyses like PCA
# Metadata is cell-level data with cells as row names and metadata variables as columns
# Reductions hold dimensionality reduction results

# Click the magnifying glass next to your object in environment
# Can also access these slots using the @ operator, eg. 

head(pbmc@meta.data) # printed in console 
View(pbmc@meta.data) # a spreadsheet like view



# For this demonstration, we have the pbmc3k dataset in the variable `pbmc`.
# Let's rename it to "sc_object" for consistency.

sc_object <- pbmc

#################################################
# 2) Quality Control (QC) & Initial Filtering
#################################################

# Typically, we look at:
# 1. Number of detected genes (features) per cell
# 2. Number of counts (UMIs) per cell
# 3. Mitochondrial gene percentage
# 
# Cells with few detected genes likely represent empty droplets or dying cells.
# Cells with extremely high gene counts may be doublets or multiplets.
# High mitochondrial gene content can indicate stressed or dying cells.

# First, let's visualize the QC metrics that are already part of pbmc3k metadata.
head(sc_object@meta.data)


# We see:
#                     orig.ident nCount_RNA nFeature_RNA
# AAACATACAACCAC-1     pbmc3k       2419          779
# AAACATTGAGCTAC-1     pbmc3k       4903         1352
# AAACATTGATCAGC-1     pbmc3k       3147         1129
# AAACCGTGCTTCCG-1     pbmc3k       2639          960
# AAACCGTGTATGCG-1     pbmc3k        980          521
# AAACGCACTGGTAC-1     pbmc3k       2163          781

# The dataset doesn't have a mitochondrial content metadata column, so we can conclude
# that mitochondrial content hasn't been set.
# mitochondrial genes are always prefaced with 'MT-'
# keep in mind gene naming convention however, eg zebrafish genes would be lowercase 'mt-'

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
sc_object[["percent.mt"]] <- PercentageFeatureSet(sc_object, pattern = "^MT-") # mt- for zebrafish, Mt- for mouse?
#how to show metadata columns 


head(sc_object@meta.data)
# we see our metadata column

# We can plot some distributions to check the ranges:
VlnPlot(sc_object, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# Filtering (thresholds):
# - Cells with fewer than 200 features are typically considered low quality.
# - Cells with extremely high gene counts may be doublets or multiplets.
# - Cells with high mitochondrial percentage (eg., >10%, >5% ) are often removed.
#   The threshold depends on your tissue type and typical mitochondrial content.

# Let's define thresholds (these can be adjusted based on your data specifics)
# and create a subset of cells that pass these QC filters:

sc_object <- subset(sc_object, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

# Lets look at our data again
VlnPlot(sc_object, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# Why adjust defaults?
# - If your dataset is known to have high mitochondrial content (e.g., certain tissues),
#   you might choose a higher max_mito threshold (like 10% or 15%).
# - If your dataset has generally lower sequencing depth, you might choose a lower 
#   min_features threshold.

# Check how many cells remain after filtering:
dim(pbmc) ## this was our unaltered data - 2700 cells
# [1] 13714  2700

dim(sc_object) ## this is our data subset - 2638
#[1] 13714  2638

##########################
# 3) Normalization
##########################

# Seurat offers multiple normalization methods. The most common approach is "LogNormalize".
# Alternatively, you can use "SCTransform" which models technical noise and can improve 
# downstream analysis. Here we demonstrate "LogNormalize":

sc_object <- NormalizeData(sc_object, 
                           normalization.method = "LogNormalize", 
                           scale.factor = 10000)

# Explanation:
# - "LogNormalize": Each cell is normalized by the total expression,
#   multiplied by a scale factor (10,000 by default), and then log-transformed.
# - scale.factor: can be adjusted if you want to keep expression within a certain range.
#   The default 10,000 is usually fine.

# LogNormalize relies on an assumption that each cell originally contains the same number of RNA molecules
# For a method that does not rely on this assumption, use SCTransform
# Seurat vignette: https://satijalab.org/seurat/articles/sctransform_vignette

# For SCTransform (which can often yield better results for complex data), you would do:
# sc_object <- SCTransform(sc_object, vars.to.regress = "percent.mt", verbose = FALSE)
# 'vars.to.regress = "percent.mt"' means that we're removing mitochondrial
# content as a factor that might cause unwanted source of noise. We could also regress out for 
# example, cell cycle genes
# 
# Note that the SCTransform approach changes the recommended workflow for 
# downstream steps. SCTransform replaces the need to run NormalizeData, 
# FindVariableFeatures, or ScaleData


# Ultimately, Normalization makes sure that expression values are comparable
# across cells by correcting for technical variation

##########################
# 4) Identify Highly Variable Features (HVG) -- USING THE LOGNORMALIZE TECHNIQUE
##########################

# Identifying HVGs focuses analysis on the most informative genes.
# For "LogNormalize" data, we typically do:

sc_object <- FindVariableFeatures(sc_object, 
                                  selection.method = "vst", 
                                  nfeatures = 2000)

# Explanation:
# - nfeatures: Number of top variable genes to focus on. 2,000 is a common default.
#   You can adjust if you want to be more or less stringent (e.g., 3,000 or 1,500).
# These genes will be used to drive our PCA later on 

# We can visualize the top variable features:
top10 <- head(VariableFeatures(sc_object), 10)
VariableFeaturePlot(sc_object) %>% LabelPoints(plot = ., points = top10, repel = TRUE)

##########################
# 5) Scaling the Data -- ONLY IF USING THE LOGNORMALIZE TECHNIQUE
##########################

# Explanations
# We scale the data so that each gene has mean = 0 and variance = 1, 
# allowing PCA to focus on relative differences rather than absolute gene expression levels.
# Different genes naturally vary in their expression levels and variability. Without
# scaling, genes that are highly expressed or have large variability can dominate the analysis. 
# This allows comparable 'units' of gene expression for downstream analyses like PCA, which 
# assumes that the data are on comparable scales -- eg. gene variance is measured in the same 'unit'
# This allows PCA to correctly capture the structure and major sources of variation in the dataset

all.genes <- rownames(sc_object)
sc_object <- ScaleData(sc_object, features = all.genes)

# Further info. 
# - If dataset is large, scaling all genes can be time-consuming. 
#   Often, you scale only your variable features to save time.
# - You can also regress out unwanted sources of variation (e.g., percent.mt, cell cycle scores).
#   e.g., ScaleData(sc_object, vars.to.regress = "percent.mt")


# Ultimately, scaling standardizes each gene's expression across cells to that 
# evey gene has equal weight in analyses like PCA or clustering


##########################
# 6) Principal Component Analysis (PCA) -- WHETHER YOU USED SCTRANSFORM OR LOGNORMALIZE
##########################

# The goal is to visualize and group cells in a way that captures major patterns 
# (which cells are similar or different) without being overwhelmed by noise 
# or irrelevant variation


# PCA is an unsupervised dimension reduction technique to make sense of the complex, high-
# dimensional data generated by scRNAseq. 
# Each cell is measured across thousands of genes, creating a dataset with a huge number
# of dimensions. PCA reduces this complexity by identifying a smaller set of new variables
# called principal components that capture the most important patterns in the data. 
# Documentation: https://satijalab.org/seurat/reference/runpca

sc_object <- RunPCA(sc_object, features = VariableFeatures(object = sc_object))

# Let's examine the results for the first 5 PCs :
print(sc_object[["pca"]], dims = 1:5, nfeatures = 5)  # top genes driving each PC

# Another method to visualise
VizDimLoadings(sc_object, dims = 1:5, reduction = "pca")

# We can visualize the variance explained by each PC to decide how many to use.
# Documentation:  https://satijalab.org/seurat/reference/elbowplot
ElbowPlot(sc_object)

# Explanation:
# - The "elbow" in the plot typically indicates a reasonable cutoff for significant PCs.
# - You might choose PCs 1 through 10 (for example), but the exact number depends on 
#   the dataset variance structure. -- in this example I am choosing 10, but 9 would be fine too

## See also JackStraw for statistics-based dimensionality determination

##########################
# 7) Clustering
##########################

# Seurat clustering typically involves two steps:
# 1. Construct a Shared Nearest Neighbor (SNN) graph for a chosen number of PCs using FindNeighbors
# the default number of neighbours is 20
# Documentation: https://satijalab.org/seurat/reference/findneighbors
sc_object <- FindNeighbors(sc_object, dims = 1:10)

# 2. Using an algorithm to cluster cells based on the snn graph with FindClusters
# several different algorithms to use here, seurat defaults to Louvain
# Documentation: https://satijalab.org/seurat/reference/findclusters
sc_object <- FindClusters(sc_object, resolution = 0.5)

# Explanation:
# - dims: we choose how many PCs to use in constructing the graph to capture appropriate variation.
#   Typically this is informed by ElbowPlot or other statistical method. 

# - resolution: controls the granularity of the clustering. Higher resolution = more clusters.
#   E.g., 0.8 might yield more clusters than 0.5.

# Check the identified clusters:
head(Idents(sc_object))
# Levels: 0 1 2 3 4 5 6 7 8
# Looks like this created 9 clusters 

##########################
# 8) Non-linear Dimensional Reduction (UMAP)
##########################

# UMAP or t-SNE can be used for visualization of the high-dimensional data 
# in 2D or 3D space. 
# UMAP works to maintain both local relationships (distance between cells) and 
# global structure (how different clusters of cell relate to eachother). The distances
# between clusters however are typically best suited for intuitive interpretation
# rather than actual analysis. 
# tSNE is focused primarily on presenving local relationships, but typically distorts
# global relationships. . 
# UMAP is generally faster and preserves global structure well.
# Documentation: https://satijalab.org/seurat/reference/runumap

sc_object <- RunUMAP(sc_object, dims = 1:10)

# Explanation:
# - dims: same number of PCs used for clustering is typically used for UMAP.
# - You can adjust "n.neighbors", "min.dist" in the RunUMAP() function 
#   for different visualization granularities.

# Let's see the 2D UMAP plot with cluster identities:
DimPlot(sc_object, reduction = "umap", label = TRUE)
# we can now visualize our 9 clusters

##########################
# 9) Differential Expression & Marker Identification
##########################

# To find differentially expressed genes (markers) for each cluster:
cluster_markers <- FindAllMarkers(sc_object, 
                                  only.pos = TRUE, ## change this to false if you want to include downregulated genes
                                  min.pct = 0.25, 
                                  logfc.threshold = 0.25)
# Explanation:
# - only.pos = TRUE: only return positive markers for each cluster.
# - min.pct = 0.25: only test genes detected in at least 25% of cells in that cluster.
# - logfc.threshold = 0.25: only return genes with log2 fold-change > 0.25.

cluster_markers
# p_val: unadjusted p-value to determine if a gene is differentially expressed
# avg_log2FC: log2 fold change in expression between a group of interest and
# the reference group. 
# pct.1: the percentage of cells within the cluster of interest that express the gene
# pct.2: the percentage of cells in the comparison group (eg. all cells outside of the 
# cluster of interest) that express the gene. 
### Comparing pct.1 and pct.2  provides context for whether the gene is uniquely or
# preferentially expressed in one group
# p_val_adj: Boneferroni-djusted p-value after correcting for multiple comparisons
# helps control for the increased chance of false positives when testing many genes simultaenously 


# Let's see the top markers per cluster: 
cluster_markers %>% # the %>% operator is part of the dplyr package
  group_by(cluster) %>% 
  slice_max(n = 5, order_by = avg_log2FC) %>% 
  print(n=45)



write.csv(cluster_markers, "cluster_markers.csv")


# If you want to do pairwise differential expression (e.g., cluster 1 vs cluster 2):
# markers_1_vs_2 <- FindMarkers(sc_object, 
#                               ident.1 = 1, 
#                               ident.2 = 2, 
#                               min.pct = 0.25)

##########################
# 10) Visualization & Figure Creation
##########################

# Seurat provides many plotting functions. Some commonly used:
# - DimPlot() for clusters, metadata visualization
DimPlot(sc_object, reduction = "umap", label = TRUE) 

# - VlnPlot() for gene expression distribution
VlnPlot(sc_object, features = "CD3D")

# - FeaturePlot() for gene expression over UMAP
FeaturePlot(sc_object, features = "CD3D")

# - RidgePlot() or DotPlot() for gene expression in a cluster context
DotPlot(sc_object, features = "CD3D")

# - DoHeatmap() for a heatmap of gene expression.
# Lets make a list of the top 5 genes per cluster first:
cluster_markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 5) %>%
  ungroup() -> top5

DoHeatmap(sc_object, features = top5$gene) + NoLegend()

# Example: Combine multiple plots with patchwork
p1 <- DimPlot(sc_object, reduction = "umap", label = TRUE) + ggtitle("UMAP Clusters")
p2 <- VlnPlot(sc_object, features = "CD3D")
p1 
p2
p1 + p2

# For a publication-quality figure, you might use: -- We'll get into this in a later module
# 1. Customize color palettes
# 2. Adjust font sizes, legends, and layout
# 3. Use ggplot2 theming functions (theme_minimal, theme_classic, etc.)
# 4. Export high-resolution images

# Example: Saving a plot
ggsave("UMAP_with_clusters.png", plot = p2, width = 6, height = 5, dpi = 300)

##########################
# 11) Saving & Exporting Data
##########################

# It's often useful to save your Seurat object for future analysis:
saveRDS(sc_object, file = "sc_object_postQC.rds")

# Later, you can load it with:
# sc_object <- readRDS("sc_object_postQC.rds")


##########################
# 12) Advanced Techniques (Brief Intro)
##########################

# 1) Integration of multiple datasets:
#    - If you have multiple samples or batches, you can use Seurat's integration workflow
#      (FindIntegrationAnchors, IntegrateData) or SCTransform-based integration.
#
# 2) Trajectory / Pseudotime Analysis:
#    - Tools like Monocle, Slingshot, or the Seurat Wrappers can help reconstruct
#      lineage relationships among cells if you suspect a continuous transition.
#
# 3) Accessing Public scRNA-seq Data:
#    - Repositories like GEO, ArrayExpress, or Single Cell Portal often provide
#      count matrices or processed data. The same pipeline can be applied after import.
#
# 4) Additional QC steps (Doublet detection):
#    - Tools like DoubletFinder or Scrublet can be integrated before clustering.

##########################
# 13) Session Info & Troubleshooting
##########################

# If you encounter errors, it’s often helpful to check your R session info:
sessionInfo()

# Common pitfalls:
# - Ensure packages are up-to-date.
# - Keep an eye on memory usage for large datasets.
# - If clustering yields too many or too few clusters, adjust the resolution parameter.
# - If your PC selection doesn't capture biology well, consider changing dims in FindNeighbors/UMAP.

###############################################################
# 
# At this point, you have:
#   1. Cleaned your data (QC)
#   2. Normalized & identified variable features
#   3. Performed PCA, clustering, and visualization with UMAP
#   4. Identified cluster-specific markers and performed differential expression
#   5. Learned how to save your results and produce publication-quality figures
# 
###############################################################