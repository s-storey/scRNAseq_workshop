###############################################################
# Title: Seurat Workflow with Annotation, DE, and Subsetting
# Author: [Your Name]
# Date: [Workshop Date]
#
# Description:
#   This script demonstrates:
#     1. Cell type annotation using known marker genes
#     2. Differential expression (cluster vs. cluster or condition vs. condition)
#     3. Subsetting specific cell types of interest
#     4. Re-processing subsets for higher resolution analyses
#
#   We build on a Seurat object that has already undergone:
#     - QC filtering
#     - Normalization
#     - Variable feature selection
#     - PCA, clustering, UMAP
#
#   Each step is heavily commented to clarify purpose and parameter choices.
###############################################################

##########################
# 0) Environment Setup
##########################



# If not installed, uncomment to install:
# install.packages("Seurat")
# install.packages("patchwork")
# install.packages("dplyr")
# install.packages("ggplot2")

library(Seurat)
library(patchwork)
library(dplyr)
library(ggplot2)

##########################
# 1) Data Loading
##########################

sc_object <- readRDS("sc_object_postQC.rds")

# For demonstration, we know that the object is already:
#  - Quality controlled
#  - Normalized (LogNormalize)
#  - Identified variable features
#  - Scaled
#  - PCA done
#  - Clusters found
#  - UMAP or t-SNE performed
#


##########################
# 2) Basic Cell Type Annotation
##########################


# If you're coing from the first module, you will still have the below cluster_markers
# variable saved in your environment, which means that you don't have to remake it
# cluster_markers <- FindAllMarkers(sc_object, 
#                                   only.pos = TRUE, ## change this to false if you want to include downregulated genes
#                                   min.pct = 0.25, 
#                                   logfc.threshold = 0.25)

cluster_markers %>% 
  group_by(cluster) %>% 
  slice_max(n = 5, order_by = avg_log2FC) %>% 
  print(n=45)


# We'll use canonical marker genes to identify cell types. For human PBMCs, 
# typical marker genes might include:
# - IL7R, CCR7 (Naive CD4+ T)
# - CD14, LYZ (CD14+ Mono)
# - IL7R, S100A4  (Memory CD4+)
# - MS4A1 (B cells)
# - CD8A (CD8 T cells)
# - FCGR3A, MS4A7 (FCGR3A+ Mono)
# - GNLY / NKG7 (NK cells)
# - FCER1A, CST3 (DC)
# - PPBP (Platelets)
# 
# We can start by checking the expression of known markers in each cluster.

FeaturePlot(sc_object, features = c("IL7R", "CD14", "S100A4", "MS4A1", "CD8A", "FCGR3A", "GNLY", "FCER1A", "PPBP"), 
            reduction = "umap", pt.size = 0.5, ncol = 3)

# Optionally, look at a dot plot which shows average expression across clusters:
DotPlot(sc_object, 
        features = c("IL7R", "CD14", "S100A4", "MS4A1", "CD8A", "FCGR3A", "GNLY", "FCER1A", "PPBP"), 
        cols = c("lightgrey", "blue")) + 
  RotatedAxis()



## I'll show two methods of annotation

###################################
# Method 1 
###################################
# From visual inspection and knowledge of typical PBMC markers, 
# we can assign cluster identities. Let's see how many clusters we have:
levels(Idents(sc_object))

# We'll create a vector named new.cluster.ids to relabel them.
# Note: The actual cluster IDs for pbmc3k might differ from run to run 
# due to random seed or parameters, so adapt as needed.

# In a fresh object processing, clusters are labeled numerically in the 
# 'seurat_clusters' column 
head(sc_object@meta.data$seurat_clusters)
# and these are what the current ident is set as
head(Idents(sc_object))

# If for whatever reason you want to store current cluster IDs, 
# you can do so with:
current.cluster.ids <- levels(sc_object)

# Suppose we've identified the clusters with the following broad labels
# (this mapping is approximate; confirm with your own analysis):
new.cluster.ids <- c("Naive CD4 T",     # Cluster 0
                     "CD14+ Mono",      # Cluster 1
                     "Memory CD4+",     # Cluster 2
                     "B",               # Cluster 3
                     "CD8+ T",          # Cluster 4
                     "FCGR3A+ Mono",    # Cluster 5
                     "NK",              # Cluster 6
                     "DC",              # Cluster 7
                     "Platelet")        # Cluster 8

new.cluster.ids

# Make sure the length of new.cluster.ids matches the number of current clusters:
length(new.cluster.ids) 
# [1] 9
length(current.cluster.ids)
# [1] 9

# Then rename the identities:
names(new.cluster.ids) <- current.cluster.ids
sc_object <- RenameIdents(sc_object, new.cluster.ids)
levels(Idents(sc_object))

# Lets visualize the data with the annotated clusters
DimPlot(sc_object, reduction = "umap", label = TRUE, pt.size = 0.5) + 
  ggtitle("UMAP with annotated cell types")


######### Method 2 ###############
# Here I'll eset the idents back to seurat_clsuters for demonstrative purposes
Idents(sc_object) <- "seurat_clusters"

sc_object <- RenameIdents(sc_object, '0' = "Naive CD4 T",'1' = "CD14+ Mono",
                          '2' = "Memory CD4+", '3' = 'B', '4' = 'CD8+ T', '5' = 'FCGR3A+ Mono'
                          , '6' = 'NK', '7' = 'DC', '8' = 'Platelets')

# Let's visualize again with annotated clusters:
DimPlot(sc_object, reduction = "umap", label = TRUE, pt.size = 0.5) + 
  ggtitle("UMAP with annotated cell types")



# After annotating our clusters we can store these identities 
# in the metadata if we want a new column that doesn't overwrite
# the "seurat_clusters" column:
sc_object[["celltype"]] <- Idents(sc_object)

# Now you can set the ident to your metadata column 
Idents(sc_object) <- "celltype"
head(sc_object@meta.data)
# or View(sc_object@meta.data)


##########################
# 3) Differential Expression (DE)
##########################

# A) Differential expression bewteen two clusters (e.g., identify markers of "NK cells" vs. B Cells)
# This would be useful for comparing two very similar or funcitonally related
# cell types. Eg. we might want to know how NK cells differ from B cells, rather
# that from a mix of other cell types. 
# We can set the identity to the cell types we just assigned if you haven't already:
Idents(sc_object) <- "celltype"

# Compare differential expression between NK cells and B cells.
# This analysis aims to identify genes that are differentially expressed in NK cells (group of interest)
# compared to B cells (reference group). Adjusting parameters like min.pct and logfc.threshold helps focus on
# genes that are both robustly and significantly different between the two groups.
markers_NK_vs_B <- FindMarkers(
  object = sc_object,
  ident.1 = "NK",          # Group of interest: NK cells
  ident.2 = "B",           # Reference group: B cells
  min.pct = 0.25,          # Only test genes expressed in at least 25% of cells in either group
  logfc.threshold = 0.25   # Consider genes with a minimum log2 fold change of 0.25
)



# Inspect the top markers based on adjusted p-value (most statistically significant differences)
top_markers_pval <- markers_NK_vs_B[order(markers_NK_vs_B$p_val_adj), ]
head(top_markers_pval, 10)

# Inspect the top markers with the highest average log2 fold change (most upregulated in NK cells)
top_markers_up <- markers_NK_vs_B[order(-markers_NK_vs_B$avg_log2FC), ]
head(top_markers_up, 10)

# Inspect the markers with the lowest average log2 fold change (most downregulated in NK cells compared to B cells)
top_markers_down <- markers_NK_vs_B[order(markers_NK_vs_B$avg_log2FC), ]
head(top_markers_down, 10)


# B) Or identify markers for each cell type (like "FindAllMarkers" by celltype).
all_celltype_markers <- FindAllMarkers(
  sc_object,
  only.pos = TRUE, 
  min.pct = 0.25, 
  logfc.threshold = 0.25
)

# Check top 5 markers per cell type:
all_celltype_markers %>%
  group_by(cluster) %>%
  slice_max(n = 5, order_by = avg_log2FC) %>%
  print(n = 45)

# C) If you have experimental groups (e.g., treatment vs. control), you can 
# switch the object identity to that group’s metadata column and run FindMarkers.
# Note: 

# Since PBMC 3k doesn't have a built-in "condition", let's simulate a random one:
sc_object$condition <- sample(c("Control", "Treatment"), size = ncol(sc_object), replace = TRUE)

# Now we can do:
DimPlot(sc_object, reduction = "umap", split.by = "condition", ncol = 2) +
  ggtitle("UMAP faceted by condition")

Idents(sc_object) <- "condition"

markers_treated_vs_control <- FindMarkers(sc_object, ident.1 = "Treatment", ident.2 = "Control")

head(markers_treated_vs_control[order(markers_treated_vs_control$p_val_adj), ], 10)
head(markers_treated_vs_control[order(markers_treated_vs_control$avg_log2FC), ], 10)
tail(markers_treated_vs_control[order(markers_treated_vs_control$avg_log2FC), ], 10)
# These are all about as insignificant as you can get because the 'conditions' are simulated
# and pretty well exactly the same

# If you change your ident (like we have here to 'condition'), make sure to change
# it back if you are continuing on with your analysis
Idents(sc_object) <- "celltype"

##########################
# 4) Subsetting a Cell Type (e.g., T cells) and Re-processing
##########################

levels(sc_object)

# Often, after broad clustering/annotation, you want to focus on one cell type 
# (like T cells) and sub-cluster to see more nuanced populations.

# Let's subset T cells (CD4 and CD8). We can do this by selecting the relevant identity:
# For example, we suspect "Naive CD4 T" and "CD8+ T" are in the object; let's subset them:

tcell_subset <- subset(sc_object, idents = c("Naive CD4 T", "CD8+ T"))

# Confirm the dimension (number of cells/genes) in the subset
dim(tcell_subset)

# Now we want to re-run the typical workflow for just these cells.
# We'll do a short pipeline (Normalize -> FindVariableFeatures -> Scale -> PCA -> Clustering -> UMAP).
# If we had SCTransform from the beginning, we’d also do SCTransform here, 
# but let's keep consistent with a LogNormalize approach.

# Re-normalize (LogNormalize):
tcell_subset <- NormalizeData(tcell_subset, 
                              normalization.method = "LogNormalize", 
                              scale.factor = 10000)

# Find variable features:
tcell_subset <- FindVariableFeatures(tcell_subset, selection.method = "vst", nfeatures = 2000)

# Scale data:
# Optionally, we can regress out percent.mt or cell cycle scores if we want:
all.genes <- rownames(tcell_subset)
tcell_subset <- ScaleData(tcell_subset, features = all.genes)

# PCA:
tcell_subset <- RunPCA(tcell_subset, features = VariableFeatures(tcell_subset))
ElbowPlot(tcell_subset)

# Clustering:
tcell_subset <- FindNeighbors(tcell_subset, dims = 1:10)
tcell_subset <- FindClusters(tcell_subset, resolution = 0.4) 
tcell_subset <- RunUMAP(tcell_subset, dims = 1:10)

# Visualize:
DimPlot(tcell_subset, reduction = "umap", label = TRUE, pt.size = 0.5) +
  ggtitle("T cell Subset UMAP")

##########################
# 5) Refining Annotation within the Subset
##########################

# Now we might identify subtypes (e.g., regulatory T cells, naive vs. memory T cells, etc.)
Idents(tcell_subset) <- "seurat_clusters"

subset_markers <- FindAllMarkers(
  tcell_subset,
  only.pos = TRUE, 
  min.pct = 0.25, 
  logfc.threshold = 0.25
)

# Check top 5 markers per cell type:
subset_markers %>%
  group_by(cluster) %>%
  slice_max(n = 5, order_by = avg_log2FC) %>%
  print(n = 45)


# We can do so by looking at known markers:


p1 <- FeaturePlot(tcell_subset, 
                  features = c("FGFBP2", "GZMH", "NCR3", "CCR7"), 
                  reduction = "umap")
p2 <- DimPlot(tcell_subset, reduction = "umap", label = TRUE, pt.size = 0.5)
p1 / p2 # patchwork

# C0 CCR7
# C1 NCR3
# C2 GZMH

# If we find distinct expression patterns, we can label sub-clusters accordingly. 
# E.g., let's see how many sub-clusters we have:
levels(Idents(tcell_subset))

# Suppose cluster 0 is naive T cells, cluster 1 is memory T cells, cluster 2 is regulatory T, etc.
# We'll rename as an example:
current.t.subcluster.ids <- levels(Idents(tcell_subset))


DimPlot(tcell_subset, reduction = "umap", label = TRUE, pt.size = 0.5) +
  ggtitle("Detailed T cell Subset Annotation")

new.t.subcluster.ids <- c("CCR7+ cells", # cluster 0
                          "NCR3+ cells", # cluster 1
                          "GZMH+ cells")  # cluster 2

# Make sure lengths match:
length(current.t.subcluster.ids)
length(new.t.subcluster.ids)

names(new.t.subcluster.ids) <- current.t.subcluster.ids
tcell_subset <- RenameIdents(tcell_subset, new.t.subcluster.ids)

DimPlot(tcell_subset, reduction = "umap", label = TRUE, pt.size = 0.5) +
  ggtitle("Detailed T cell Subset Annotation")



# You'll notice that our condition element was preserved
DimPlot(tcell_subset, reduction = "umap", label = TRUE, pt.size = 0.5, split.by = 'condition') 


# We could further do differential expression among these T sub-clusters or 
# use advanced T cell annotation references.

##########################
# 6) Saving and Reintegrating Subset Annotations
##########################

# A) Save the T cell subset as a separate object:
saveRDS(tcell_subset, file = "tcell_subset_seurat.rds")

# B) If you want to integrate your new T cell annotations back into the original object,
# you could add a metadata column. For example, for cells in the T cell subset, we store 
# the subcluster identity. For all other cells, we keep their original annotation.

# First, create a data frame with cell barcodes and new identities:
tcell_ann_df <- data.frame(
  cell_barcode = colnames(tcell_subset),
  tcell_subcluster = Idents(tcell_subset)
)

# Next, we can add a new metadata column to the original object:
# We'll initialize a new column (e.g. "TcellSubtype") in sc_object as "NA"
sc_object$TcellSubtype <- NA

# For cells in the T cell subset, fill in the subtype:
common_barcodes <- intersect(tcell_ann_df$cell_barcode, colnames(sc_object))
sc_object$TcellSubtype[common_barcodes] <- tcell_ann_df$tcell_subcluster[match(common_barcodes, tcell_ann_df$cell_barcode)]

# Now sc_object has a column "TcellSubtype" that includes sub-cluster annotation for T cells 
# and NA for non-T cells.

##########################
# 7) Optional: Overlapping / Combined Visualizations
##########################

# You can highlight your new T cell subtype metadata on the original object's UMAP:
Idents(sc_object) <- "TcellSubtype"
DimPlot(sc_object, reduction = "umap", label = TRUE, pt.size = 0.5) +
  ggtitle("Original UMAP Colored by T Cell Subtypes (Others = NA)")

# Or you might want to keep the broad celltype annotation:
Idents(sc_object) <- "celltype"
DimPlot(sc_object, reduction = "umap", label = TRUE, pt.size = 0.5) +
  ggtitle("Original UMAP with Broad Cell Type Annotation")

##########################
# 8) Additional Notes / Troubleshooting
##########################

# - Always verify that your marker-based annotations align with biological expectations 
#   (e.g., known marker genes, public references, or SingleR comparisons).
# - For differential expression, watch out for confounding factors like 
#   library size, batch effects, etc.
# - When subsetting, remember to re-run normalization steps if you plan to recluster. 
#   The new library sizes and distributions can differ from the full dataset.
# - If sub-clustering yields many or few clusters, adjust the resolution parameter.

##########################
# 9) Session Info
##########################

sessionInfo()

###############################################################
# End of Script
# 
# Key Takeaways:
# 1. Once your data is clustered, use canonical markers or known references 
#    to assign biological identities to clusters.
# 2. Differential expression can be done cluster vs cluster or 
#    condition vs condition.
# 3. Subsetting allows deeper resolution on a particular cell type. 
#    Re-run the workflow (normalization, HVG, scaling, PCA, clustering) 
#    on the subset for nuanced populations.
# 4. Always document your annotation logic and cross-check with known databases 
#    or references for final labeling.
###############################################################