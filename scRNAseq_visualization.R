###############################################################
# Title: Advanced Visualization & Analysis Techniques in scRNA-seq
# Author: [Your Name]
# Date: [Workshop Date]
#
# Description:
#   A comprehensive script demonstrating the most valuable data visualization
#   and analysis methods in scRNA-seq, explaining *when* and *why* each is used.
###############################################################

##########################
# 0) Environment Setup
##########################

# If not installed, uncomment to install:
# install.packages("Seurat")
# install.packages("patchwork")
# install.packages("ggplot2")
# install.packages("dplyr")

library(Seurat)
library(patchwork)
library(ggplot2)
library(dplyr)

##########################
# 1) Load/Prepare Data
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

###############################################################
# Section A: Whole-transcriptome & Cluster-level Visualization
###############################################################

##########################
# 2) DimPlot / Cluster Map
##########################

# WHEN to use:
#   - To visualize overall cell clustering in a reduced 2D or 3D space (e.g., UMAP, t-SNE, PCA).
#   - To quickly see the major clusters or group identities (e.g., cell types, conditions).
#
# WHY:
#   - Provides an intuitive overview of high-dimensional data.
#   - Essential first step in scRNA-seq data presentation.

DimPlot(sc_object, reduction = "umap", group.by = "seurat_clusters", label = TRUE) +
  ggtitle("UMAP of PBMC 3k (Clusters)")

# Variation: color by another metadata column (e.g., cell type if annotated)
#   DimPlot(sc_object, reduction = "umap", group.by = "celltype", label = TRUE)

##########################
# 3) FeaturePlot
##########################

# WHEN to use:
#   - To visualize expression of one or more genes across all cells on a 2D embedding (UMAP, t-SNE).
#
# WHY:
#   - Quickly see how well a gene (or genes) segregate the data.
#   - Useful for verifying marker gene expression within clusters.

# Can plot one or multiple genes

# Plotting one gene
FeaturePlot(sc_object, features = "CD3D", reduction = "umap", 
            pt.size = 0.5) 


FeaturePlot(sc_object, features = c("MS4A1", "CD3D"), reduction = "umap", 
            pt.size = 0.5) 

# Variation:
#   - Add columns if you want to customize how the multi-gene plot is laid out with 'ncol = X'

FeaturePlot(sc_object, features = c("MS4A1", "CD3D", "LTB"), reduction = "umap", 
            pt.size = 0.5) ### automatically assigns placement
#vs 
FeaturePlot(sc_object, features = c("MS4A1", "CD3D", "LTB"), reduction = "umap", 
            pt.size = 0.5, ncol = 3) 


#   - Adjust point size or color scale as needed.

##########################
# 4) DotPlot
##########################

# WHEN to use:
#   - To compare expression of multiple genes across multiple clusters (or conditions).
#   - Particularly powerful for showcasing marker panels or pathways across clusters.
#
# WHY:
#   - Summarizes both the percentage of cells expressing a gene and the average expression level.
#   - Great for presenting many genes in a concise figure.

DotPlot(sc_object, 
        features = c("MS4A1", "CD3D", "CD14", "GNLY", "FCGR3A", "CD8A"), 
        cols = c("lightgrey", "blue")) +
  RotatedAxis() +
  ggtitle("DotPlot of key markers across clusters")

# Variation:
#   - Customize color scales, dot sizes, grouping variables.

##########################
# 5) Violin Plot
##########################

# WHEN to use:
#   - To visualize the distribution of a gene's expression across clusters or conditions.
#
# WHY:
#   - Shows the full distribution (density) of expression values. 
#   - Useful for demonstrating variability within each cluster.

VlnPlot(sc_object, features = c("CD14"), pt.size = 0.1) +
  ggtitle("Violin Plot of CD14 expression")

# Variation:
#   - Plot multiple genes side-by-side: `VlnPlot(sc_object, features = c("CD14", "MS4A1"))`
#   - Adjust or remove `pt.size` for clarity (large datasets can produce over-plotted points).


###############################################################
# Section B: Marker Gene & Cluster Relationship Visualizations
###############################################################

##########################
# 6) Heatmap of Marker Genes
##########################

# WHEN to use:
#   - After finding markers (e.g., with `FindAllMarkers`), you can visualize top genes per cluster.
#   - Great for summarizing cluster-specific signatures.
#
# WHY:
#   - Heatmaps provide a clear overview of how each cluster expresses selected genes.
#   - Ideal for cluster annotation or comparing expression patterns across multiple markers.

# Let's identify markers for demonstration. We limit for speed:
cluster_markers <- FindAllMarkers(sc_object, only.pos = TRUE, 
                                  min.pct = 0.25, logfc.threshold = 0.25)
top_markers <- cluster_markers %>% group_by(cluster) %>% slice_max(n = 5, order_by = avg_log2FC)

# We can use DoHeatmap to visualize expression of top markers across clusters:
DoHeatmap(sc_object, features = top_markers$gene) +
  scale_fill_gradientn(colors = c("blue", "white", "red")) + # here I've just customized the colours
  ggtitle("Heatmap of Top Cluster Markers")

# Variation:
#   - Add clustering or annotation to rows or columns.
#   - For larger gene sets, consider external packages like ComplexHeatmap for more customization.


###############################################################
# Section C: Condition-level & Multi-group Visualizations
###############################################################

##########################
# 7) Split or Facet Plots
##########################

# WHEN to use:
#   - If you have multiple conditions (e.g., Control vs. Treatment), 
#     or time points (Day0, Day7, Day14), you can compare them side by side.
#
# WHY:
#   - Provides a direct visual comparison of cluster structure or expression across conditions.

# Example: If sc_object has a "condition" metadata column, we can do:
# DimPlot(sc_object, reduction = "umap", split.by = "condition", ncol = 2)

# Since PBMC 3k doesn't have a built-in "condition", let's simulate a random one:

sc_object$condition <- sample(c("Control", "Treatment"), size = ncol(sc_object), replace = TRUE)

# Now we can do:
DimPlot(sc_object, reduction = "umap", split.by = "condition", ncol = 2) +
  ggtitle("UMAP faceted by condition")

#or

DimPlot(sc_object, reduction = "umap", group.by = "condition") +
  ggtitle("UMAP faceted by condition")

##########################
# 8) Combined Patchwork Figures
##########################

# WHEN to use:
#   - To present multiple plots together in a single figure panel.
#
# WHY:
#   - Tells a cohesive story in a publication or presentation 
#   - Minimizes the overhead of multiple figure windows.

p1 <- DimPlot(sc_object, reduction = "umap", label = TRUE) + ggtitle("Clusters")
p2 <- FeaturePlot(sc_object, features = "MS4A1") + ggtitle("MS4A1 expression")

# Combine plots side-by-side (default horizontal layout)
combined_horizontal <- p1 + p2
combined_horizontal

# Combine plots vertically
combined_vertical <- p1 / p2
combined_vertical

# Combine three plots with two columns
p3 <- VlnPlot(sc_object, features = "MS4A1") + ggtitle("MS4A1 Violin Plot")

combined_layout <- p1 + p2 + p3 + plot_layout(ncol = 2)
combined_layout

combined_annotated <- p1 + p2 +
  plot_annotation(
    title = "Combined UMAP and Feature Plot",
    subtitle = "Visualizing Clusters and MS4A1 Expression"
  ) &
  theme(
    plot.title = element_text(hjust = 0.5),      # Center the title
    plot.subtitle = element_text(hjust = 0.5)      # Center the subtitle
  )

combined_annotated


##########################
# 9) Example using ggplot for data representation
##########################


## Here we're going to make a stacked bar plot to represent
# cellt ype proportions by condition (remember we invented the condition)
library(scales)

Idents(sc_object) <- "celltype"
# Extract metadata from your Seurat object
metadata <- sc_object@meta.data

# Summarize counts per condition and cell type
cell_data <- metadata %>%
  group_by(condition, celltype) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(condition) %>%
  mutate(proportion = count / sum(count)) %>%
  ungroup()

# Create a stacked bar plot of cell type proportions by condition
stacked_bar <- ggplot(cell_data, aes(x = condition, y = proportion, fill = celltype)) +
  geom_bar(stat = "identity") +
  scale_y_continuous(labels = percent_format()) +
  labs(
    title = "Cell Type Proportions by Condition",
    x = "Condition",
    y = "Proportion",
    fill = "Cell Type"
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    legend.position = "right"
  )

print(stacked_bar)



#### We could also use ggplot to show data in other ways

# Just inventing some data here 
df <- data.frame(
  condition = rep(c("Control", "Treatment"), each = 100),
  gene_expression = c(rnorm(100, mean = 5), rnorm(100, mean = 7)),
  cell_type = rep(c("A", "B"), 100)
)

# Basic scatter plot with advanced customizations:
p <- ggplot(df, aes(x = condition, y = gene_expression, color = cell_type)) +
  geom_jitter(width = 0.2, size = 2, alpha = 0.7) +
  labs(
    title = "Gene Expression Across Conditions",
    subtitle = "Example of advanced ggplot customization",
    x = "Condition",
    y = "Expression Level",
    caption = "Data simulated for demonstration purposes",
    color = "Cell Type"
  ) +
  scale_color_brewer(palette = "Set1") +            # Use a brewer palette for colors
  theme_minimal(base_size = 14) +                    # Clean minimal theme
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "top",
    panel.grid.major = element_line(color = "grey90"),
    panel.grid.minor = element_blank()
  ) +
  facet_wrap(~ cell_type)                           # Facet by cell type for comparison

# Display the plot
print(p)



##########################
# 10) Color Palettes & Aesthetics
##########################

# WHEN to use:
#   - Whenever you want a publication-ready figure or improved clarity for color-blind audiences.
#
# WHY:
#   - Good color choices emphasize data patterns and increase figure readability.

# Example: using RColorBrewer:
library(RColorBrewer)
display.brewer.all()


DimPlot(sc_object, reduction = "umap", label = TRUE) + 
  scale_color_brewer(palette = "Set2") +
  ggtitle("UMAP with Set2 Palette")

DimPlot(sc_object, reduction = "umap", label = TRUE) + 
  scale_color_brewer(palette = "Dark2") +
  ggtitle("UMAP with Set2 Palette")



# Can also use hexidecimal codes or just type what you want
# Caveat: this method you have to specify the exact number of colours to your clusters

DimPlot(tcell_subset, reduction = "umap", label = TRUE) +
  scale_color_manual(values = c("#984ea3", "pink", "#4daf4a")) +
  ggtitle("UMAP with Custom Color Vector")


### You can also preset colours by ident
# Remember, you can view your idents by:
levels(sc_object)

# Next create your named vector
my_colors <- c(
  "Naive CD4 T"   = "#e41a1c",  # red
  "CD14+ Mono"    = "#377eb8",  # blue
  "Memory CD4+"   = "#4daf4a",  # green
  "B"             = "#984ea3",  # purple
  "CD8+ T"        = "#ff7f00",  # orange
  "FCGR3A+ Mono"  = "#ffff33",  # yellow
  "NK"            = "#a65628",  # brown
  "DC"            = "#f781bf",  # pink
  "Platelets"     = "#999999"   # gray
)

DimPlot(sc_object, reduction = "umap", label = TRUE) +
  scale_color_manual(values = my_colors) +
  ggtitle("UMAP with Custom Manual Colors")


###############################################################
# Section E: Practical Tips & Saving Outputs
###############################################################

##########################
# 11) Downloading High Quality Plots for Publications
##########################

# Example plot creation (replace 'my_plot' with your ggplot object, e.g., combined_annotated)
my_plot <- DimPlot(sc_object, reduction = "umap", label = TRUE) +
  ggtitle("UMAP Clusters")

# Export as a PNG for high-resolution raster output (300 dpi)
ggsave("UMAP_clusters.png", plot = my_plot, width = 6, height = 5, dpi = 300)

# Export as a PDF (vector-based format)
ggsave("UMAP_clusters.pdf", plot = my_plot, width = 6, height = 5)


# Tips:
# - Use vector formats (PDF, SVG) for large figures to ensure they remain sharp when zoomed.
# - Ensure consistent labeling, legends, and color themes across all figures.

##########################
# 12) Summary & Next Steps
##########################

# - DimPlot & FeaturePlot are your go-to for cluster overviews and gene expression mapping.
# - DotPlot, Violin, and RidgePlots highlight distributions and group differences
# - Heatmaps help with multi-gene patterns and cluster marker validation.
# - Facet/Split or patchwork: essential for multi-condition or multi-plot narratives.
# - Always tailor color palettes and theming for clarity and accessibility.




