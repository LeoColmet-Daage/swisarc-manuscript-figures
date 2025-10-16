#!/usr/bin/env Rscript

# -----------------------------------------------------------------------------
# Metadata
# -----------------------------------------------------------------------------
# Project: SWISARC manuscript figures
# Figure: 1H - UMAPs
# Author: Leo Colmet Daage
# Purpose: Build UMAPs
# Date: 2024-02-07
#
# Workflow:
# 1. Load pre-computed UMAP coordinates from Seurat object
# 2. Extract and rename sample identifiers for publication-ready labels
# 3. Define tumor cluster annotations based on cluster identities
# 4. Create color palettes for sample and tumor cluster visualization
# 5. Generate UMAP plots colored by sample ID and tumor cluster
# 6. Generate UMAP plot colored by tumor cluster with meaningful labels
# 7. Extract legends from UMAP plots
# 8. Combine legends and UMAP plots
# 9. Export figure as PDF and PNG
#
# Requirements:
# - Pre-computed Seurat object with UMAP embeddings
# - R packages: ggplot2, Seurat, ggpubr, here, cowplot


# -----------------------------------------------------------------------------
# Package loading: Load required libraries with error handling
# -----------------------------------------------------------------------------
# Note: Ensure all packages are installed before running this script
# The "ggplot2", "Seurat", "ggpubr", "here" and "cowplot" packages are available on CRAN:
# install.packages(c("ggplot2", "Seurat", "ggpubr", "here", "cowplot"))

required_packages <- c("ggplot2", "Seurat", "ggpubr", "here", "cowplot")

# Check if packages are available and load them
for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(paste("Package", pkg, "is required but not installed. Please install it first."))
  }
  suppressPackageStartupMessages(library(pkg, character.only = TRUE))
}

# -----------------------------------------------------------------------------
# Configuration: Key parameters for analysis
# -----------------------------------------------------------------------------
# Analysis parameters
POINT_SIZE <- 0.2     # Size of points in UMAP plots
LEGEND_POINT_SIZE <- 6  # Size of points in legend

# File paths (relative to project root)
UMAP_FILE <- here::here("UMAP", "merge_tumor_harmony_SCTransform_pca_40_0.2.rda")

# -----------------------------------------------------------------------------
# Reproducibility
# -----------------------------------------------------------------------------
set.seed(37)
# Ensure the working directory is set to the project root
setwd(here::here())

# -----------------------------------------------------------------------------
# Load pre-computed UMAP data from Seurat object
# -----------------------------------------------------------------------------
# Load the Seurat object containing UMAP embeddings and metadata
load(UMAP_FILE)

# -----------------------------------------------------------------------------
# Extract and rename sample identifiers for publication-ready labels
# -----------------------------------------------------------------------------
# Extract original sample identifiers from Seurat object metadata
SAMPLE_ID <- as.factor(sobj@meta.data$orig.ident)

# Rename sample identifiers to publication-ready format (P##)
# Maps internal sample names to patient identifiers for manuscript
sample_rename_map <- c(
  "BLI_PIE_GE_10112021_GE" = "P36",
  "STE_THI_GE_16022021_GE" = "P16",
  "SAN_MAR_GE_22072021_GE" = "P34",
  "FER_MAR_GE_08062021_GE" = "P23",
  "MAR_FRA_GE_21102021_GE" = "P6B",
  "RAG_JEA_GE_17032021_GE" = "P35",
  "SMARCB1-GR2_GE"         = "P12",
  "SARC_EPI_GR1_GE"        = "P6A"
)
levels(SAMPLE_ID) <- ifelse(levels(SAMPLE_ID) %in% names(sample_rename_map),
                           sample_rename_map[levels(SAMPLE_ID)],
                           levels(SAMPLE_ID))

# Add renamed sample IDs to Seurat object metadata
sobj@meta.data$SAMPLE_ID <- SAMPLE_ID

# -----------------------------------------------------------------------------
# Define color palettes for sample and tumor cluster visualization
# -----------------------------------------------------------------------------
# Define distinct colors for each sample
SAMPLE_ID_COLORS <- c(
  "P36" = "#ff71a8",
  "P16" = "#9d2425",
  "P34" = "#4da3ff",
  "P6A" = "#006739",
  "P6B" = "#6add83",
  "P23" = "#967200",
  "P35" = "#8826ae",
  "P12" = "#ed7d01"
)
 
# Set factor levels for consistent sample ordering in plots
sobj@meta.data$SAMPLE_ID <- factor(sobj@meta.data$SAMPLE_ID, levels=c("P36","P16","P34","P6A","P6B","P23","P35","P12"))

# Define color palette for tumor cluster visualization
TUMOR_CLUSTER_COLORS <- c(
  "0" = "#8dd3c7",
  "1" = "#ffd700",
  "2" = "#bebada",
  "3" = "#80b1d3",
  "4" = "#fdb462",
  "5" = "#fb8072",
  "6" = "#b3de69",
  "7" = "#fccde5",
  "8" = "#d9d9d9",
  "9" = "#bc80bd"
)

# Define color palettes for tumor cluster labeled legend
TUMOR_LABELED_COLORS <- c(
  "Cell migration_0"               = "#8dd3c7",
  "Apoptotic signaling_1"          = "#ffd700",
  "Regulation of growth_2"         = "#bebada",
  "Immune response_3"              = "#80b1d3",
  "Ribosomal_IFN response_4"       = "#fdb462",
  "EMT_5"                          = "#fb8072",
  "Cycling cells_6"                = "#b3de69",
  "LncRNA_7"                       = "#fccde5",
  "Circulatory system development_8"= "#d9d9d9",
  "Response to metal ion_9"        = "#bc80bd"
)

# -----------------------------------------------------------------------------
# Define tumor cluster annotations based on cluster identities
# -----------------------------------------------------------------------------
# Extract cluster identities from Seurat object active identity
TUMOR_CLUSTER <- sobj@active.ident

# Map cluster numbers to biologically meaningful tumor cluster names
# Based on marker gene expression and known tumor cluster characteristics
tumor_cluster_map <- c(
  "0" = "Cell migration_0",
  "1" = "Apoptotic signaling_1",
  "2" = "Regulation of growth_2",
  "3" = "Immune response_3",
  "4" = "Ribosomal_IFN response_4",
  "5" = "EMT_5",
  "6" = "Cycling cells_6",
  "7" = "LncRNA_7",
  "8" = "Circulatory system development_8",
  "9" = "Response to metal ion_9"
)
levels(TUMOR_CLUSTER) <- tumor_cluster_map[levels(TUMOR_CLUSTER)]

# Add tumor cluster annotations to Seurat object metadata
sobj@meta.data$TUMOR_CLUSTER <- TUMOR_CLUSTER

# -----------------------------------------------------------------------------
# Generate UMAP plots colored by sample ID and tumor cluster
# -----------------------------------------------------------------------------
# Create UMAP plot colored by sample ID
umap_sample_fig1h <- DimPlot(
  object = sobj, 
  group.by = "SAMPLE_ID", 
  pt.size = POINT_SIZE, 
  cols = SAMPLE_ID_COLORS
) + 
  ggtitle("Samples") + 
  theme(
    axis.title = element_text(size = 20), 
    legend.title = element_text(size = 20, face = "bold"), 
    legend.text = element_text(size = 12), 
    legend.position = "top"
  ) +
  guides(
    colour = guide_legend(
      override.aes = list(size = LEGEND_POINT_SIZE), 
      nrow = 2
    )
  ) +
  xlab("UMAP_1") +
  ylab("UMAP_2")

# Create UMAP plot colored by tumor cluster with numeric labels and no legend
umap_tumor_cluster_fig1h <- DimPlot(
  object = sobj, 
  group.by = "ident", 
  pt.size = POINT_SIZE, 
  cols = TUMOR_CLUSTER_COLORS
) + 
  ggtitle("") + 
  theme(
    axis.title = element_text(size = 20)
  ) + 
  xlab("UMAP_1") + 
  ylab("") + 
  NoLegend() 

# Add label to plot
umap_tumor_cluster_fig1h <- LabelClusters(
  plot = umap_tumor_cluster_fig1h,
  id = "ident",
  box = FALSE,
  size = 8
)
 
# Create UMAP plot colored by tumor cluster with meaningful labels
umap_tumor_cluster_fig1h_labeled <- DimPlot(
  object = sobj,
  group.by = "TUMOR_CLUSTER",
  pt.size = POINT_SIZE,
  cols = TUMOR_LABELED_COLORS
) +
  ggtitle("") +
  theme(
    axis.title = element_text(size = 20),
    legend.title = element_text(size = 16, face = "bold"),
    legend.text = element_text(size = 12),
    legend.position = "top"
  ) +
  ylab("") +
  guides(
    colour = guide_legend(
      title.position = "top",
      override.aes = list(size = LEGEND_POINT_SIZE),
      nrow = 5
    )
  ) +
  xlab("UMAP_1") +
  labs(color = 'Malignant cell subclusters')

# Extract legend from UMAP plot
legend_sample_fig1h <- as_ggplot(get_legend(umap_sample_fig1h))
legend_tumor_cluster_fig1h <- as_ggplot(get_legend(umap_tumor_cluster_fig1h_labeled))

# Combine legends side by side
legends<-ggarrange(legend_sample_fig1h,legend_tumor_cluster_fig1h,ncol = 2)

# Combine both UMAP plots side by side for publication figure
umap_fig1h <- ggarrange(umap_sample_fig1h + NoLegend(), umap_tumor_cluster_fig1h, ncol = 2)

# -----------------------------------------------------------------------------
# Export figure to PDF and PNG (Figure 1H)
# -----------------------------------------------------------------------------
# Save as PDF
pdf(file.path(here::here("Figures", "Figure_1H_umaps.pdf")), width = 14, height = 7)
print(plot_grid(legends,umap_fig1h,ncol = 1,nrow = 2,rel_heights = c(0.4, 0.6)))
dev.off()

# Save as PNG
png(file.path(here::here("Figures", "Figure_1H_umaps.png")), width = 4200, height = 2100, res = 300)
print(plot_grid(legends,umap_fig1h,ncol = 1,nrow = 2,rel_heights = c(0.4, 0.6)))
dev.off()