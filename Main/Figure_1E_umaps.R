#!/usr/bin/env Rscript

# -----------------------------------------------------------------------------
# Metadata
# -----------------------------------------------------------------------------
# Project: SWISARC manuscript figures
# Figure: 1E - UMAPs
# Author: Leo Colmet Daage
# Purpose: Build UMAPs
# Date: 2024-02-07
#
# Workflow:
# 1. Load pre-computed UMAP coordinates from Seurat object
# 2. Extract and rename sample identifiers for publication-ready labels
# 3. Define cell type annotations based on cluster identities
# 4. Create color palettes for sample and cell type visualization
# 5. Generate UMAP plots colored by sample ID and cell type
# 6. Combine plots and export as PDF and PNG for publication
#
# Requirements:
# - Pre-computed Seurat object with UMAP embeddings
# - R packages: ggplot2, Seurat, hues, ggpubr, here


# -----------------------------------------------------------------------------
# Package loading: Load required libraries with error handling
# -----------------------------------------------------------------------------
# Note: Ensure all packages are installed before running this script
# The "ggplot2", "Seurat", "hues", "ggpubr" and "here" packages are available on CRAN:
# install.packages(c("ggplot2", "Seurat", "hues", "ggpubr", "here"))

required_packages <- c("ggplot2", "Seurat", "hues", "ggpubr", "here")

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
UMAP_FILE <- here::here("UMAP", "SWISARC_all_samples_nokeep_SCTransform_pca_33_0.5_umap.rda")

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
# Define color palettes for sample visualization
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

# -----------------------------------------------------------------------------
# Define cell type annotations based on cluster identities
# -----------------------------------------------------------------------------
# Extract cluster identities from Seurat object active identity
CELL_TYPE <- sobj@active.ident

# Map cluster numbers to biologically meaningful cell type names
# Based on marker gene expression and known cell type characteristics
celltype_map <- c(
  "0"  = "Myeloid cells",
  "1"  = "Malignant cells",
  "2"  = "Malignant cells",
  "3"  = "CAFs",
  "4"  = "Malignant cells",
  "5"  = "Malignant cells",
  "6"  = "T/NK cells",
  "7"  = "Myeloid cells",
  "8"  = "T/NK cells",
  "9"  = "Myeloid cells",
  "10" = "Malignant cells",
  "11" = "CAFs",
  "12" = "Malignant cells",
  "13" = "Malignant cells",
  "14" = "Endothelial cells",
  "15" = "Malignant cells",
  "16" = "Myeloid cells",
  "17" = "Malignant cells",
  "18" = "Malignant cells",
  "19" = "Myeloid cells",
  "20" = "Malignant cells",
  "21" = "Myeloid cells",
  "22" = "Myeloid cells",
  "23" = "Malignant cells",
  "24" = "Malignant cells",
  "25" = "CAFs",
  "26" = "Malignant cells",
  "27" = "B cells",
  "28" = "B cells"
)
levels(CELL_TYPE) <- celltype_map[levels(CELL_TYPE)]

# Add cell type annotations to Seurat object metadata
sobj@meta.data$CELL_TYPE <- CELL_TYPE

# -----------------------------------------------------------------------------
# Generate UMAP plots colored by sample ID and cell type
# -----------------------------------------------------------------------------
# Create UMAP plot colored by sample ID
umap_sample_fig1e <- DimPlot(
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

# Create UMAP plot colored by cell type
umap_cell_type_fig1e <- DimPlot(
  object = sobj, 
  group.by = "CELL_TYPE", 
  pt.size = POINT_SIZE
) + 
  scale_color_iwanthue(
    hmin = 0, hmax = 360,
    cmin = 40, cmax = 70,
    lmin = 15, lmax = 85,
    random = FALSE
  ) + 
  scale_fill_iwanthue(
    hmin = 0, hmax = 360,
    cmin = 30, cmax = 80,
    lmin = 0, lmax = 100,
    random = FALSE
  ) + 
  ggtitle("Cell types") + 
  theme(
    axis.title = element_text(size = 20), 
    legend.position = "top", 
    axis.text.y = element_blank(), 
    legend.title = element_text(size = 20, face = "bold"), 
    legend.text = element_text(size = 12)
  ) +
  guides(
    colour = guide_legend(
      override.aes = list(size = LEGEND_POINT_SIZE),
      nrow = 2
    )
  ) +
  xlab("UMAP_1") +
  ylab("")

# Combine both UMAP plots side by side for publication figure
umap_fig1e <- ggarrange(umap_sample_fig1e, umap_cell_type_fig1e, ncol = 2)

# -----------------------------------------------------------------------------
# Export figure to PDF and PNG (Figure 1E)
# -----------------------------------------------------------------------------
# Save as PDF
pdf(file.path(here::here("Figures", "Figure_1E_umaps.pdf")), width = 14, height = 7)
print(umap_fig1e)
dev.off()

# Save as PNG
png(file.path(here::here("Figures", "Figure_1E_umaps.png")), width = 4200, height = 2100, res = 300)
print(umap_fig1e)
dev.off()