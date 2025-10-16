#!/usr/bin/env Rscript

# -----------------------------------------------------------------------------
# Metadata
# -----------------------------------------------------------------------------
# Project: SWISARC manuscript figures
# Figure: 1C - MCP counter
# Author: Leo Colmet Daage
# Purpose: Build MCP counter immune cell deconvolution heatmap
# Date: 2024-05-06
#
# Workflow:
# 1. Load clinical metadata and TPM expression matrices from multiple batches
# 2. Filter samples (exclude poorly differentiated chordoma and EERT)
# 3. Merge TPM matrices and apply MCPcounter immune cell deconvolution
# 4. Standardize scores and perform statistical testing
# 5. Generate annotated heatmap with clinical variables and p-values
# 6. Export as PDF and PNG for publication
#
# Requirements:
# - Clinical metadata file with required columns
# - TPM matrices from salmon quantification
# - R packages: immunedeconv, dplyr, ComplexHeatmap, tidyr, here, genefilter

# -----------------------------------------------------------------------------
# Package loading: Load required libraries with error handling
# -----------------------------------------------------------------------------
# Note: Ensure all packages are installed before running this script
# The "immunedeconv" package is not available on CRAN; install from GitHub using remotes:
# if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")
# remotes::install_github("icbi-lab/immunedeconv")
# The "dplyr", "tidyr" and "here" packages are available on CRAN:
# install.packages(c("dplyr", "tidyr", "here"))
# The "ComplexHeatmap" and "genefilter" packages are available on Bioconductor:
# if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
# BiocManager::install("ComplexHeatmap")
# BiocManager::install("genefilter")

required_packages <- c("immunedeconv", "dplyr", "ComplexHeatmap", "tidyr", "here", "genefilter")

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
DECONVOLUTION_METHOD <- "mcp_counter"  # Immune cell deconvolution method

# File paths (relative to project root)
CLINICAL_FILE <- here::here("Metadata", "RNAseq_sample_plan_SWISARC.txt")
TPM_COUNT_FILES <- c(
  here::here("TPM", "salmon.merged.gene_tpm_1.tsv"),
  here::here("TPM", "salmon.merged.gene_tpm_2.tsv"),
  here::here("TPM", "salmon.merged.gene_tpm_3.tsv")
)

# -----------------------------------------------------------------------------
# Reproducibility
# -----------------------------------------------------------------------------
set.seed(37)
# Ensure the working directory is set to the project root
setwd(here::here())

# -----------------------------------------------------------------------------
# Load clinical metadata and TPM expression matrices
# -----------------------------------------------------------------------------
# Load clinical metadata
clinical_metadata <- read.csv(CLINICAL_FILE, header = TRUE, sep = "\t", row.names = 1)

# -----------------------------------------------------------------------------
# Load TPM expression matrices (TPM-normalized RNA-seq gene counts) from multiple experimental batches
# -----------------------------------------------------------------------------
tpm_count_files <- setNames(TPM_COUNT_FILES, paste0("batch", 1:length(TPM_COUNT_FILES)))

# Function to read a count matrix for a given batch,
# handling possible duplicated gene ID columns or extra annotation columns
read_count_batch <- function(file) {
  dat <- read.csv(file, header = TRUE, sep = "\t", row.names = 1)
  # Remove "gene_name" column (if present), or a duplicate ID column (if present)
  if ("gene_name" %in% colnames(dat)) {
    dat <- dat[, !(colnames(dat) == "gene_name"), drop = FALSE]
  } else if (colnames(dat)[1] == rownames(dat)[1]) {
    dat <- dat[, -1, drop = FALSE]
  }
  dat
}

# Read TPM count matrices for each batch
tpm_count_matrix_list <- lapply(tpm_count_files, read_count_batch)

# -----------------------------------------------------------------------------
# Merge TPM expression matrices from all experimental batches (by gene ID) and filter samples
# -----------------------------------------------------------------------------
# Ensure all matrices have proper and comparable rownames (gene IDs)
# Use Reduce with merge, using row names as keys for all
merge_on_rownames <- function(x, y) {
  merged <- merge(x, y, by = "row.names", all = FALSE)
  rownames(merged) <- merged$Row.names
  merged <- merged[ , -1, drop = FALSE]
  return(merged)
}
tpm_merged_count_matrix <- Reduce(merge_on_rownames, tpm_count_matrix_list)

# Exclude poorly differentiated chordoma and EERT cases, retain only samples passing QC
samples_to_keep <- clinical_metadata$histological.subtype != "Poorly differenciated chordoma" & 
                   clinical_metadata$histological.subtype != "EERT"

# Subset TPM count matrix to keep only selected samples
tpm_filtered_count_matrix <- tpm_merged_count_matrix[, rownames(clinical_metadata)[samples_to_keep]]

# Subset clinical metadata to align with filtered samples
clinical_to_keep <- clinical_metadata[samples_to_keep, ]

# Fix gene name mismatches to ensure compatibility with MCPcounter's expected gene symbols
rownames(tpm_filtered_count_matrix)[rownames(tpm_filtered_count_matrix) == "SDPR"] <- "C"
rownames(tpm_filtered_count_matrix)[rownames(tpm_filtered_count_matrix) == "CXorf36"] <- "DIPK2B"

# -----------------------------------------------------------------------------
# Define color palettes for sample annotations and heatmaps
# -----------------------------------------------------------------------------
col_Transcriptomic_subtype <- c(
  "Proximal-like EpS" = "#3D426B",
  "Distal-like EpS" = "#CD534C"
)

col_Histological_subtype <- c(
  "Proximal EpS"  = "#2965e7",
  "Distal EpS"    = "#b3154f",
  "Hybrid EpS"    = "#b214ba"
)

col_list <- list(
  "Transcriptomic subtype"  = col_Transcriptomic_subtype,
  "Histological subtype"    = col_Histological_subtype
)

heatmap_colors <- colorRampPalette(c("navy", "white", "firebrick"))(255)

# -----------------------------------------------------------------------------
# Deconvolution and standardization of immune cell scores using MCPcounter
# -----------------------------------------------------------------------------
# Apply MCPcounter deconvolution to the TPM-normalized count matrix
mcp_score_matrix <- deconvolute(tpm_filtered_count_matrix, DECONVOLUTION_METHOD)

# Rename columns to include "cell_type" and sample AnnoIDs
colnames(mcp_score_matrix) <- c("cell_type", clinical_to_keep$AnnoID)

# Remove redundant "Monocyte" row (duplicate of "Macrophage/Monocyte")
mcp_score_matrix <- mcp_score_matrix[mcp_score_matrix$cell_type != "Monocyte", ]

# Standardize cell type names for readability and consistency
cell_type_replacements <- c(
  "T cells",
  "CD8 T cells",
  "Cytotoxic lymphocytes",
  "NK cells",
  "B lineage",
  NA,                              # Keep original name (Macrophage/Monocyte)
  "Myeloid dendritic cells",
  "Neutrophils",
  "Endothelial cells",
  "Fibroblasts"
)

indices_to_replace <- which(!is.na(cell_type_replacements))
mcp_score_matrix$cell_type[indices_to_replace] <- cell_type_replacements[indices_to_replace]

# -----------------------------------------------------------------------------
# Construct cell type by sample matrix and perform z-score normalization
# -----------------------------------------------------------------------------
# Define standardization function
standardize <- function(x) (x - mean(x)) / sd(x)

# Z-score standardization and matrix creation
matrix_fig1c <- mcp_score_matrix %>%
  select(-cell_type) %>%           # Remove cell type column
  apply(1, standardize) %>%        # Apply z-score standardization row-wise
  t()                              # Transpose to get cell types as rows, samples as columns
rownames(matrix_fig1c) <- mcp_score_matrix$cell_type

# -----------------------------------------------------------------------------
# Prepare clinical annotation and ordering for heatmap
# -----------------------------------------------------------------------------
# Order samples by transcriptomic subtype (descending) and histological subtype
order_idx <- with(clinical_to_keep, rev(order(Transcriptomic.subtype, histological.subtype)))

# Reorder matrix columns according to sample ordering
matrix_fig1c <- matrix_fig1c[, order_idx]

# Create annotation data frame for heatmap, ordered consistently with the matrix
clinical_to_keep_fig1c <- clinical_to_keep[order_idx, c("histological.subtype", "Transcriptomic.subtype")]
rownames(clinical_to_keep_fig1c) <- clinical_to_keep$AnnoID[order_idx]

# Rename annotation columns for display
colnames(clinical_to_keep_fig1c) <- c(
  "Histological subtype",
  "Transcriptomic subtype"
)

# Convert annotation columns to factors with set levels for consistent coloring
clinical_to_keep_fig1c$`Transcriptomic subtype` <- factor(
  clinical_to_keep_fig1c$`Transcriptomic subtype`,
  levels = c("Proximal-like EpS", "Distal-like EpS")
)
clinical_to_keep_fig1c$`Histological subtype` <- factor(
  clinical_to_keep_fig1c$`Histological subtype`,
  levels = c("Proximal EpS", "Distal EpS", "Hybrid EpS")
)

# -----------------------------------------------------------------------------
# Statistical testing and annotation (p-values for immune cells)
# -----------------------------------------------------------------------------
# Perform t-tests to compare immune cell scores between transcriptomic subtypes
ttest_fig1c <- rowttests(
  matrix_fig1c,
  factor(
    clinical_to_keep$Transcriptomic.subtype,
    levels = c("Proximal-like EpS", "Distal-like EpS")
  )
)

# Create row annotation showing p-values from t-tests
row_annotation_fig1c <- rowAnnotation(
  "P value" = anno_numeric(
    ttest_fig1c$p.value,
    rg = c(min(ttest_fig1c$p.value), 1),    # Range from minimum p-value to 1
    x_convert = function(x) -log10(x),       # Convert to -log10 scale
    labels_format = function(x) sprintf("%.2e", x)  # Format as scientific notation
  ),
  annotation_name_rot = 0                   # Keep annotation name horizontal
)

# Draw heatmap with clinical annotations
heatmap_fig1c <- ComplexHeatmap::pheatmap(
  matrix_fig1c,
  cluster_rows            = FALSE,
  cluster_cols            = FALSE,
  row_names_side          = "left",
  scale                   = "row",
  fontsize_col            = 10,
  fontsize_row            = 14,
  fontsize                = 12,
  annotation_col          = clinical_to_keep_fig1c,
  annotation_colors       = col_list,
  color                   = heatmap_colors,
  gaps_col                = 13,
  border_color            = NA,
  angle_col               = "45",
  treeheight_col          = 0,
  treeheight_row          = 0,
  show_colnames           = TRUE,
  show_rownames           = TRUE,
  heatmap_legend_param    = list(title = "Row Z-score")
)


# -----------------------------------------------------------------------------
# Export heatmap to PDF and PNG (Figure 1C)
# -----------------------------------------------------------------------------
# Save as PDF
pdf(file.path(here::here("Figures", "Figure_1C_mpc_counter.pdf")), width = 12, height = 4)
print(heatmap_fig1c+row_annotation_fig1c)
dev.off()

# Save as PNG
png(file.path(here::here("Figures", "Figure_1C_mpc_counter.png")), width = 3600, height = 1200, res = 300)
print(heatmap_fig1c+row_annotation_fig1c)
dev.off()