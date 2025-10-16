#!/usr/bin/env Rscript

# -----------------------------------------------------------------------------
# Metadata
# -----------------------------------------------------------------------------
# Project: SWISARC manuscript figures
# Figure: 1B - Hierarchical clustering
# Author: Leo Colmet Daage
# Purpose: Build sample-to-sample distance heatmap with clinical annotations
# Date: 2024-02-07
#
# Workflow:
# 1. Load clinical metadata and RNA-seq count matrices from multiple batches
# 2. Filter samples (exclude poorly differentiated chordoma)
# 3. Merge count matrices and apply DESeq2 variance-stabilizing transformation
# 4. Select top variable genes for clustering
# 5. Compute sample-to-sample distances and hierarchical clustering
# 6. Generate annotated heatmap with clinical variables
# 7. Export as PDF and PNG for publication
#
# Requirements:
# - Clinical metadata file with required columns
# - Count matrices from salmon quantification
# - R packages: DESeq2, pheatmap, here, vegan, ggplot2

# -----------------------------------------------------------------------------
# Package loading: Load required libraries with error handling
# -----------------------------------------------------------------------------
# Note: Ensure all packages are installed before running this script
# The "pheatmap", "here", "vegan", and "ggplot2" packages are available on CRAN:
# install.packages(c("pheatmap", "here", "vegan", "ggplot2"))
# The "DESeq2" package is available on Bioconductor:
# if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
# BiocManager::install("DESeq2")

required_packages <- c("DESeq2", "pheatmap", "here", "vegan", "ggplot2")

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
MIN_COUNT_THRESHOLD <- 10          # Minimum total counts per gene for filtering
TOP_VARIABLE_PCT <- 0.25           # Percentage of most variable genes to use
CLUSTERING_METHOD <- "ward.D2"     # Hierarchical clustering method
DISTANCE_METHOD <- "pearson"       # Distance method for correlation

# File paths (relative to project root)
CLINICAL_FILE <- here::here("Metadata", "RNAseq_sample_plan_SWISARC.txt")
COUNT_FILES <- c(
  here::here("Counts", "salmon.merged.gene_counts_1.tsv"),
  here::here("Counts", "salmon.merged.gene_counts_2.tsv"),
  here::here("Counts", "salmon.merged.gene_counts_3.tsv")
)

# -----------------------------------------------------------------------------
# Reproducibility
# -----------------------------------------------------------------------------
set.seed(37)
# Ensure the working directory is set to the project root
setwd(here::here())

# -----------------------------------------------------------------------------
# Load clinical metadata and RNA-seq count matrices
# -----------------------------------------------------------------------------
# Load clinical metadata
clinical_metadata <- read.csv(CLINICAL_FILE, header = TRUE, sep = "\t", row.names = 1)

# -----------------------------------------------------------------------------
# Load expression matrices (raw RNA-seq gene counts) from multiple experimental batches
# -----------------------------------------------------------------------------
count_files <- setNames(COUNT_FILES, paste0("batch", 1:length(COUNT_FILES)))

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

# Read count matrices for each batch
count_matrix_list <- lapply(count_files, read_count_batch)

# -----------------------------------------------------------------------------
# Merge expression count matrices from all experimental batches (by gene ID) and filter samples
# -----------------------------------------------------------------------------
# Ensure all matrices have proper and comparable rownames (gene IDs)
# Use Reduce with merge, using row names as keys for all
merge_on_rownames <- function(x, y) {
  merged <- merge(x, y, by = "row.names", all = FALSE)
  rownames(merged) <- merged$Row.names
  merged <- merged[ , -1, drop = FALSE]
  return(merged)
}
merged_count_matrix <- Reduce(merge_on_rownames, count_matrix_list)

# Exclude poorly differentiated chordoma samples and keep those passing QC
samples_to_keep <- clinical_metadata$histological.subtype != "Poorly differenciated chordoma"
filtered_count_matrix <- merged_count_matrix[, rownames(clinical_metadata)[samples_to_keep]]    
clinical_to_keep <- clinical_metadata[samples_to_keep, ]

# -----------------------------------------------------------------------------
# Define color palettes for sample annotations and heatmaps
# -----------------------------------------------------------------------------
col_Histological_subtype <- c(
  "EERT"          = "#FFFE11",
  "Proximal EpS"  = "#2965e7",
  "Distal EpS"    = "#b3154f",
  "Hybrid EpS"    = "#b214ba"
)
col_Tumor_status <- c(
  "Pre-treated primary" = "#26A2F0",
  "Untreated primary"   = "#DB40C9",
  "Local recurrence"    = "#BD672B",
  "Metastasis"          = "#099B8C"
)

col_Age_at_diagnosis <- c(
  "[16-25[" = "#FEFB01",
  "[25-37[" = "#CEFB02",
  "[37-48[" = "#00ED01",
  "[48-72[" = "#029e19"
)

col_Primary_anatomical_site <- c(
  "Distal lower limb" = "#cb4f42",
  "Distal upper limb" = "#c65c8a",
  "Head and neck" = "#c98443",
  "Perineal/Genital" = "#9b9c3b",
  "Proximal lower limb" = "#648ace",
  "Proximal upper limb" = "#58a865",
  "Internal trunk" = "#a361c7"
)

col_Transcriptomic_subtype <- c(
  "EERT" = "#b4943e",
  "Proximal-like EpS" = "#3D426B",
  "Distal-like EpS" = "#CD534C"
)

col_list <- list(
  "Transcriptomic subtype"  = col_Transcriptomic_subtype,
  "Histological subtype"    = col_Histological_subtype,
  "Tumor status"            = col_Tumor_status,
  "Age at diagnosis (y)"    = col_Age_at_diagnosis,
  "Primary anatomical site" = col_Primary_anatomical_site
)

heatmap_colors <- colorRampPalette(c("navy", "white", "firebrick"))(255)

# -----------------------------------------------------------------------------
# DESeq2 object creation and variance stabilization
# -----------------------------------------------------------------------------
# Create DESeq2 dataset using an intercept-only design (unsupervised workflow)
dds <- DESeqDataSetFromMatrix(round(filtered_count_matrix), colData = clinical_to_keep, design = ~1)
colnames(dds) <- colData(dds)$AnnoID

# Low-count filtering and Variance-stabilizing transform 
keep_high_count <- rowSums(counts(dds)) >= MIN_COUNT_THRESHOLD
dds_filtered <- dds[keep_high_count,]
vsd_filtered <- vst(dds_filtered, blind = TRUE)

# -----------------------------------------------------------------------------
# Select top most variable genes for downstream clustering/heatmap
# -----------------------------------------------------------------------------
# Select the top variable genes based on variance
n_genes <- nrow(assay(vsd_filtered))
top_n <- ceiling(n_genes * TOP_VARIABLE_PCT)
gene_variances <- rowVars(assay(vsd_filtered))
top_genes_idx <- head(order(gene_variances, decreasing = TRUE), top_n)
vsd_top_variable <- vsd_filtered[top_genes_idx, ]

# -----------------------------------------------------------------------------
# Compute sample-to-sample distances and hierarchical clustering (Fig 1B)
# -----------------------------------------------------------------------------
# Extract expression values (variance-stabilized counts) for the top variable genes
normalized_vsd_top_variable <- assay(vsd_top_variable)

# Use 1 - correlation as sample distance (computed on top-variable genes)
sampleDists <- as.dist(1 - cor(normalized_vsd_top_variable, method = DISTANCE_METHOD))
sample_hclust <- hclust(sampleDists, method = CLUSTERING_METHOD)
sample_dendrogram <- as.dendrogram(sample_hclust)
sample_dendrogram <- reorder(sample_dendrogram, wts = rev(sample_hclust$order))
sampleDistMatrix <- as.matrix(sampleDists)

# Ensure consistent histotype order in annotations and set rownames to AnnoID
clinical_to_keep$histological.subtype <- factor(
  clinical_to_keep$histological.subtype,
  levels = c("EERT", "Proximal EpS", "Distal EpS", "Hybrid EpS")
)

rownames(clinical_to_keep) <- clinical_to_keep$AnnoID
clinical_to_keep$AnnoID <- NULL

# Rename columns of clinical_to_keep to match annotation labels for heatmap
colnames(clinical_to_keep) <- c(
  "Primary anatomical site",
  "Age at diagnosis (y)",
  "Tumor status",
  "Histological subtype",
  "Transcriptomic subtype"
)

# Draw heatmap with clinical annotations
heatmap_fig1b <- pheatmap(
  sampleDistMatrix,
  cluster_rows       = as.hclust(sample_dendrogram),
  cluster_cols       = as.hclust(sample_dendrogram),
  fontsize_col       = 16,
  fontsize           = 14,
  annotation_col     = clinical_to_keep,
  color              = heatmap_colors,
  annotation_colors  = col_list,
  treeheight_col     = 100,
  treeheight_row     = 0,
  cellheight         = 20,
  cellwidth          = 20,
  show_colnames      = TRUE,
  show_rownames      = FALSE,
  angle_col          = "45",
  border_color       = NA,
  na_col             = "grey"
)

# -----------------------------------------------------------------------------
# Export heatmap to PDF and PNG (Figure 1B)
# -----------------------------------------------------------------------------
# Save as PDF
pdf(file.path(here::here("Figures", "Figure_1B_Hierarchical_clustering.pdf")), width = 14, height = 14)
print(heatmap_fig1b)
dev.off()

# Save as PNG
png(file.path(here::here("Figures", "Figure_1B_Hierarchical_clustering.png")), width = 4200, height = 4200, res = 300)
print(heatmap_fig1b)
dev.off()