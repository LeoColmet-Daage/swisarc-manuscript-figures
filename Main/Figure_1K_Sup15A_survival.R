#!/usr/bin/env Rscript

# -----------------------------------------------------------------------------
# Metadata
# -----------------------------------------------------------------------------
# Project: SWISARC manuscript figures
# Figure: 1K & Sup15A - Survival analysis
# Author: Leo Colmet Daage
# Purpose: Generate survival analysis plots for EMT_5 cluster signature
# Date: 2024-02-07
#
# Workflow:
# 1. Load Seurat object with UMAP embeddings and cluster annotations
# 2. Load survival metadata and clinical information
# 3. Load and merge RNA-seq count matrices from multiple batches
# 4. Create DESeq2 object and perform variance-stabilizing transformation
# 5. Run differential gene expression analysis to identify cluster signatures
# 6. Extract top 100 marker genes for each tumor cluster as signatures
# 7. Calculate signature scores using Z-score arithmetic mean
# 8. Merge survival data with signature scores
# 9. Perform survival analysis with optimal cutpoint determination
# 10. Generate Kaplan-Meier survival curves for overall survival and metastasis-free survival
# 11. Export figures as PDF and PNG for publication
#
# Requirements:
# - Pre-computed Seurat object with UMAP embeddings
# - Survival metadata table with clinical follow-up data
# - RNA-seq count matrices from multiple experimental batches
# - R packages: dplyr, ggplot2, Seurat, here, survival, DESeq2, survminer


# -----------------------------------------------------------------------------
# Package loading: Load required libraries with error handling
# -----------------------------------------------------------------------------
# Note: Ensure all packages are installed before running this script
# The "dplyr", "ggplot2", "Seurat", "here", "survival", "DESeq2" and "survminer" packages are available on CRAN:
# install.packages(c("dplyr","ggplot2", "Seurat", "here", "survival", "DESeq2", "survminer"))

required_packages <- c("dplyr", "ggplot2", "Seurat", "here", "survival", "DESeq2", "survminer")

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

# File paths (relative to project root)
UMAP_FILE <- here::here("UMAP", "merge_tumor_harmony_SCTransform_pca_40_0.2.rda")
SURVIVAL_FILE <- here::here("Metadata", "Survival_table_SWISARC.txt")
CLINICAL_FILE <- here::here("Metadata", "RNAseq_sample_plan_SWISARC.txt")
COUNT_FILES <- c(
  here::here("Counts", "salmon.merged.gene_counts_1.tsv"),
  here::here("Counts", "salmon.merged.gene_counts_2.tsv"),
  here::here("Counts", "salmon.merged.gene_counts_3.tsv")
)

# -----------------------------------------------------------------------------
# Reproducibility
# -----------------------------------------------------------------------------
set.seed(27112022)
# Ensure the working directory is set to the project root
setwd(here::here())

# -----------------------------------------------------------------------------
# Load pre-computed UMAP data from Seurat object
# -----------------------------------------------------------------------------
# Load the Seurat object containing UMAP embeddings and metadata
load(UMAP_FILE)

# Map cluster numbers to biologically meaningful tumor cluster names
# Based on marker gene expression and known tumor cluster characteristics
levels(sobj@active.ident)[1]<-"Cell migration_0"
levels(sobj@active.ident)[levels(sobj@active.ident)=="1"]<-"Apoptotic signaling_1"
levels(sobj@active.ident)[levels(sobj@active.ident)=="2"]<-"Regulation of growth_2"
levels(sobj@active.ident)[levels(sobj@active.ident)=="3"]<-"Immune response_3"
levels(sobj@active.ident)[levels(sobj@active.ident)=="4"]<-"Ribosomal_IFN response_4"
levels(sobj@active.ident)[levels(sobj@active.ident)=="5"]<-"EMT_5"
levels(sobj@active.ident)[levels(sobj@active.ident)=="6"]<-"Cycling_cells_6"
levels(sobj@active.ident)[levels(sobj@active.ident)=="7"]<-"LncRNA_7"
levels(sobj@active.ident)[levels(sobj@active.ident)=="8"]<-"Circulatory system development_8"
levels(sobj@active.ident)[levels(sobj@active.ident)=="9"]<-"Response to metal ion_9"

sobj@meta.data$active.ident<-sobj@active.ident

# -----------------------------------------------------------------------------
# Load survival metadata and preprocess
# -----------------------------------------------------------------------------
# Load survival metadata
survival_metadata <- read.csv(SURVIVAL_FILE, sep = "\t", dec = ",")

# Convert dates to Date objects
survival_metadata$Date.of.diagnosis <- as.Date(survival_metadata$Date.of.diagnosis, format = "%d/%m/%Y")
survival_metadata$Date.of.last.FU <- as.Date(survival_metadata$Date.of.last.FU, format = "%d/%m/%Y")
survival_metadata$Date.of.metastasis.last.FU <- as.Date(survival_metadata$Date.of.metastasis.last.FU, format = "%d/%m/%Y")

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
# DESeq2 object creation and variance stabilization
# -----------------------------------------------------------------------------
# Create DESeq2 dataset using an intercept-only design (unsupervised workflow)
dds <- DESeqDataSetFromMatrix(round(filtered_count_matrix), colData = clinical_to_keep, design = ~1)
colnames(dds) <- colData(dds)$AnnoID

# Variance-stabilizing transformation of the count data
vsd <- vst(dds, blind = TRUE)

# -----------------------------------------------------------------------------
# Run differential gene expression to get top 100 signatures
# -----------------------------------------------------------------------------
markers_merge <- FindAllMarkers(
  object = sobj, 
  only.pos = TRUE,
  min.pct = 0.25, 
  logfc.threshold = 0.25
)

# Delete redundant genes between clusters
markers_merge_unique <- markers_merge[!duplicated(markers_merge$gene),]

# -----------------------------------------------------------------------------
# Get top 100 unique marker genes for each tumor cluster (as signature genes)
# -----------------------------------------------------------------------------
# Extract gene and cluster columns
signatures <- markers_merge_unique[, c("gene", "cluster")]

# Split by cluster, extract top 100 (or fewer if not enough) genes for each cluster
signatures_list <- split(signatures$gene, signatures$cluster)
signatures_top100_list <- lapply(signatures_list, function(x) head(unique(x), 100))

# Create a data.frame filling missing values with NA, each column is a cluster signature
max_length <- max(sapply(signatures_top100_list, length))
signatures_top100 <- as.data.frame(lapply(signatures_top100_list, function(x) {
  length(x) <- max_length
  x
}))
colnames(signatures_top100) <- names(signatures_top100_list)

# -----------------------------------------------------------------------------
# vsd transformed data -> Z score -> arithmetic mean
# -----------------------------------------------------------------------------
# Extract vsd transformed data as a data.frame
vsd_matrix <- as.data.frame(assay(vsd))
colnames(vsd_matrix) <- colData(vsd)$AnnoID

# For each tumor cluster, extract the top 100 marker genes, z-score them (row-wise), and calculate the mean score per sample
# Get scores for each cluster (column of signatures_top100)
cluster_scores <- lapply(seq_len(ncol(signatures_top100)), function(i) {
  genes <- signatures_top100[, i]
  genes <- genes[!is.na(genes)] # Remove NAs (in case less than 100)
  gene_data <- vsd_matrix[rownames(vsd_matrix) %in% genes, , drop = FALSE]
  # Row-wise z-score, replace NA from constant genes with zero
  z_data <- t(scale(t(gene_data)))
  z_data[is.na(z_data)] <- 0
  # Calculate mean z-score per sample
  score <- colMeans(z_data)
  # Return as data.frame with correct cluster name
  cluster_name <- paste0("top100_", names(signatures_top100)[i])
  setNames(data.frame(score), cluster_name)
})

# Merge scores into one data frame
top100_score_df <- do.call(cbind, cluster_scores)
top100_score_df$Patient <- colnames(vsd_matrix)

# -----------------------------------------------------------------------------
# Survival analysis
# -----------------------------------------------------------------------------
# Helper function to process survival data
process_survival_data <- function(data, time_start, time_end, status_col, exclude_types = c("ECRT", "Chordoma")) {
  df <- data %>%
    select(Patient, Subtype) %>%
    mutate(
      time = as.numeric(data[[time_end]] - data[[time_start]]),
      month = round(time / 30.417, 0),
      year = round(time / 365.25, 2),
      status = as.numeric(data[[status_col]])
    ) %>%
    filter(complete.cases(.), !Subtype %in% exclude_types)
  df
}

# Define survival input parameter list
survival_inputs <- list(
  OS  = list(time_start = "Date.of.diagnosis", time_end = "Date.of.last.FU",              status_col = "OS"),
  MFS = list(time_start = "Date.of.diagnosis", time_end = "Date.of.metastasis.last.FU",   status_col = "MFS")
)

# Process, join and store all survival datasets
df_list   <- list()
full_list <- list()
for (name in names(survival_inputs)) {
  params <- survival_inputs[[name]]
  df     <- process_survival_data(survival_metadata, params$time_start, params$time_end, params$status_col)
  full   <- inner_join(df, top100_score_df, by = "Patient") %>%
              mutate(year = as.numeric(year), month = as.numeric(month))
  df_list[[name]]   <- df
  full_list[[name]] <- full
}

# -----------------------------------------------------------------------------
# Helper function to calculate p-values for survival analysis
# -----------------------------------------------------------------------------
calculate_survival_pvalue <- function(fit, data, time_var = "year", event_var = "status", group_var) {
  # Perform log-rank test using survdiff
  surv_diff <- survdiff(
    as.formula(fit$call),
    data = data
  )
  
  # Extract p-value from chi-square test
  p_value <- 1 - pchisq(surv_diff$chisq, length(surv_diff$n) - 1)
  return(p_value)
}

# -----------------------------------------------------------------------------
# Find optimal cutpoint and fit survival models for each signature (OS and MFS)
# -----------------------------------------------------------------------------
cutpoint_category <- list()
fit_list <- list()
pvalue_list <- list()
for (outcome in c("OS", "MFS")) {
  cutpoint_category[[outcome]] <- list()
  fit_list[[outcome]] <- list()
  pvalue_list[[outcome]] <- list()
  score_cols <- colnames(full_list[[outcome]][, grepl("top100", names(full_list[[outcome]]))])
  for (i in score_cols) {
    cutpoint <- surv_cutpoint(
      full_list[[outcome]],
      time = "year",
      event = "status",
      variables = i,
      minprop = 0.1,
      progressbar = TRUE
    )
    cutpoint_category[[outcome]][[i]] <- surv_categorize(cutpoint, variables = NULL, labels = c("low", "high"))
    
    fit <- survfit(
      Surv(year, status) ~ cutpoint_category[[outcome]][[i]][, i],
      data = cutpoint_category[[outcome]][[i]]
    )
    # Clean up strata names
    names(fit$strata) <- sub("cutpoint_category\\[\\[outcome\\]\\]\\[\\[i\\]\\]\\[\\, i\\]\\=", "", names(fit$strata))
    fit_list[[outcome]][[i]] <- fit
    
    # Calculate p-value using log-rank test
    pvalue_list[[outcome]][[i]] <- calculate_survival_pvalue(
      fit = fit,
      data = cutpoint_category[[outcome]][[i]],
      group_var = i
    )
  }
}

# -----------------------------------------------------------------------------
# Plot EMT_5 survival curves
# -----------------------------------------------------------------------------
# Helper function to plot survival curves for a given outcome and signature
plot_signature_survival <- function(
  fit_list, 
  cutpoint_category, 
  pvalue_list,
  outcome = c("OS", "MFS"), 
  signature = "top100_EMT_5", 
  ylab_str = "Survival rate (%)",
  add_pvalue = FALSE,
  annotate_text = NULL,
  annotate_pos = c(x = 1.6, y = 0.03)
) {
  outcome <- match.arg(outcome)
  
  # y-label depending on outcome
  if (missing(ylab_str) || is.null(ylab_str)) {
    ylab_str <- ifelse(outcome == "OS", "Overall survival rate (%)", "Metastasis-free survival rate (%)")
  }
  
  cluster_name <- gsub("top100_", "", signature)
  
  # Get the calculated p-value
  calculated_pvalue <- pvalue_list[[outcome]][[signature]]
  
  # Use calculated p-value if no custom annotation is provided
  if (is.null(annotate_text)) {
    annotate_text <- bquote(italic("P") * "=" * .(sprintf("%.3f", calculated_pvalue)))
  }
  
 # Plot survival curves
  plot <- ggsurvplot(
    fit_list[[outcome]][[signature]],
    data = cutpoint_category[[outcome]][[signature]],
    risk.table = TRUE,
    conf.int = TRUE,
    pval = add_pvalue,
    break.time.by = 2,
    break.y = 0.2,
    fontsize = 6,
    legend.title = "Patient subgroup",
    tables.theme = theme(
      axis.text.x = element_text(size = 18),
      axis.text.y = element_text(size = 18),
      axis.title.x = element_text(size = 20),
      title = element_text(size = 18)
    ),
    title = cluster_name
  ) + 
    ylab(ylab_str)

  # Plot survival curves with custom theme
  plot$plot <- plot$plot +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = 20), 
      axis.text.x = element_text(size = 18), 
      axis.text.y = element_text(size = 18),
      axis.title.x = element_blank(),
      axis.title.y = element_text(size = 20),
      legend.text = element_text(size = 16),
      legend.title = element_text(size = 18),
      legend.direction = "vertical",
      legend.title.position = "left"
    ) +
    scale_y_continuous(
      breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1),
      labels = c("0", "20", "40", "60", "80", "100")
    ) +
    scale_color_manual(
      labels = c(
        bquote(italic(.(cluster_name)) * "-high"),
        bquote(italic(.(cluster_name)) * "-low")
      ),
      values = c("#F8766D", "#00BFC4")
    ) +
    scale_fill_manual(
      labels = c(
        bquote(italic(.(cluster_name)) * "-high"),
        bquote(italic(.(cluster_name)) * "-low")
      ),
      values = c("#F8766D", "#00BFC4"),
      drop = FALSE
    )

  # Plot risk table with custom theme
  plot$table <- plot$table +
    ylab("") +
    xlab("Time after diagnosis (years)") +
    theme(
      plot.title = element_text(size = 18),
      axis.text.y.left = element_text(vjust = 0.25)
    ) +
    scale_y_discrete(
      labels = c(
        paste0("<br><i>", cluster_name, "</i>-low"),
        paste0("<br><i>", cluster_name, "</i>-high")
      )
    )
  
  # Optionally annotate with p-value or other text
  if (!is.null(annotate_text)) {
    plot$plot <- plot$plot +
      annotate(
        "text",
        x = annotate_pos['x'],
        y = annotate_pos['y'],
        label = annotate_text,
        size = 8
      )
  }
  
  return(plot)
}

# Generate survival plots for EMT_5 signature for OS and MFS
plot_OS_EMT5  <- plot_signature_survival(
  fit_list = fit_list,
  cutpoint_category = cutpoint_category,
  pvalue_list = pvalue_list,
  outcome = "OS",
  signature = "top100_EMT_5",
  ylab_str = "Overall survival rate (%)"
)
plot_MFS_EMT5 <- plot_signature_survival(
  fit_list = fit_list,
  cutpoint_category = cutpoint_category,
  pvalue_list = pvalue_list,
  outcome = "MFS",
  signature = "top100_EMT_5",
  ylab_str = "Metastasis-free survival rate (%)"
)

# -----------------------------------------------------------------------------
# Export figures to PDF and PNG (Figure 1K & Sup15A)
# -----------------------------------------------------------------------------
# Export each plot separately to PDF
pdf(file.path(here::here("Figures", "Figure_1K_survival_OS.pdf")), width = 7, height = 7)
print(plot_OS_EMT5)
dev.off()

pdf(file.path(here::here("Figures", "Figure_Sup15A_survival_MFS.pdf")), width = 7, height = 7)
print(plot_MFS_EMT5)
dev.off()

# Export each plot separately to PNG
png(file.path(here::here("Figures", "Figure_1K_survival_OS.png")), width = 2000, height = 2000, res = 300)
print(plot_OS_EMT5)
dev.off()

png(file.path(here::here("Figures", "Figure_Sup15A_survival_MFS.png")), width = 2000, height = 2000, res = 300)
print(plot_MFS_EMT5)
dev.off()