# SWISARC Manuscript Figures

This repository contains the R scripts and data processing pipelines used to generate figures for the SWISARC manuscript. The repository provides reproducible code for single-cell RNA sequencing analysis, immune cell deconvolution, survival analysis, and visualization of RNAseq epithelioid sarcoma data.

## 📋 Table of Contents

- [Overview](#overview)
- [Repository Structure](#repository-structure)
- [Requirements](#requirements)
- [Installation](#installation)
- [Usage](#usage)
- [Figure Descriptions](#figure-descriptions)
- [Data Requirements](#data-requirements)
- [Contributing](#contributing)
- [License](#license)

## 🔬 Overview

This repository contains computational pipelines for analyzing single-cell and bulk RNA sequencing data from epitheliod sarcoma samples. The analysis includes:

- **Single-cell RNA-seq analysis** with UMAP visualization and cell type annotation
- **Immune cell deconvolution** using MCPcounter
- **Hierarchical clustering** of transcriptomic subtypes
- **Survival analysis** with Kaplan-Meier curves
- **Differential gene expression** analysis and signature scoring

## 📁 Repository Structure

```
swisarc-manuscript-figures/
├── Main/                           # Main figure generation scripts
│   ├── Figure_1B_Hierarchical_clustering.R
│   ├── Figure_1C_mcp_counter.R
│   ├── Figure_1E_umaps.R
│   ├── Figure_1H_umaps.R
│   └── Figure_1K_Sup15A_survival.R
├── Figures/                        # Generated figure outputs (PDF/PNG)
├── UMAP/                          # Pre-computed UMAP embeddings
├── TPM/                           # TPM-normalized expression matrices
├── Counts/                        # Raw count matrices
└── README.md
```

## 🔧 Requirements

### R Packages

The following R packages are required and can be installed from CRAN or Bioconductor:

**CRAN packages:**
```r
install.packages(c(
  "dplyr", "ggplot2", "Seurat", "here", "survival",
  "ggpubr", "hues", "cowplot", "tidyr", "pheatmap",
  "vegan", "survminer"
))
```

**Bioconductor packages:**
```r
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
BiocManager::install(c(
  "ComplexHeatmap", "genefilter", "DESeq2"
))
```

**GitHub packages:**
```r
if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}
remotes::install_github("icbi-lab/immunedeconv")
```

### System Requirements

- **R version**: ≥ 4.0.0
- **Memory**: ≥ 16GB RAM recommended
- **Storage**: ≥ 5GB free space for data and outputs

## 🚀 Installation

1. **Clone the repository:**
```bash
git clone https://github.com/your-username/SWISARC-MANUSCRIPT-FIGURES.git
cd SWISARC-MANUSCRIPT-FIGURES
```

2. **Install required R packages** (see Requirements section above)

3. **Prepare data directories:**
```bash
mkdir -p Figures UMAP TPM Counts
```

## 📊 Usage

### Running Individual Figure Scripts

Each figure script is self-contained and can be run independently:

```bash
# Run from the repository root directory
Rscript Main/Figure_1B_Hierarchical_clustering.R
Rscript Main/Figure_1C_mcp_counter.R
Rscript Main/Figure_1E_umaps.R
Rscript Main/Figure_1H_umaps.R
Rscript Main/Figure_1K_Sup15A_survival.R
```

### Expected Outputs

Each script generates publication-ready figures in both PDF and PNG formats:

- **PDF files**: Vector format for publication (scalable)
- **PNG files**: High-resolution raster format for presentations

## 📈 Figure Descriptions

### Figure 1B - Hierarchical Clustering
- **Script**: `Figure_1B_Hierarchical_clustering.R`
- **Purpose**: Hierarchical clustering analysis of transcriptomic subtypes
- **Output**: `Figure_1B_hierarchical_clustering.pdf/png`

### Figure 1C - MCP Counter
- **Script**: `Figure_1C_mcp_counter.R`
- **Purpose**: Immune cell deconvolution heatmap using MCPcounter
- **Output**: `Figure_1C_mpc_counter.pdf/png`

### Figure 1E - UMAPs (All Cells)
- **Script**: `Figure_1E_umaps.R`
- **Purpose**: UMAP visualization colored by sample ID and cell type
- **Output**: `Figure_1E_umaps.pdf/png`

### Figure 1H - UMAPs (Tumor Cells)
- **Script**: `Figure_1H_umaps.R`
- **Purpose**: UMAP visualization of tumor cell subclusters
- **Output**: `Figure_1H_umaps.pdf/png`

### Figure 1K & Sup15A - Survival Analysis
- **Script**: `Figure_1K_Sup15A_survival.R`
- **Purpose**: Kaplan-Meier survival curves for EMT_5 signature
- **Outputs**:
  - `Figure_1K_survival_OS.pdf/png` (Overall Survival)
  - `Figure_Sup15A_survival_MFS.pdf/png` (Metastasis-Free Survival)

## 📂 Data Requirements

### Required Data Files

The scripts expect the following data files in their respective directories:

**UMAP/ directory:**
- `SWISARC_all_samples_nokeep_SCTransform_pca_33_0.5_umap.rda`
- `merge_tumor_harmony_SCTransform_pca_40_0.2.rda`

**TPM/ directory:**
- `salmon.merged.gene_tpm_1.tsv`
- `salmon.merged.gene_tpm_2.tsv`
- `salmon.merged.gene_tpm_3.tsv`

**Counts/ directory:**
- `salmon.merged.gene_counts_1.tsv`
- `salmon.merged.gene_counts_2.tsv`
- `salmon.merged.gene_counts_3.tsv`

**Metadata/ directory:**
- `RNAseq_sample_plan_SWISARC.txt`
- `Survival_table_SWISARC.txt`

### Data Format Requirements

- **Expression matrices**: Tab-separated files with genes as rows and samples as columns
- **Metadata files**: Tab-separated files with sample information
- **Seurat objects**: RDS/RDA files containing pre-computed UMAP embeddings

## 🔬 Analysis Workflow

### 1. Data Preprocessing
- Load and merge expression matrices from multiple experimental batches
- Filter samples based on quality control criteria
- Apply variance-stabilizing transformations

### 2. Single-cell Analysis
- Load pre-computed UMAP embeddings
- Annotate cell types and tumor clusters
- Generate publication-ready visualizations

### 3. Immune Cell Deconvolution
- Apply MCPcounter algorithm for immune cell scoring
- Perform statistical testing between groups
- Generate annotated heatmaps

### 4. Survival Analysis
- Calculate signature scores using differential gene expression
- Determine optimal cutpoints for survival stratification
- Generate Kaplan-Meier survival curves

## 🛠️ Troubleshooting

### Common Issues

1. **Package installation errors**: Ensure you have the latest version of R and all dependencies
2. **Memory issues**: Increase R memory limit or run scripts on a machine with more RAM
3. **File path errors**: Ensure all data files are in the correct directories as specified in the scripts

### Getting Help

If you encounter issues:
1. Check that all required data files are present
2. Verify R package installations
3. Review the script comments for specific requirements
4. Open an issue on GitHub with error details

## 📝 Script Documentation

Each script includes comprehensive documentation:
- **Metadata section**: Project information, workflow, and requirements
- **Configuration section**: Key parameters and file paths
- **Inline comments**: Detailed explanations of each analysis step
- **Export section**: Figure generation and saving procedures

## 🤝 Contributing

Contributions are welcome! Please:

1. Fork the repository
2. Create a feature branch
3. Make your changes with appropriate documentation
4. Test your changes thoroughly
5. Submit a pull request

## 📄 License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## 📞 Contact

For questions about this repository or the SWISARC project, please contact:

- **Author**: Leo Colmet Daage
- **Institution**: Institut Gustave Roussy (IGR)
- **Project**: SWISARC (Epithelioid Sarcoma Research)

## 📚 Citation

If you use this code in your research, please cite the SWISARC manuscript (citation details to be added upon publication).

---

**Note**: This repository contains research code for manuscript figure generation. Ensure you have appropriate permissions and follow institutional guidelines when using this code.
