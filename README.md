Spatial Transcriptomic Analysis of Visium Human Colorectal Cancer (FFPE)

Overview

This repository contains the R script and workflow for spatial transcriptomic analysis of Visium Human Colorectal Cancer (FFPE) data from 10x Genomics. The project aims to:
Perform quality control and normalization of spatial transcriptomics data.
Apply dimensionality reduction techniques (PCA, UMAP, t-SNE).
Identify spatially variable genes and marker genes.
Visualize gene expression across tissue architecture.

Dataset
Source: 10x Genomics Visium platform
Sample Type: Human Colorectal Cancer, FFPE (11 mm Capture Area)
Format: Spatial gene expression matrix, tissue images, and metadata files.

Folder Structure

├── data/                 # Raw 10x Genomics data (filtered_feature_bc_matrix, spatial)
├── scripts/              # R scripts for analysis
├── results/              # Outputs: plots, marker genes, processed data
├── Visium_Analysis.R     # Main R script for spatial transcriptomics analysis
└── README.md             # Project documentation

Installation

Ensure you have R (>= 4.0.0) installed.

# Install required packages
install.packages(c("Seurat", "ggplot2", "dplyr", "patchwork"))
BiocManager::install(c("TENxVisiumData", "SpatialExperiment"))

# Load libraries
library(Seurat)
library(SpatialExperiment)
library(ggplot2)
library(dplyr)
library(patchwork)

Analysis Workflow

Load Data: Import 10x Genomics spatial transcriptomic data.
Quality Control: Filter low-quality spots based on gene count and mitochondrial content.
Normalization: Apply SCTransform for normalization.
Dimensionality Reduction: Perform PCA, UMAP, and t-SNE.
Clustering: Identify clusters and spatial domains.
Marker Gene Identification: Discover key genes defining clusters.
Visualization: Spatial plots for clusters and gene expression.



Contributing

Contributions are welcome! Please fork the repository, make changes, and submit a pull request.

License

This project is licensed under the MIT License.
