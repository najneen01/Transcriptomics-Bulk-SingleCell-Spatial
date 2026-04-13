 # Single-Cell Transcriptomic Analysis using Seurat (R)

## Project Overview

This project performs single-cell RNA sequencing (scRNA-seq) analysis using the Seurat package in R. The workflow includes data preprocessing, quality control, dimensionality reduction, clustering, and marker visualization.

The analysis is based on publicly available data from GEO.

## Data Source
Dataset: GSE275330
Source: NCBI Gene Expression Omnibus (GEO)
Link: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE275330

## Tools & Packages

The analysis was conducted in RStudio using the following packages:

Seurat
dplyr
Matrix
ggplot2
patchwork

## Workflow Summary
### 1. Data Loading
Imported matrix, barcodes, and features files
Constructed Seurat object
### 2. Quality Control (QC)
Filtered cells based on:
Number of features (genes)
Total counts
Mitochondrial gene percentage
### 3. Normalization & Feature Selection
Log normalization
Identification of highly variable genes
### 4. Dimensionality Reduction
Principal Component Analysis (PCA)
t-SNE and UMAP visualization
### 5. Clustering
Cell clustering based on selected principal components
### 6. Marker Analysis
Visualization of:
Neural markers (e.g., MKI67, NES, DCX)
Lymphoma markers (DLBCL, MCL)
