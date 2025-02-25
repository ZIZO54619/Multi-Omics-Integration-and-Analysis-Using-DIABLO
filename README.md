# DIABLO Multi-Omics Analysis

## Overview
This project performs multi-omics data integration using the `mixOmics` package in R. It follows a step-by-step approach to preprocess, analyze, and visualize multi-omics datasets using the DIABLO (Data Integration Analysis for Biomarker discovery using Latent variable approaches for Omics studies) framework.

## Prerequisites
Ensure you have the following software and dependencies installed:
- **R (>= 4.0)**
- **Bioconductor**
- Required R packages:
  - `mixOmics`
  - `ggplot2`
  - `tools`
  - `ragg`

## Installation
Run the following commands in R to install the required packages:

```r
# Install BiocManager if not already installed
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# Install mixOmics package
BiocManager::install("mixOmics", update = FALSE, ask = FALSE)

# Install ragg package
install.packages("ragg", dependencies = TRUE, update = FALSE)
```

## Data Preparation
1. **Set the working directory:** The dataset should be stored in `../Dataset/TCGA/`.
2. **Ensure data format:** Each subdirectory should contain CSV files with omics data.
3. **Data loading process:** The script automatically reads the CSV files, processes them into matrices, and extracts subtype labels if available.

## Analysis Workflow
### 1. Load Required Libraries
```r
library(ragg)
library(mixOmics)
library(tools)
library(ggplot2)
```

### 2. Load Data
- Reads multi-omics datasets from CSV files
- Converts them into matrices
- Extracts subtype information

### 3. Define DIABLO Model Components
- Extracts training data (`X`), subtype data (`Y`)
- Defines a **design matrix** to control correlation strength between data types

### 4. Perform Unsupervised PLS Analysis
- Computes **Partial Least Squares (PLS)** correlations between different omic layers

### 5. Build DIABLO Model
- Constructs a **block PLS-DA model**
- Visualizes model structure and performance

### 6. Model Tuning and Validation
- Performs **cross-validation** to assess model accuracy
- Selects optimal **number of components (ncomp)**
- Tunes **keepX parameters** (number of selected features per omic layer)

### 7. Export Results
- Saves **loadings for each omic dataset** as CSV files
- Prints **design matrix** and important selected variables

### 8. Visualization
- **Sample plots:** `plotDiablo()`, `plotIndiv()`, `plotArrow()`
- **Variable plots:** `plotVar()`, `circosPlot()`, `network()`
- **Loading plots:** `plotLoadings()`

## Running the Analysis
To execute the workflow, simply run the R script in an R environment:
```r
source("your_script.R")
```
Make sure the dataset is correctly placed in the specified directory before running the script.

## Output
- Model performance evaluation
- Selected variables from each omic dataset
- Visualizations of sample and variable distributions
- Loadings exported as CSV files

## Contact
For questions or issues, please reach out to the project contributors or refer to the `mixOmics` documentation: [https://mixomics.org](https://mixomics.org).

