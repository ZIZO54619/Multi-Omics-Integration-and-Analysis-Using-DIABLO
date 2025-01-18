# Install BiocManager if not already installed
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# Install mixOmics package
BiocManager::install("mixOmics", update = FALSE, ask = FALSE)

# Set CRAN repository and install ragg package
options(repos = c(CRAN = "https://cran.r-project.org"))
install.packages("ragg", dependencies = TRUE, update = FALSE)

# Load Library
library(ragg)         # Graphics package
library(mixOmics)    # For multivariate analysis, including DIABLO
library(tools)        # For file path manipulation
library(ggplot2)     # For creating plots

