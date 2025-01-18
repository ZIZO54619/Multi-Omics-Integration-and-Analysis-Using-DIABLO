---
title: 'Diaplo Model'
author: "El_zoz"
date: "2024-08-13"
output: html_document
---

## Install and Load Required Packages

-   This section installs and loads the necessary R packages for the analysis.

-   Package installation is done only once per package on your system.

```{r install_packages, message=FALSE, warning=FALSE}
# Install BiocManager if not already installed
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# Install mixOmics package
BiocManager::install("mixOmics", update = FALSE, ask = FALSE)

# Set CRAN repository and install ragg package
options(repos = c(CRAN = "https://cran.r-project.org"))
install.packages("ragg", dependencies = TRUE, update = FALSE)
```

## Load the libraries.

-   This needs to be done in every R session where you want to use the packages.

```{r load_libraries, message=FALSE, warning=FALSE}
library(ragg)         # Graphics package
library(mixOmics)    # For multivariate analysis, including DIABLO
library(tools)        # For file path manipulation
library(ggplot2)     # For creating plots
```

## Load and Prepare Data

-   This section loads the data from CSV files, preprocesses it, and structures it for the DIABLO model.

```{r load_data, message=FALSE, warning=FALSE}
# Set the working directory
setwd(r"{D:\Multi-Omics\TCGA}")

# List directories in the working directory
folders <- list.dirs(path = ".", full.names = TRUE, recursive = FALSE)
data_set <- list() # Initialize an empty list to store data

# Loop through each folder (representing a data type)
for (folder in folders) {
  folder_name <- basename(folder)  # Extract folder name
  folder_files <- list.files(path = folder, pattern = "\\.csv$", full.names = TRUE)  # List all CSV files in the folder
  folder_data <- list()  # Initialize list to hold data from each file
  subtype <- NULL # Initialize the subtype vector
  
  # Loop through each CSV file
  for (file in folder_files) {
    df <- read.csv(file, header = TRUE, row.names = 1) # Read csv data into dataframe
    name <- tools::file_path_sans_ext(basename(file)) # Extract file name without extension
    
    # Check if there is a 'Label' column, if so, convert to factor and remove it from dataframe
    if ("Label" %in% colnames(df)) {
        Label <- df$Label 
        subtype <- factor(Label)
        df <- df[, -which(names(df) == "Label")]
        df <- sapply(df, as.numeric) # Convert all values in df to numeric
    }
    row_names <- rownames(df)
    df_matrix <- as.matrix(df)
    rownames(df_matrix) <- row_names # set row names for matrix
    folder_data[[name]] <- df_matrix # add matrix to folder_data list
  }
  folder_data$subtype <- subtype # add subtype factor
  data_set[[folder_name]] <- folder_data # add folder data to main data_set list
}

names(data_set) # Print the names of the data_set list
```

## Extract the training data from the loaded data set.

```{r extract_training_data}
# Determine the number of omic views
number_of_view <- length(data_set$data.train) - 1
X <- list() # Initialize a list to store the omic data
# Loop through each view and store it in list
for (i in 1:number_of_view) {
  view_name <- names(data_set$data.train[i])
  X[[view_name]] <- data_set$data.train[[view_name]]
}

names(X) # Print the names of the views
```

## Extract and summarize the subtype data from the training set.

```{r extract_subtype_data}
Y <- data_set$data.train$subtype # Get subtype factor
summary(Y) # Provide a summary of the subtype factor
```

## Define the design matrix for DIABLO, controlling the level of supervision.

```{r define_design_matrix}
# Set design matrix for DIABLO model, this supervises the correlation. A higher value leads to stronger supervision.
design <- matrix(0.1, ncol = length(X), nrow = length(X),
                 dimnames = list(names(X), names(X)))
# Set diagonal to zero as self-connections are not considered
diag(design) <- 0
design # Print the design matrix
```

## This section explores the correlation structure between omic views using unsupervised PLS.

```{r unsupervised_pls}
# Unsupervised PLS between mRNA and protein
res1.pls.tcga <- pls(X$mrna, X$protein, ncomp = 1)
cor(res1.pls.tcga$variates$X, res1.pls.tcga$variates$Y) # Check correlation of components

# Unsupervised PLS between mRNA and miRNA
res2.pls.tcga <- pls(X$mrna, X$mirna, ncomp = 1)
cor(res2.pls.tcga$variates$X, res2.pls.tcga$variates$Y) # Check correlation of components

# Unsupervised PLS between miRNA and protein
res3.pls.tcga <- pls(X$mirna, X$protein, ncomp = 1)
cor(res3.pls.tcga$variates$X, res3.pls.tcga$variates$Y) # Check correlation of components
```

## Build the DIABLO model using the prepared data and specified parameters.

```{r build_diablo_model}
# Create a block.plsda model
diablo.tcga <- block.plsda(X, Y, ncomp = 5, design = design)
plot(diablo.tcga) # plot of the diablo model
```

## Evaluate the DIABLO model performance through cross-validation.

```{r evaluate_model_performance}
# Set seed for reproducibility
set.seed(123)
# Model performance using cross-validation
perf.diablo.tcga = perf(diablo.tcga, validation = 'Mfold', folds = 10, nrepeat = 10)
```

## Visualize error rates based on weighted vote.

```{r plot_error_rates}
plot(perf.diablo.tcga) # Plot of the error rates
```

## Print the optimal number of components chosen.

```{r optimal_ncomp}
perf.diablo.tcga$choice.ncomp$WeightedVote # output the weighted vote error
```

## Automatically select the optimal number of components based on cross-validation results.

```{r select_ncomp}
# Select the optimal number of components based on the performance assessment
ncomp <- perf.diablo.tcga$choice.ncomp$WeightedVote["Overall.BER", "centroids.dist"]
ncomp # print the optimal number of components
```

## Tune the number of variables to select per omic view.

```{r tune_keepX}
# Define the range of variables to test
set.seed(123)
test.keepX <- list(mrna = c(5:9, seq(10, 25, 5)),
                   mirna = c(5:9, seq(10, 20, 2)),
                   protein = c(seq(5, 25, 5)))

# Tune number of variables for each data type
tune.diablo.tcga <- tune.block.splsda(X, Y, ncomp = 2, 
                              test.keepX = test.keepX, design = design,
                              validation = 'Mfold', folds = 10, nrepeat = 1, 
                              BPPARAM = BiocParallel::SnowParam(workers = 2),
                              dist = "centroids.dist")
list.keepX <- tune.diablo.tcga$choice.keepX # Select optimal features
list.keepX # Print the optimal features
```

## Hard code the chosen ncomp and keepX parameters

```{r custom_keepX}
# Define custom keepX parameters if you don't want to use the tuned ones
list.keepX <- list( mrna = c(8, 25), mirna = c(14,5), protein = c(10, 5))
list.keepX # Print the customized keepX values
```

## Build the final DIABLO model using the optimal parameters.

```{r final_diablo_model}
# Build the final model using the optimal number of components and features selected
diablo.tcga <- block.splsda(X, Y, ncomp = ncomp, 
                            keepX = list.keepX, design = design)
```

## Output the loadings from the final DIABLO model into CSV files.

```{r export_loadings}
# Output the loadings for each omic layer
write.csv(diablo.tcga$loadings$mRNA,"loadings-mRNA")
write.csv(diablo.tcga$loadings$miRNA,"loadings-miRNA")
write.csv(diablo.tcga$loadings$protein,"loadings-protein")
write.csv(diablo.tcga$loadings$Y,"loadings-Y")
```

## Print the design matrix for the DIABLO model.

```{r print_design_matrix}
diablo.tcga$design # Print the design matrix
```

## Extract and display the selected variables for each component.

```{r select_vars_by_comp}
# Select and display variables of the mrna data for component 2
selectVar(diablo.tcga, block = 'mrna', comp = 2)
```

# Model Visualization

## Visualize sample distributions to understand correlations and subtype discrimination.

```{r sample_plots}
# Plot samples using first two components
plotDiablo(diablo.tcga, ncomp = 1) # sample plot

# Plot samples using first two components
plotIndiv(diablo.tcga, ind.names = FALSE, legend = TRUE, 
          title = 'TCGA, DIABLO comp 1 - 2') # Plot sample distribution with subtypes

# Plot arrows showing sample position in each view with sample centroid
plotArrow(diablo.tcga, ind.names = FALSE, legend = TRUE, 
          title = 'TCGA, DIABLO comp 1 - 2') # Plot arrows showing data set correlation
```

## Visualize variable contributions and correlations through correlation circles.

```{r variable_plots}
# Plot variables using first two components
plotVar(diablo.tcga, var.names = FALSE, style = 'graphics', legend = TRUE, 
        pch = c(16, 17, 15), cex = c(2,2,2), 
        col = c('darkorchid', 'brown1', 'lightgreen'),
        title = 'TCGA, DIABLO comp 1 - 2')# plot variables

# Circular correlation plot
circosPlot(diablo.tcga, cutoff = 0.7, line = TRUE, 
           color.blocks = c('darkorchid', 'brown1', 'lightgreen'),
           color.cor = c("chocolate3","grey20"), size.labels = 1.5) # correlation circle plot

# Variable network plot
network(diablo.tcga, blocks = c(1,2,3), 
        cutoff = 0.4,
        color.node = c('darkorchid', 'brown1', 'lightgreen'),
        save ='png', name.save = 'diablo-network.png'
) # network plot
plotLoadings(diablo.tcga, comp = 1, contrib = 'max', method = 'median')#loading plot
```

## Generate and save a heatmap of DIABLO results.

```{r heatmap_plot}
pdf("cimDiablo_plot.pdf") # Plot the heatmap of model results and save it into a PDF file
cimDiablo(diablo.tcga, color.blocks = c('darkorchid', 'brown1', 'lightgreen'),
          comp = 1, margin = c(8, 20), legend.position = "right")
dev.off() # Close the PDF device
```

# Model Performance Evaluation and Prediction

## Model Performance Evaluation

## Evaluate the model performance with cross-validation.

```{r performance_evaluation}
set.seed(123)
perf.diablo.tcga <- perf(diablo.tcga,  validation = 'Mfold', folds = 10, 
                         nrepeat = 10, dist = 'centroids.dist')
```

# Output Error Rates

## Display the majority vote error rate and the weighted vote error rate.

```{r error_rates}
perf.diablo.tcga$MajorityVote.error.rate # output MajorityVote error rate
```

```{r weighted_error_rates}
perf.diablo.tcga$WeightedVote.error.rate # output WeightedVote error rate
```

# Calculate AUC

## Calculate the AUC for the specified block and component.

```{r calculate_auc}
auc.diablo.tcga <- auroc(diablo.tcga, roc.block = "protein", roc.comp = 1,
                   print = FALSE) # Calculate the AUC for protein block in component 1
```

# Prepare Test Data

## Prepare the test data from the loaded data set.

```{r prepare_test_data}
# Extract the testing data
data.test.tcga <- list(mrna = data_set$data.test$mrna, 
                       mirna = data_set$data.test$mirna)
```

# Predict Using the Test Data

## Use the DIABLO model to predict class labels for the test data set.

```{r predict_test_data}
predict.diablo.tcga <- predict(diablo.tcga, newdata = data.test.tcga) # predict test data set using the final diablo model
```

# Confusion Matrix

## Display the confusion matrix for the model performance.

```{r confusion_matrix}
# calculate confusion matrix
confusion.mat.tcga <- get.confusion_matrix(truth = data_set$data.test$subtype, 
                                           predicted = predict.diablo.tcga$WeightedVote$centroids.dist[,2])
confusion.mat.tcga # output the confusion matrix
```

# Balanced Error Rate

## Display the Balanced Error Rate (BER) for the model.

```{r balanced_error_rate}
get.BER(confusion.mat.tcga) # Calculate and output the BER
```
