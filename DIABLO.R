## ----load_libraries, message=FALSE, warning=FALSE---------------------------
required_packages <- c("ragg", "mixOmics", "ggplot2")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages) > 0) {
  stop(
    sprintf(
      "Missing required package(s): %s. Install them before running this script.",
      paste(missing_packages, collapse = ", ")
    ),
    call. = FALSE
  )
}

library(ragg)         # Graphics package
library(mixOmics)    # For multivariate analysis, including DIABLO
library(tools)        # For file path manipulation
library(ggplot2)     # For creating plots


## ----load_data, message=FALSE, warning=FALSE--------------------------------
resolve_data_root <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  data_root_arg <- args[grepl("^--data-root=", args)]
  if (length(data_root_arg) > 0) {
    return(sub("^--data-root=", "", data_root_arg[[1]]))
  }
  return(file.path(getwd(), "Dataset", "TCGA"))
}

extract_subtype <- function(folder_data, folder_name) {
  subtype <- NULL
  subtype_file <- file.path(folder_name, "subtype.csv")

  if ("subtype" %in% names(folder_data)) {
    subtype_df <- as.data.frame(folder_data$subtype)
    label_col <- intersect(c("Label", "subtype", "x"), colnames(subtype_df))
    if (length(label_col) == 0) {
      stop(sprintf("'%s/subtype.csv' exists but has no supported label column (Label/subtype/x).", folder_name), call. = FALSE)
    }
    subtype <- factor(subtype_df[[label_col[[1]]]])
    folder_data$subtype <- NULL
  } else {
    for (nm in names(folder_data)) {
      block_df <- as.data.frame(folder_data[[nm]])
      if ("Label" %in% colnames(block_df)) {
        subtype <- factor(block_df$Label)
        folder_data[[nm]] <- as.matrix(block_df[, setdiff(colnames(block_df), "Label"), drop = FALSE])
        break
      }
    }
  }

  if (is.null(subtype)) {
    stop(
      sprintf(
        "Could not resolve subtype labels for '%s'. Add '%s' with one of [Label, subtype, x] columns or include a 'Label' column in one block file.",
        folder_name,
        subtype_file
      ),
      call. = FALSE
    )
  }

  folder_data$subtype <- subtype
  folder_data
}

data_root <- resolve_data_root()
if (!dir.exists(data_root)) {
  stop(sprintf("Data root not found: %s", data_root), call. = FALSE)
}
required_folders <- c("data.train", "data.test")
missing_folders <- required_folders[!dir.exists(file.path(data_root, required_folders))]
if (length(missing_folders) > 0) {
  stop(
    sprintf(
      "Missing required folder(s) under '%s': %s",
      data_root,
      paste(missing_folders, collapse = ", ")
    ),
    call. = FALSE
  )
}

# List directories in the working directory
folders <- list.dirs(path = data_root, full.names = TRUE, recursive = FALSE)
data_set <- list() # Initialize an empty list to store data

# Loop through each folder (representing a data type)
for (folder in folders) {
  folder_name <- basename(folder)  # Extract folder name
  folder_files <- list.files(path = folder, pattern = "\\.csv$", full.names = TRUE)  # List all CSV files in the folder
  folder_data <- list()  # Initialize list to hold data from each file
  # Loop through each CSV file
  for (file in folder_files) {
    df <- read.csv(file, header = TRUE, row.names = 1) # Read csv data into dataframe
    name <- tools::file_path_sans_ext(basename(file)) # Extract file name without extension

    row_names <- rownames(df)
    df_matrix <- as.matrix(df)
    rownames(df_matrix) <- row_names # set row names for matrix
    folder_data[[name]] <- df_matrix # add matrix to folder_data list
  }
  folder_data <- extract_subtype(folder_data, folder_name)
  data_set[[folder_name]] <- folder_data # add folder data to main data_set list
}

names(data_set) # Print the names of the data_set list

required_blocks <- c("mrna", "mirna", "protein")
missing_train_blocks <- setdiff(required_blocks, names(data_set$data.train))
missing_test_blocks <- setdiff(required_blocks, names(data_set$data.test))
if (length(missing_train_blocks) > 0) {
  stop(sprintf("Missing training block(s): %s", paste(missing_train_blocks, collapse = ", ")), call. = FALSE)
}
if (length(missing_test_blocks) > 0) {
  stop(sprintf("Missing testing block(s): %s", paste(missing_test_blocks, collapse = ", ")), call. = FALSE)
}


## ----extract_training_data--------------------------------------------------
# Determine the number of omic views
number_of_view <- length(data_set$data.train) - 1
X <- list() # Initialize a list to store the omic data
# Loop through each view and store it in list
for (i in 1:number_of_view) {
  view_name <- names(data_set$data.train[i])
  X[[view_name]] <- data_set$data.train[[view_name]]
}

names(X) # Print the names of the views


## ----extract_subtype_data---------------------------------------------------
Y <- data_set$data.train$subtype # Get subtype factor
summary(Y) # Provide a summary of the subtype factor


## ----define_design_matrix---------------------------------------------------
# Set design matrix for DIABLO model, this supervises the correlation. A higher value leads to stronger supervision.
design <- matrix(0.1, ncol = length(X), nrow = length(X),
                 dimnames = list(names(X), names(X)))
# Set diagonal to zero as self-connections are not considered
diag(design) <- 0
design # Print the design matrix


## ----unsupervised_pls-------------------------------------------------------
# Unsupervised PLS between mRNA and protein
res1.pls.tcga <- pls(X$mrna, X$protein, ncomp = 1)
cor(res1.pls.tcga$variates$X, res1.pls.tcga$variates$Y) # Check correlation of components

# Unsupervised PLS between mRNA and miRNA
res2.pls.tcga <- pls(X$mrna, X$mirna, ncomp = 1)
cor(res2.pls.tcga$variates$X, res2.pls.tcga$variates$Y) # Check correlation of components

# Unsupervised PLS between miRNA and protein
res3.pls.tcga <- pls(X$mirna, X$protein, ncomp = 1)
cor(res3.pls.tcga$variates$X, res3.pls.tcga$variates$Y) # Check correlation of components


## ----build_diablo_model-----------------------------------------------------
# Create a block.plsda model
diablo.tcga <- block.plsda(X, Y, ncomp = 5, design = design)
plot(diablo.tcga) # plot of the diablo model


## ----evaluate_model_performance---------------------------------------------
# Set seed for reproducibility
set.seed(123)
# Model performance using cross-validation
perf.diablo.tcga = perf(diablo.tcga, validation = 'Mfold', folds = 10, nrepeat = 10)


## ----plot_error_rates-------------------------------------------------------
plot(perf.diablo.tcga) # Plot of the error rates


## ----optimal_ncomp----------------------------------------------------------
perf.diablo.tcga$choice.ncomp$WeightedVote # output the weighted vote error


## ----select_ncomp-----------------------------------------------------------
# Select the optimal number of components based on the performance assessment
ncomp <- perf.diablo.tcga$choice.ncomp$WeightedVote["Overall.BER", "centroids.dist"]
ncomp # print the optimal number of components


## ----tune_keepX-------------------------------------------------------------
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


## ----custom_keepX-----------------------------------------------------------
# Define custom keepX parameters if you don't want to use the tuned ones
list.keepX <- list( mrna = c(8, 25), mirna = c(14,5), protein = c(10, 5))
list.keepX # Print the customized keepX values


## ----final_diablo_model-----------------------------------------------------
# Build the final model using the optimal number of components and features selected
diablo.tcga <- block.splsda(X, Y, ncomp = ncomp, 
                            keepX = list.keepX, design = design)


## ----export_loadings--------------------------------------------------------
# Output the loadings for each omic layer
for (block_name in names(diablo.tcga$loadings)) {
  write.csv(diablo.tcga$loadings[[block_name]], paste0("loadings-", block_name, ".csv"))
}


## ----print_design_matrix----------------------------------------------------
diablo.tcga$design # Print the design matrix


## ----select_vars_by_comp----------------------------------------------------
# Select and display variables of the mrna data for component 2
selectVar(diablo.tcga, block = 'mrna', comp = 2)


## ----sample_plots-----------------------------------------------------------
# Plot samples using first two components
plotDiablo(diablo.tcga, ncomp = 1) # sample plot

# Plot samples using first two components
plotIndiv(diablo.tcga, ind.names = FALSE, legend = TRUE, 
          title = 'TCGA, DIABLO comp 1 - 2') # Plot sample distribution with subtypes

# Plot arrows showing sample position in each view with sample centroid
plotArrow(diablo.tcga, ind.names = FALSE, legend = TRUE, 
          title = 'TCGA, DIABLO comp 1 - 2') # Plot arrows showing data set correlation


## ----variable_plots---------------------------------------------------------
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


## ----heatmap_plot-----------------------------------------------------------
pdf("cimDiablo_plot.pdf") # Plot the heatmap of model results and save it into a PDF file
cimDiablo(diablo.tcga, color.blocks = c('darkorchid', 'brown1', 'lightgreen'),
          comp = 1, margin = c(8, 20), legend.position = "right")
dev.off() # Close the PDF device


## ----performance_evaluation-------------------------------------------------
set.seed(123)
perf.diablo.tcga <- perf(diablo.tcga,  validation = 'Mfold', folds = 10, 
                         nrepeat = 10, dist = 'centroids.dist')


## ----error_rates------------------------------------------------------------
perf.diablo.tcga$MajorityVote.error.rate # output MajorityVote error rate


## ----weighted_error_rates---------------------------------------------------
perf.diablo.tcga$WeightedVote.error.rate # output WeightedVote error rate


## ----calculate_auc----------------------------------------------------------
auc.diablo.tcga <- auroc(diablo.tcga, roc.block = "protein", roc.comp = 1,
                   print = FALSE) # Calculate the AUC for protein block in component 1


## ----prepare_test_data------------------------------------------------------
# Extract the testing data
data.test.tcga <- list(mrna = data_set$data.test$mrna, 
                       mirna = data_set$data.test$mirna,
                       protein = data_set$data.test$protein)


## ----predict_test_data------------------------------------------------------
expected_prediction_blocks <- names(diablo.tcga$X)
provided_prediction_blocks <- names(data.test.tcga)
missing_prediction_blocks <- setdiff(expected_prediction_blocks, provided_prediction_blocks)
unexpected_prediction_blocks <- setdiff(provided_prediction_blocks, expected_prediction_blocks)

if (length(missing_prediction_blocks) > 0 || length(unexpected_prediction_blocks) > 0) {
  stop(
    sprintf(
      paste0(
        "Prediction block mismatch. Expected blocks: [%s]. Provided blocks: [%s]. ",
        "Missing: [%s]. Unexpected: [%s]."
      ),
      paste(expected_prediction_blocks, collapse = ", "),
      paste(provided_prediction_blocks, collapse = ", "),
      paste(missing_prediction_blocks, collapse = ", "),
      paste(unexpected_prediction_blocks, collapse = ", ")
    ),
    call. = FALSE
  )
}

predict.diablo.tcga <- predict(diablo.tcga, newdata = data.test.tcga) # predict test data set using the final diablo model


## ----confusion_matrix-------------------------------------------------------
# calculate confusion matrix
confusion.mat.tcga <- get.confusion_matrix(truth = data_set$data.test$subtype, 
                                           predicted = predict.diablo.tcga$WeightedVote$centroids.dist[,2])
confusion.mat.tcga # output the confusion matrix


## ----balanced_error_rate----------------------------------------------------
get.BER(confusion.mat.tcga) # Calculate and output the BER
