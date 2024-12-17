## We first create a directory for downloading the data
dir.create("./Final_project_data")R.versioninstall.packages("tidyverse")
install.packages("TCGAbiolinks")
install.packages("UCSCXenaTools")
install.packages("ExpressionNormalizationWorkflow")
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("DESeq2")install.packages("sqldf")
install.packages("dplyr")
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("TCGAbiolinks")
# library(tidyverse); library(TCGAbiolinks); library(sqldf); library(UCSCXenaTools)
# library(ExpressionNormalizationWorkflow); 

library("caret")

library(DESeq2)

cancer.list <- c("TCGA-LUAD", "TCGA-LUSC")

# Query GDC API to check for available miRNA expression data
query.miRNA <- NULL

for (i in cancer.list) {                          
  query.miRNA[[i]] <- GDCquery(
    project = i,
    data.category = "Transcriptome Profiling",
    data.type = "miRNA Expression Quantification",
    experimental.strategy = "miRNA-Seq",
    sample.type = c("Primary Tumor", "Solid Tissue Normal")
  )
}

# Summarize the number of tumor and normal samples in the miRNA datasets (before filtering)
miRNA_samples <- lapply(query.miRNA, function(x) {
  table(x[[1]][[1]]$sample_type)
})

# Convert sample summary into a data frame
miRNA_samples <- data.frame(
  names(miRNA_samples), 
  data.table::rbindlist(lapply(miRNA_samples, as.data.frame.list), fill = TRUE)
)

# Rename the first column for clarity
colnames(miRNA_samples)[1] <- "cancer.list"

# Display the sample summary
miRNA_samples

miRNA_tumor_purity <- NULL

for (i in miRNA_samples[, 1]) {
  miRNA_tumor_purity[[i]] <- TCGAtumor_purity(
    getResults(query.miRNA[[i]])$cases,
    estimate = 0,
    absolute = 0,
    lump = 0,
    ihc = 0,
    cpe = 0.6
  )
}

## 4. Filter miRNA Data by Tumor Purity ----

# Query GDC API to download only samples with available tumor purity estimates
query.miRNA.download <- NULL

for (i in cancer.list) {                          
  query.miRNA.download[[i]] <- GDCquery(
    project = i,
    data.category = "Transcriptome Profiling",
    data.type = "miRNA Expression Quantification",
    experimental.strategy = "miRNA-Seq",
    sample.type = c("Primary Tumor", "Solid Tissue Normal"),
    barcode = c(miRNA_tumor_purity[[i]]$pure_barcodes, miRNA_tumor_purity[[i]]$filtered)
  )
}

query.miRNA.download

lapply(query.miRNA.download, function(x) {
  GDCdownload(x, directory = ".")
})

GDCprepare(query.miRNA.download[['TCGA-LUAD']], 
             directory = ".", 
             save = TRUE, 
             summarizedExperiment = TRUE,
             save.filename = "TCGA-LUADTranscriptome_Profiling.RData",
             remove.files.prepared = TRUE)

GDCprepare(query.miRNA.download[['TCGA-LUSC']], 
             directory = ".", 
             save = TRUE, 
             summarizedExperiment = TRUE,
             save.filename = "TCGA-LUSCTranscriptome_Profiling.RData",
             remove.files.prepared = TRUE)

# Move all files containing "TCGA" in their name to the "data" directory
file.rename(from = list.files(pattern = "TCGA"), 
            to = file.path("Final_project_data", list.files(pattern = "TCGA")))

# Load .RData files from the "Final_project_data/" directory into a list of environments
dat.mir <- lapply(list.files("Final_project_data/", pattern = "RData", full.names = TRUE), 
                  function(x) {
                    new.env <- new.env()
                    load(x, envir = new.env)
                    new.env
                  })

# Extract the "data" object from each environment and store it in the list
dat.mir <- lapply(dat.mir, function(env) env[["data"]])

# Assign cancer types (from `cancer.list`) as names to the list
names(dat.mir) <- cancer.list

# Filter miRNA expression data to retain only miRNA IDs and count columns
for (i in cancer.list) {
  dat.mir[[i]] <- subset(dat.mir[[i]], 
                         select = c("miRNA_ID", 
                                    colnames(dat.mir[[i]])[grep("count", 
                                                                colnames(dat.mir[[i]]))]))
  colnames(dat.mir[[i]]) <- gsub("read_count_", "", colnames(dat.mir[[i]]))
}

# Count the number of tumor and normal samples in each cancer dataset (after filtering)
miRNA_samples <- lapply(query.miRNA.download, function(x) {
  table(x[[1]][[1]]$sample_type)
})

miRNA_samples <- data.frame(cancer.list, 
                            data.table::rbindlist(lapply(miRNA_samples, as.data.frame.list), fill = TRUE))

miRNA_samples

as.data.frame(lapply(dat.mir, function(x) {
  dim(x)[1]
}))

# Filter out miRNAs with low read counts (row sum < 10)
row.sum <- lapply(dat.mir, function(x) {
  # Exclude the first column (miRNA_ID) from the row sum
  which((rowSums(x[, -1]) >= 10) == TRUE) 
})

# Apply filtering to retain only miRNAs with sufficient read counts
for (i in cancer.list) {
  dat.mir[[i]] <- dat.mir[[i]][row.sum[[i]], ]
}

# Count the total number of miRNAs remaining after pre-filtering
as.data.frame(lapply(dat.mir, function(x) {
  dim(x)[1]
}))

dat.mir <- lapply(dat.mir, function(x){
  # Set the first column (miRNA_ID) as row names
  rownames(x) <- x[,1]  
  # Remove the first column from the count data
  x <- x[,-1]  
})


# Prepare colData by creating a data frame with sample names and extracting the sample type 
condition <- lapply(dat.mir, function(x){
  data.frame(colnames(x), 
             # Extract sample type from the column names (e.g., "01" or "11")
             substr(colnames(x), start = 14, stop = 15))  
})

# Clean and modify the condition data to assign appropriate sample types (Tumor/Normal)
condition <- lapply(condition, function(x){
  # Name the columns as 'cases' and 'sample_type'
  i <- setNames(object = x, nm = c("cases", "sample_type"))  
  # Label Tumor as "Tumor" and others as "Normal"
  i <- mutate(i, sample_type = ifelse(sample_type == "01", "Tumor", "Normal"))

    
})

#specify experiment
condition[[1]]$experiment<- "TCGA-LUAD"
condition[[2]]$experiment<- "TCGA-LUSC"


# Create a list of DESeqDataSet objects using the count data and sample information
dds <- NULL
for (i in 1:length(dat.mir)) {
  dds[[i]] <- DESeqDataSetFromMatrix(countData = dat.mir[[i]],  
                                     colData = condition[[i]],  
                                     design = ~sample_type)  
}

# Perform normalization to account for sequencing depth and library composition 
vsd <- NULL
for (i in 1:length(dds)) {
  vsd[[i]] <- varianceStabilizingTransformation(dds[[i]], blind=F)  
}

# Save dat.mir (assuming dat.mir is a list of data frames)
for (i in names(dat.mir)) {
  write.csv(dat.mir[[i]], file = paste0("Final_project_data/dat_mir_", i, ".csv"), row.names = TRUE)
}

# Save condition (assuming condition is a list of data frames)
for (i in names(condition)) {  # Loop over names of 'condition' instead of 'dat.mir'
  write.csv(condition[[i]], file = paste0("Final_project_data/condition_", i, ".csv"), row.names = FALSE)
}


# Download the human genome GFF3 file from miRBase containing miRNA annotations
# Use curl to download the file
system("curl -o Final_project_data/hsa.gff3 https://www.mirbase.org/download/hsa.gff3")

# Read the GFF3 file into a data frame
hsa <- read.table("Final_project_data/hsa.gff3")

# Filter the data to keep only the rows related to "miRNA_primary_transcript"
# and select relevant columns for further analysis
hsa <- hsa %>% filter(V3 == "miRNA_primary_transcript") %>%
  select(V1, V4, V5, V7, V9)

# Clean up the 'V9' column to extract the miRNA names
hsa <- hsa %>% mutate(V9 = gsub(".*Name=", "", V9))

# Rename columns for easier interpretation
colnames(hsa) <- c("Chr", "Start", "End", "Strand", "MiRNA")

# Compute the inter-miRNA distance and create flags to identify clusters based on a 10Kb distance
hsa <- hsa %>% group_by(Chr) %>%
  mutate(difference_1 = Start - lag(End, default = dplyr::first(Start)),
         difference_2 = abs(End - lead(Start, default = dplyr::first(End)))) %>%
  # Flag based on proximity to previous miRNA
  # Flag based on proximity to next miRNA
  mutate(flag_1 = ifelse(difference_1 < 10000, T, F),  
         flag_2 = ifelse(difference_2 < 10000, T, F))  

# Filter out miRNAs that do not meet the clustering criteria (distance < 10Kb)
hsa <- hsa %>% filter(!(difference_1 == 0 & flag_1 == T & flag_2 == F))

# Reorder the data by combining non-clustered miRNAs and clustered miRNAs
hsa <- rbind(hsa %>% filter(flag_1 == F & lead(flag_1) == T),
             hsa %>% filter(flag_1 == T))

# Arrange the miRNAs by chromosome and start position for easier processing
hsa <- hsa %>% arrange(Chr, Start)

# Using SQL to assign a cluster number to each miRNA, considering proximity to adjacent miRNAs
hsa <- sqldf::sqldf("SELECT hs1.*, SUM(flag)
                    OVER (PARTITION BY Chr ORDER BY Start) as cluster
                    FROM (SELECT hs1.*, ((Start - LAG(Start, 1, Start)
                    OVER (PARTITION BY Chr ORDER BY Start) > 10000)) as flag
                    FROM hsa hs1) hs1;")

# Add a unique cluster label based on the chromosome and cluster number
hsa <- hsa %>%
  mutate(cluster = cluster + 1,  # Increment cluster number
         cluster = paste(Chr, "_", cluster, sep="")) %>%
  select(c(Chr, Start, End, Strand, MiRNA, cluster))

# Save the clustered miRNA data to a text file
write.table(hsa, "Final_project_data/hsa_cluster.txt", sep = "\t")

#Label DESeq objects according to dataset
LUAD <- dds[[1]]
LUSC <- dds[[2]]

print(LUAD)

#set factorlevel
LUAD$sample_type <- relevel(LUAD$sample_type, ref = "Normal")
LUSC$sample_type <- relevel(LUSC$sample_type, ref = "Normal")

LUAD <- DESeq(LUAD)

res_luad <- results(LUAD)
res_luad

summary(res_luad)

plotMA(res_luad)

#Convert into df and order by lowest p adjusted value
LUAD_df <- as.data.frame(res_luad)
LUAD_df <- LUAD_df[order(LUAD_df$padj),]
LUAD_df

#filter based on p adjusted value
filter_LUAD <- LUAD_df %>% filter(LUAD_df$padj <0.05)
filter_LUAD <- filter_LUAD %>% filter(abs(filter_LUAD$log2FoldChange)>1)
filter_LUAD

LUSC <- DESeq(LUSC)

res_lusc <- results(LUSC)
res_lusc

summary(res_lusc)

plotMA(res_lusc)

#Convert into df and order by lowest p adjusted value
LUSC_df <- as.data.frame(res_lusc)
LUSC_df <- LUSC_df[order(LUSC_df$padj),]
LUSC_df

#filter based on p adjusted value
filter_LUSC <- LUSC_df %>% filter(LUSC_df$padj <0.05)
filter_LUSC <- filter_LUSC %>% filter(abs(filter_LUSC$log2FoldChange)>1)
filter_LUSC

## 1. Feature Selection for miRNA Clusters in LUAD and LUSC Datasets----
# Create a list where each element contains miRNAs for a specific cluster
miRNA_clusters <- hsa %>%
  group_by(cluster) %>%
  summarise(MiRNA_list = list(MiRNA)) %>%
  deframe()  

# Create lists of LUAD and LUSC datasets for each cluster
tcga.luad <- list()
tcga.lusc <- list()

for (i in seq_along(miRNA_clusters)) {
  tcga.luad[[i]] <- dat.mir[[1]] %>% filter(rownames(dat.mir[[1]]) %in% 
                                              miRNA_clusters[[i]])
  tcga.lusc[[i]] <- dat.mir[[2]] %>% filter(rownames(dat.mir[[2]]) %in% 
                                              miRNA_clusters[[i]])
}

# Assign cluster names to the lists
names(tcga.luad) <- names(miRNA_clusters)
names(tcga.lusc) <- names(miRNA_clusters)

# Filter out clusters with fewer than 4 miRNAs
tcga.luad <- tcga.luad[sapply(tcga.luad, nrow) > 2]
tcga.lusc <- tcga.lusc[sapply(tcga.lusc, nrow) > 2]

## 2. Split the data into training and test data----

# Function to create training and testing datasets based on condition data
create_train_test <- function(dat, condition_data) {
  # Set the seed for reproducibility
  set.seed(123)
  
  # Create a partition based on "sample_type"
  train_index <- createDataPartition(condition_data$sample_type, p = 0.7, list = FALSE)
  
  # Extract training and testing cases
  train_condition <- condition_data[train_index, ]
  test_condition <- condition_data[-train_index, ]
  
  # Get the cases for training and testing
  train_cases <- train_condition$cases
  test_cases <- test_condition$cases
  
  # Initialize lists for the training and testing datasets
  train_data <- list()
  test_data <- list()
  
  # Loop through each cluster and filter based on cases
  for (i in seq_along(dat)) {
    cluster_data <- dat[[i]]
    train_data[[i]] <- cluster_data[, colnames(cluster_data) %in% train_cases]
    test_data[[i]] <- cluster_data[, colnames(cluster_data) %in% test_cases]
  }
  
  return(list(train = train_data, test = test_data))
}

# Apply the function to tcga.luad (condition[[1]]) and tcga.lusc (condition[[2]])
train_test_luad <- create_train_test(tcga.luad, condition[[1]])
train_test_lusc <- create_train_test(tcga.lusc, condition[[2]])

# Access the training and testing data for tcga.luad and tcga.lusc
train_tcga_luad <- train_test_luad$train
test_tcga_luad <- train_test_luad$test
train_tcga_lusc <- train_test_lusc$train
test_tcga_lusc <- train_test_lusc$test

# Function to extract the target variable (sample_type) from condition data based on cases
extract_target <- function(condition_data, cases) {
  return(condition_data$sample_type[condition_data$cases %in% cases])
}
# Initialize an empty list to store target variables
train_luad_targets <- list()
test_luad_targets <- list()

for (i in seq_along(train_tcga_luad)) {
  # Extract the cases for the current cluster
  train_cases <- colnames(train_tcga_luad[[i]])
  test_cases <- colnames(test_tcga_luad[[i]])

  # Extract the target variables for training and testing data
  train_luad_targets[[i]] <- extract_target(condition[[1]], train_cases)
  test_luad_targets[[i]] <- extract_target(condition[[1]], test_cases)
}

# Initialize an empty list to store target variables
train_lusc_targets <- list()
test_lusc_targets <- list()

for (i in seq_along(train_tcga_lusc)) {
  # Extract the cases for the current cluster
  train_cases <- colnames(train_tcga_lusc[[i]])
  test_cases <- colnames(test_tcga_lusc[[i]])

  # Extract the target variables for training and testing data
  train_lusc_targets[[i]] <- extract_target(condition[[2]], train_cases)
  test_lusc_targets[[i]] <- extract_target(condition[[2]], test_cases)
}

prepare_logistic_regression_data <- function(gene_expression_data, target) {
  # Transpose the gene expression data to have samples as rows
  logistic_data <- t(gene_expression_data)
  logistic_data <- data.frame(logistic_data)
  
  # Add the target variable (sample_type) as a column
  logistic_data$target <- target
  
  return(logistic_data)
}

# Initialize lists to store the prepared data for logistic regression
train_data_luad <- list()
test_data_luad <- list()

# Loop through each cluster and prepare the logistic regression data
for (i in seq_along(train_tcga_luad)) {
  train_data_luad[[i]] <- prepare_logistic_regression_data(train_tcga_luad[[i]], train_luad_targets[[i]])
  test_data_luad[[i]] <- prepare_logistic_regression_data(test_tcga_luad[[i]], test_luad_targets[[i]])
}

# Initialize lists to store the prepared data for logistic regression
train_data_lusc <- list()
test_data_lusc <- list()

# Loop through each cluster and prepare the logistic regression data
for (i in seq_along(train_tcga_lusc)) {
  train_data_lusc[[i]] <- prepare_logistic_regression_data(train_tcga_lusc[[i]], train_lusc_targets[[i]])
  test_data_lusc[[i]] <- prepare_logistic_regression_data(test_tcga_lusc[[i]], test_lusc_targets[[i]])
}

# for (i in seq_along(train_data_luad)) {
# train_data_luad[[i]]$target <- ifelse(train_data_luad[[i]]$target == "Tumor", 1, 0)}

# for (i in seq_along(train_data_lusc)) {
# train_data_lusc[[i]]$target <- ifelse(train_data_lusc[[i]]$target == "Tumor", 1, 0)}

# for (i in seq_along(test_data_luad)) {
# test_data_luad[[i]]$target <- ifelse(test_data_luad[[i]]$target == "Tumor", 1, 0)}

# for (i in seq_along(test_data_lusc)) {
# test_data_lusc[[i]]$target <- ifelse(test_data_lusc[[i]]$target == "Tumor", 1, 0)}


#f1 analysis
calculate_f1_score <- function(predictions, actual) {
  cm <- confusionMatrix(predictions, actual)
  precision <- cm$byClass["Precision"]
  recall <- cm$byClass["Recall"]
  f1 <- 2 * (precision * recall) / (precision + recall)
  return(f1)
}

install.packages("kernlab")
library(kernlab)
library(e1071)  
library(pROC)

# Function to perform K-fold cross-validation for SVM and calculate accuracy and F1 score
perform_svm_cv_and_evaluate <- function(train_data_list, test_data_list,list_of_clusters, k=5) {

    
    #have list of clusters as list
    cluster <- as.list(names(list_of_clusters))
  # Initialize an empty data frame to store the results (model, accuracy, F1 score)
  results <- data.frame(Model = character(),
                        Accuracy = numeric(),
                        F1_Score = numeric(),
                        stringsAsFactors = FALSE)
  
  # Loop through each cluster in the data list
  for (i in seq_along(train_data_list)) {
    
    # Get the current training data
    train_data <- train_data_list[[i]]
    
    # Ensure target is a factor for classification
    train_data$target <- factor(train_data$target, levels = c("Normal", "Tumor"))  # 0 = Normal, 1 = Tumor
    
    # Define training control for K-fold cross-validation (k = 5)
    train_control <- trainControl(method = "cv", number = k, 
                                  summaryFunction = twoClassSummary, 
                                  classProbs = TRUE, 
                                  savePredictions = "final")
    
    # Train SVM model (SVM with linear kernel)
    svm_model <- train(target ~ ., data = train_data, method = "svmLinear", 
                       trControl = train_control, 
                       metric = "ROC")
    
    # Print the summary of results (accuracy, ROC, etc.)
    print(paste("Results for Cluster", cluster[[i]]))
    print(svm_model$results)
    
    # Get the predictions for the test data
    test_data <- test_data_list[[i]]
    
    # Ensure test data has the same preprocessing (target as factor)
    test_data$target <- factor(test_data$target, levels = c("Normal", "Tumor"))  # target as factor
    
    # Get predictions on the test data
    test_predictions <- predict(svm_model, newdata = test_data)

    # Test Probabilities
    test_probabilities <- predict(svm_model, newdata = test_data, type = "prob")

    # Calculate ROC for the positive class (assuming the second column is the probability of the Tumor class)
    roc_curve <- roc(test_data$target, test_probabilities[, 2])  # second column corresponds to Tumor class
    plot(roc_curve, main = paste("ROC Curve for Cluster", cluster[[i]]), col = "blue")
    
    # Calculate accuracy for the test data
    accuracy <- mean(test_predictions == test_data$target)
    
    # Calculate F1 score for the test data
    f1_score <- calculate_f1_score(test_predictions, test_data$target)
    
    # Store the results in the results table
    results <- rbind(results, data.frame(Model = paste("Cluster", cluster[[i]]),
                                         Accuracy = accuracy,
                                         F1_Score = f1_score))
  }
  
  return(results)
}

luad_svm <- perform_svm_cv_and_evaluate(train_data_luad, test_data_luad, tcga.luad)

luad_svm <-as.data.frame(luad_svm)
print(luad_svm)

luad_top_svm <- luad_svm[order(luad_svm$Accuracy, decreasing = TRUE), ]
luad_top_svm

lusc_svm <- perform_svm_cv_and_evaluate(train_data_lusc, test_data_lusc, tcga.lusc)

lusc_svm <- as.data.frame(lusc_svm) 
print(lusc_svm)

lusc_top_svm <- lusc_svm[order(lusc_svm$Accuracy, decreasing = TRUE), ]
lusc_top_svm

# Function to perform K-fold cross-validation for SVM with grid search for hyperparameters
perform_svm_cv_and_evaluate <- function(train_data_list, test_data_list, list_of_clusters, k=5) {

  # Define a list of clusters
  cluster <- as.list(names(list_of_clusters))

  # Initialize an empty data frame to store the results (model, accuracy, F1 score)
  results <- data.frame(Model = character(),
                        Accuracy = numeric(),
                        F1_Score = numeric(),
                        Best_C = numeric(),
                        Best_Sigma = numeric(),
                        stringsAsFactors = FALSE)

  # Loop through each cluster in the data list
  for (i in seq_along(train_data_list)) {

    # Get the current training data
    train_data <- train_data_list[[i]]

    # Ensure target is a factor for classification
    train_data$target <- factor(train_data$target, levels = c("Normal", "Tumor"))  # 0 = Normal, 1 = Tumor

    # Define training control for K-fold cross-validation (k = 5)
    train_control <- trainControl(method = "cv", number = k, 
                                  summaryFunction = twoClassSummary, 
                                  classProbs = TRUE, 
                                  savePredictions = "final")

    # Define the grid of hyperparameters to search (C and sigma for SVM with RBF kernel)
    tune_grid <- expand.grid(C = 2^(-5:5), sigma = 2^(-15:3))

    # Train SVM model using grid search
    svm_model <- train(target ~ ., data = train_data, method = "svmRadial", 
                       trControl = train_control, 
                       tuneGrid = tune_grid, 
                       metric = "ROC")

    # Print the best hyperparameters and results
    print(paste("Results for Cluster", cluster[[i]]))
    print(svm_model$results)
    print(paste("Best C:", svm_model$bestTune$C))
    print(paste("Best Sigma:", svm_model$bestTune$sigma))

    # Get the predictions for the test data
    test_data <- test_data_list[[i]]

    # Ensure test data has the same preprocessing (target as factor)
    test_data$target <- factor(test_data$target, levels = c("Normal", "Tumor"))  # target as factor

    # Get predictions on the test data
    test_predictions <- predict(svm_model, newdata = test_data)

    # Test Probabilities
    test_probabilities <- predict(svm_model, newdata = test_data, type = "prob")

    # Calculate ROC for the positive class (assuming the second column is the probability of the Tumor class)
    roc_curve <- roc(test_data$target, test_probabilities[, 2])  # second column corresponds to Tumor class
    plot(roc_curve, main = paste("ROC Curve for Cluster", cluster[[i]]), col = "blue")

    # Save ROC curve data to a CSV file
    roc_data <- data.frame(
      Threshold = roc_curve$thresholds,
      Sensitivity = roc_curve$sensitivities,
      Specificity = roc_curve$specificities
    )
    
    # Generate a file name for the CSV based on the cluster
    file_name <- paste0("ROC_Curve_Cluster_", cluster[[i]], ".csv")
    
    # Write the data to CSV
    write.csv(roc_data, file = file_name, row.names = FALSE)

    # Calculate accuracy for the test data
    accuracy <- mean(test_predictions == test_data$target)

    # Calculate F1 score for the test data
    f1_score <- calculate_f1_score(test_predictions, test_data$target)

    # Store the results in the results table
    results <- rbind(results, data.frame(Model = paste("Cluster", cluster[[i]]),
                                         Accuracy = accuracy,
                                         F1_Score = f1_score,
                                         Best_C = svm_model$bestTune$C,
                                         Best_Sigma = svm_model$bestTune$sigma))
  }

  return(results)
}


op_luad_svm <- perform_svm_cv_and_evaluate(train_data_luad, test_data_luad, tcga.luad)

op_luad_svm <- as.data.frame(op_luad_svm)
print(op_luad_svm)

op_luad_top_svm <- op_luad_svm[order(op_luad_svm$Accuracy, decreasing = TRUE), ]
op_luad_top_svm

write.csv(op_luad_top_svm, "luad_svm_accuracy.csv", row.names = FALSE)

op_lusc_svm <- perform_svm_cv_and_evaluate(train_data_lusc, test_data_lusc, tcga.lusc)

op_lusc_svm <- as.data.frame(op_lusc_svm)
print(op_lusc_svm)

op_lusc_top_svm <- op_lusc_svm[order(op_lusc_svm$Accuracy, decreasing = TRUE), ]
op_lusc_top_svm

write.csv(op_lusc_top_svm, "lusc_svm_accuracy.csv", row.names = FALSE)

library(rpart)
library(pROC)
library(caret)

# Function to perform Decision Tree model training with K-fold cross-validation and evaluation
perform_decision_tree_cv_and_evaluate <- function(train_data_list, test_data_list,list_of_clusters, k = 5) {
    
    
    # Define a list of clusters
  cluster <- as.list(names(list_of_clusters))
  
  # Initialize an empty data frame to store the results (model, accuracy, F1 score, AUC)
  results <- data.frame(Model = character(),
                        Accuracy = numeric(),
                        F1_Score = numeric(),
                        AUC = numeric(),
                        stringsAsFactors = FALSE)
  
  # Loop through each pair of train and test dataframes in the list
  for (i in seq_along(train_data_list)) {
    
    # Get the current training and test data
    train_data <- train_data_list[[i]]
    test_data <- test_data_list[[i]]
    
    # Ensure that 'target' is the name of the response variable in your dataframe
    # If the target variable has a different name, replace 'target' with the correct column name
    train_data$target <- factor(train_data$target)  # Convert target to factor for classification
    test_data$target <- factor(test_data$target)  # Convert target to factor for classification
    
    # Define training control for K-fold cross-validation
    train_control <- trainControl(method = "cv", number = k, 
                                  summaryFunction = twoClassSummary, 
                                  classProbs = TRUE, 
                                  savePredictions = "final")
    
    # Train the decision tree model using rpart
    decision_tree_model <- train(target ~ ., data = train_data, method = "rpart", 
                                 trControl = train_control, 
                                 metric = "ROC")
    
    # Print the summary of results (accuracy, ROC, etc.)
    print(paste("Results for Model", cluster[[i]]))
    print(decision_tree_model$results)
    
    # Get predictions on the test data
    test_predictions <- predict(decision_tree_model, newdata = test_data)
    test_predicted_probs <- predict(decision_tree_model, newdata = test_data, type = "prob")[, 2]  # Probabilities for the positive class
    
    # Calculate accuracy
    accuracy <- mean(test_predictions == test_data$target)
    
    # Calculate F1 score
    f1_score <- calculate_f1_score(test_predictions, test_data$target)
    
    # Calculate AUC (Area Under the Curve)
    roc_curve <- roc(test_data$target, test_predicted_probs)
    auc_value <- auc(roc_curve)
    
    # Append the results to the results table
    results <- rbind(results, data.frame(Model = paste("Model", cluster[[i]]),
                                         Accuracy = accuracy,
                                         F1_Score = f1_score,
                                         AUC = auc_value))
    
    # Plot ROC curve for the current model
    plot(roc_curve, main = paste("ROC Curve for Model", cluster[[i]]))
  }
  
  # Return the results table with accuracy, F1 score, and AUC
  return(results)
}


dt_luad <- perform_decision_tree_cv_and_evaluate(train_data_luad, test_data_luad, tcga.luad)

dt_luad <- as.data.frame(dt_luad)
print(dt_luad)

luad_top_svm <- luad_svm[order(luad_svm$Accuracy, decreasing = TRUE), ]
luad_top_svm

dt_lusc <- perform_decision_tree_cv_and_evaluate(train_data_lusc, test_data_lusc, tcga.lusc)

dt_lusc <- as.data.frame(dt_lusc)
print(dt_lusc)

dt_lusc <- dt_lusc[order(dt_lusc$Accuracy, decreasing = TRUE), ]
dt_lusc

perform_decision_tree_cv_and_evaluate <- function(train_data_list, test_data_list, list_of_clusters, k = 5) {
  
  # Define a list of clusters
  cluster <- as.list(names(list_of_clusters))
  
  # Initialize an empty data frame to store the results (model, accuracy, F1 score, AUC)
  results <- data.frame(Model = character(),
                        Accuracy = numeric(),
                        F1_Score = numeric(),
                        AUC = numeric(),
                        Best_CP = numeric(),
                        stringsAsFactors = FALSE)
  
  # Define the grid of hyperparameters for tuning (just cp)
  tune_grid <- expand.grid(
    cp = seq(0.01, 0.1, by = 0.01)  # Complexity parameter (cp) only
  )
  
  # Loop through each pair of train and test dataframes in the list
  for (i in seq_along(train_data_list)) {
    
    # Get the current training and test data
    train_data <- train_data_list[[i]]
    test_data <- test_data_list[[i]]
    
    # Ensure that 'target' is the name of the response variable in your dataframe
    train_data$target <- factor(train_data$target)  # Convert target to factor for classification
    test_data$target <- factor(test_data$target)  # Convert target to factor for classification
    
    # Define training control for K-fold cross-validation
    train_control <- trainControl(method = "cv", number = k, 
                                  summaryFunction = twoClassSummary, 
                                  classProbs = TRUE, 
                                  savePredictions = "final")
    
    # Train the decision tree model using rpart and grid search
    decision_tree_model <- train(target ~ ., data = train_data, method = "rpart", 
                                 trControl = train_control, 
                                 tuneGrid = tune_grid,            # Use the grid for hyperparameter tuning
                                 metric = "ROC")                  # Optimize for ROC metric
    
    # Print the summary of results (accuracy, ROC, etc.)
    print(paste("Results for Cluster", cluster[[i]]))
    print(decision_tree_model$results)
    
    # Get the best hyperparameter (cp)
    best_cp <- decision_tree_model$bestTune$cp
    
    # Get predictions on the test data
    test_predictions <- predict(decision_tree_model, newdata = test_data)
    test_predicted_probs <- predict(decision_tree_model, newdata = test_data, type = "prob")[, 2]  # Probabilities for the positive class
    
    # Calculate accuracy
    accuracy <- mean(test_predictions == test_data$target)
    
    # Calculate F1 score
    f1_score <- calculate_f1_score(test_predictions, test_data$target)
    
    # Calculate AUC (Area Under the Curve)
    roc_curve <- roc(test_data$target, test_predicted_probs)
    auc_value <- auc(roc_curve)
    
    # Append the results to the results table
    results <- rbind(results, data.frame(Model = paste("Cluster", cluster[[i]]),
                                         Accuracy = accuracy,
                                         F1_Score = f1_score,
                                         AUC = auc_value,
                                         Best_CP = best_cp))
    
    # Plot ROC curve for the current model
    plot(roc_curve, main = paste("ROC Curve for Model", cluster[[i]]))
  }
  
  # Return the results table with accuracy, F1 score, and AUC
  return(results)
}


dt_op_luad <- perform_decision_tree_cv_and_evaluate(train_data_luad, test_data_luad, tcga.luad)

dt_op_luad <- as.data.frame(dt_op_luad)
print(dt_op_luad)

dt_op_luad <- dt_op_luad[order(dt_op_luad$Accuracy, decreasing = TRUE), ]
dt_op_luad

write.csv(dt_op_luad, "luad_dt_accuracy.csv", row.names = FALSE)

dt_op_lusc <- perform_decision_tree_cv_and_evaluate(train_data_lusc, test_data_lusc, tcga.lusc)

dt_op_lusc <- as.data.frame(dt_op_lusc)
print(dt_op_lusc)

dt_op_lusc <- dt_op_lusc[order(dt_op_lusc$Accuracy, decreasing = TRUE), ]
dt_op_lusc

write.csv(dt_op_lusc, "lusc_dt_accuracy.csv", row.names = FALSE)



install.packages("stats")
library("stats")

install.packages("caTools")
library(caTools)

perform_logistic_regression_cv_and_evaluate <- function(train_data_list, test_data_list, list_of_clusters, k = 5) {
  
  # Define a list of clusters
  cluster <- as.list(names(list_of_clusters))
  
  # Initialize an empty data frame to store the results (model, accuracy, F1 score, AUC)
  results <- data.frame(Model = character(),
                        Accuracy = numeric(),
                        F1_Score = numeric(),
                        AUC = numeric(),
                        stringsAsFactors = FALSE)
  
  # Loop through each pair of train and test dataframes in the list
  for (i in seq_along(train_data_list)) {
    
    # Get the current training and test data
    train_data <- train_data_list[[i]]
    test_data <- test_data_list[[i]]
    
    # Ensure that 'target' is the name of the response variable in your dataframe
    # If the target variable has a different name, replace 'target' with the correct column name
    train_data$target <- factor(train_data$target)  # Convert target to factor for classification
    test_data$target <- factor(test_data$target)  # Convert target to factor for classification
    
    # Define training control for K-fold cross-validation
    train_control <- trainControl(method = "cv", number = k, 
                                  summaryFunction = twoClassSummary, 
                                  classProbs = TRUE, 
                                  savePredictions = "final")
    
    # Train the logistic regression model using glm (logistic regression)
    logistic_model <- train(target ~ ., data = train_data, method = "glm", 
                            family = "binomial",  # Logistic regression for binary classification
                            trControl = train_control, 
                            metric = "ROC")  # Optimize the model based on ROC (AUC)
    
    # Print the summary of results (accuracy, ROC, etc.)
    print(paste("Results for Cluster", cluster[[i]]))
    print(logistic_model$results)
    
    # Get predictions on the test data
    test_predictions <- predict(logistic_model, newdata = test_data)
    test_predicted_probs <- predict(logistic_model, newdata = test_data, type = "prob")[, 2]  # Probabilities for the positive class
    
    # Calculate accuracy
    accuracy <- mean(test_predictions == test_data$target)
    
    # Calculate F1 score
    f1_score <- calculate_f1_score(test_predictions, test_data$target)
    
    # Calculate AUC (Area Under the Curve)
    roc_curve <- roc(test_data$target, test_predicted_probs)
    auc_value <- auc(roc_curve)
    
    # Append the results to the results table
    results <- rbind(results, data.frame(Model = paste("Cluster", cluster[[i]]),
                                         Accuracy = accuracy,
                                         F1_Score = f1_score,
                                         AUC = auc_value))
    
    # Plot ROC curve for the current model
    plot(roc_curve, main = paste("ROC Curve for Cluster", cluster[[i]]))
  }
  
  # Return the results table with accuracy, F1 score, and AUC
  return(results)
}

perform_logistic_regression_cv_and_evaluate <- function(train_data_list, test_data_list, list_of_clusters, k = 5) {
  
  # Define a list of clusters
  cluster <- as.list(names(list_of_clusters))
  
  # Initialize an empty data frame to store the results (model, accuracy, F1 score, AUC)
  results <- data.frame(Model = character(),
                        Accuracy = numeric(),
                        F1_Score = numeric(),
                        AUC = numeric(),
                        stringsAsFactors = FALSE)
  
  # Loop through each pair of train and test dataframes in the list
  for (i in seq_along(train_data_list)) {
    
    # Get the current training and test data
    train_data <- train_data_list[[i]]
    test_data <- test_data_list[[i]]
    
    # Ensure that 'target' is the name of the response variable in your dataframe
    # If the target variable has a different name, replace 'target' with the correct column name
    train_data$target <- factor(train_data$target)  # Convert target to factor for classification
    test_data$target <- factor(test_data$target)  # Convert target to factor for classification
    
    # Define training control for K-fold cross-validation
    train_control <- trainControl(method = "cv", number = k, 
                                  summaryFunction = twoClassSummary, 
                                  classProbs = TRUE, 
                                  savePredictions = "final")
    
    # Train the logistic regression model using glm (logistic regression)
    logistic_model <- train(target ~ ., data = train_data, method = "glm", 
                            family = "binomial",  # Logistic regression for binary classification
                            trControl = train_control, 
                            metric = "ROC")  # Optimize the model based on ROC (AUC)
    
    # Print the summary of results (accuracy, ROC, etc.)
    print(paste("Results for Cluster", cluster[[i]]))
    print(logistic_model$results)
    
    # Get predictions on the test data
    test_predictions <- predict(logistic_model, newdata = test_data)
    test_predicted_probs <- predict(logistic_model, newdata = test_data, type = "prob")[, 2]  # Probabilities for the positive class
    
    # Calculate accuracy
    accuracy <- mean(test_predictions == test_data$target)
    
    # Calculate F1 score
    f1_score <- calculate_f1_score(test_predictions, test_data$target)
    
    # Calculate AUC (Area Under the Curve)
    roc_curve <- roc(test_data$target, test_predicted_probs)
    auc_value <- auc(roc_curve)
    
    # Append the results to the results table
    results <- rbind(results, data.frame(Model = paste("Cluster", cluster[[i]]),
                                         Accuracy = accuracy,
                                         F1_Score = f1_score,
                                         AUC = auc_value))
    
    # Plot ROC curve for the current model
    plot(roc_curve, main = paste("ROC Curve for Cluster", cluster[[i]]))
  }
  
  # Return the results table with accuracy, F1 score, and AUC
  return(results)
}

perform_logistic_regression_cv_and_evaluate <- function(train_data_list, test_data_list, list_of_clusters, k = 5) {
  
  # Define a list of clusters
  cluster <- as.list(names(list_of_clusters))
  
  # Initialize an empty data frame to store the results (model, accuracy, F1 score, AUC)
  results <- data.frame(Model = character(),
                        Accuracy = numeric(),
                        F1_Score = numeric(),
                        AUC = numeric(),
                        stringsAsFactors = FALSE)
  
  # Loop through each pair of train and test dataframes in the list
  for (i in seq_along(train_data_list)) {
    
    # Get the current training and test data
    train_data <- train_data_list[[i]]
    test_data <- test_data_list[[i]]
    
    # Ensure that 'target' is the name of the response variable in your dataframe
    # If the target variable has a different name, replace 'target' with the correct column name
    train_data$target <- factor(train_data$target)  # Convert target to factor for classification
    test_data$target <- factor(test_data$target)  # Convert target to factor for classification
    
    # Define training control for K-fold cross-validation
    train_control <- trainControl(method = "cv", number = k, 
                                  summaryFunction = twoClassSummary, 
                                  classProbs = TRUE, 
                                  savePredictions = "final")
    
    # Train the logistic regression model using glm (logistic regression)
    logistic_model <- train(target ~ ., data = train_data, method = "glm", 
                            family = "binomial",  # Logistic regression for binary classification
                            trControl = train_control, 
                            metric = "ROC")  # Optimize the model based on ROC (AUC)
    
    # Print the summary of results (accuracy, ROC, etc.)
    print(paste("Results for Cluster", cluster[[i]]))
    print(logistic_model$results)
    
    # Get predictions on the test data
    test_predictions <- predict(logistic_model, newdata = test_data)
    test_predicted_probs <- predict(logistic_model, newdata = test_data, type = "prob")[, 2]  # Probabilities for the positive class
    
    # Calculate accuracy
    accuracy <- mean(test_predictions == test_data$target)
    
    # Calculate F1 score
    f1_score <- calculate_f1_score(test_predictions, test_data$target)
    
    # Calculate AUC (Area Under the Curve)
    roc_curve <- roc(test_data$target, test_predicted_probs)
    auc_value <- auc(roc_curve)
    
    # Append the results to the results table
    results <- rbind(results, data.frame(Model = paste("Cluster", cluster[[i]]),
                                         Accuracy = accuracy,
                                         F1_Score = f1_score,
                                         AUC = auc_value))
    
    # Plot ROC curve for the current model
    plot(roc_curve, main = paste("ROC Curve for Cluster", cluster[[i]]))
  }
  
  # Return the results table with accuracy, F1 score, and AUC
  return(results)
}

log_luad<-perform_logistic_regression_cv_and_evaluate(train_data_luad, test_data_luad, tcga.luad)

log_luad <- as.data.frame(log_luad)
print(log_luad)

log_luad <- log_luad[order(log_luad$Accuracy, decreasing = TRUE), ]
log_luad

log_lusc<-perform_logistic_regression_cv_and_evaluate(train_data_lusc, test_data_lusc, tcga.lusc)

log_lusc <- as.data.frame(log_lusc)
print(log_lusc)

log_lusc <- log_lusc[order(log_lusc$Accuracy, decreasing = TRUE), ]
log_lusc

install.packages("glmnet")
library("glmnet")

install.packages("stepPlr")
library("stepPlr")

library("caTools")



perform_logit_boost_cv_and_evaluate <- function(train_data_list, test_data_list, list_of_clusters, k = 5) {
  
  # Define a list of clusters
  cluster <- as.list(names(list_of_clusters))
  
  # Initialize an empty data frame to store the results (model, accuracy, F1 score, AUC)
  results <- data.frame(Model = character(),
                        Accuracy = numeric(),
                        F1_Score = numeric(),
                        AUC = numeric(),
                        Best_CP = numeric(),
                        stringsAsFactors = FALSE)
  
  # Define the grid of hyperparameters for tuning (just cp)
  tune_grid <- expand.grid(
    nIter = c(10,50,100)  # Complexity parameter (cp) only
  )
  
  # Loop through each pair of train and test dataframes in the list
  for (i in seq_along(train_data_list)) {
    
    # Get the current training and test data
    train_data <- train_data_list[[i]]
    test_data <- test_data_list[[i]]
    
    # Ensure that 'target' is the name of the response variable in your dataframe
    train_data$target <- factor(train_data$target)  # Convert target to factor for classification
    test_data$target <- factor(test_data$target)  # Convert target to factor for classification
    
    # Define training control for K-fold cross-validation
    train_control <- trainControl(method = "cv", number = k, 
                                  summaryFunction = twoClassSummary, 
                                  classProbs = TRUE, 
                                  savePredictions = "final")
    
    # Train the decision tree model using rpart and grid search
    decision_tree_model <- train(target ~ ., data = train_data, method = "LogitBoost", 
                                 trControl = train_control, 
                                 tuneGrid = tune_grid,            # Use the grid for hyperparameter tuning
                                 metric = "ROC")                  # Optimize for ROC metric
    
    # Print the summary of results (accuracy, ROC, etc.)
    print(paste("Results for Cluster", cluster[[i]]))
    print(decision_tree_model$results)
    
    # Get the best hyperparameter (cp)
    best_iter <- decision_tree_model$bestTune$nIter
    
    # Get predictions on the test data
    test_predictions <- predict(decision_tree_model, newdata = test_data)
    test_predicted_probs <- predict(decision_tree_model, newdata = test_data, type = "prob")[, 2]  # Probabilities for the positive class
    
    # Calculate accuracy
    accuracy <- mean(test_predictions == test_data$target)
    
    # Calculate F1 score
    f1_score <- calculate_f1_score(test_predictions, test_data$target)
    
    # Calculate AUC (Area Under the Curve)
    roc_curve <- roc(test_data$target, test_predicted_probs)
    auc_value <- auc(roc_curve)
    
    # Append the results to the results table
    results <- rbind(results, data.frame(Model = paste("Cluster", cluster[[i]]),
                                         Accuracy = accuracy,
                                         F1_Score = f1_score,
                                         AUC = auc_value,
                                         n_iter = best_iter))
    
    # Plot ROC curve for the current model
    plot(roc_curve, main = paste("ROC Curve for Model", cluster[[i]]))
  }
  
  # Return the results table with accuracy, F1 score, and AUC
  return(results)
}


log_opt_luad<-perform_logit_boost_cv_and_evaluate(train_data_luad, test_data_luad, tcga.luad)

log_opt_luad <- as.data.frame(log_opt_luad)
print(log_opt_luad)

log_opt_luad <- log_opt_luad[order(log_opt_luad$Accuracy, decreasing = TRUE), ]
log_opt_luad

log_opt_lusc<-perform_logit_boost_cv_and_evaluate(train_data_lusc, test_data_lusc, tcga.lusc)

log_opt_lusc <- as.data.frame(log_opt_lusc)
print(log_opt_lusc)

log_opt_lusc <- log_opt_lusc[order(log_opt_lusc$Accuracy, decreasing = TRUE), ]
log_opt_lusc
