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
library(tidyverse); library(TCGAbiolinks); library(sqldf); library(UCSCXenaTools)
library(ExpressionNormalizationWorkflow); library(DESeq2)

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

install.packages("stats")
library("stats")

# Function to prepare logistic regression data (combine gene expression data and target)
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

for (i in seq_along(train_data_luad)) {
train_data_luad[[i]]$target <- ifelse(train_data_luad[[i]]$target == "Tumor", 1, 0)}

for (i in seq_along(train_data_lusc)) {
train_data_lusc[[i]]$target <- ifelse(train_data_lusc[[i]]$target == "Tumor", 1, 0)}

for (i in seq_along(test_data_luad)) {
test_data_luad[[i]]$target <- ifelse(test_data_luad[[i]]$target == "Tumor", 1, 0)}

for (i in seq_along(test_data_lusc)) {
test_data_lusc[[i]]$target <- ifelse(test_data_lusc[[i]]$target == "Tumor", 1, 0)}

logit_models_luad <- list()
for (i in seq_along(train_data_luad)) {
  logit_models_luad[[i]] <- glm(target ~ ., data = train_data_luad[[i]], family = 'binomial')}

logit_models_lusc <- list()
for (i in seq_along(train_data_lusc)) {
  logit_models_lusc[[i]] <- glm(target ~ ., data = train_data_lusc[[i]], family = 'binomial')}

# Initialize an empty list to store accuracy results
accuracy_results <- data.frame(Model = character(),
                               Accuracy = numeric(),
                               stringsAsFactors = FALSE)

# Loop through the models and calculate accuracy for each one
for (i in seq_along(logit_models_luad)) {
  
  # Get the current model
  model <- logit_models_luad[[i]]
  
  # Get the corresponding test data (assuming target column is in the test data)
  test_data <- test_data_luad[[i]]
  
  # Make predictions (predicting on the test data)
  predictions <- predict(model, newdata = test_data, type = "response")
  
  # Convert predictions to binary (0 or 1) based on a threshold (typically 0.5)
  predicted_class <- ifelse(predictions > 0.5, 1, 0)
  
  # Calculate accuracy: compare predicted class to actual target in test data
  actual_class <- test_data$target  # Assuming target is the true class in the test data
  
  accuracy <- mean(predicted_class == actual_class)
  
  # Append the results to the accuracy results table
  accuracy_results <- rbind(accuracy_results, data.frame(Model = paste("Cluster", i),
                                                        Accuracy = accuracy))
}

# Print the summary table of accuracy results
print(accuracy_results)

accuracy_results <- accuracy_results[order(accuracy_results$Accuracy, decreasing = TRUE), ] %>% filter(Accuracy >= 0.90)
accuracy_results

clusters <- as.list(names(tcga.lusc))
model <- as.list(c(8,18,32,19,14))
for (i in model){
    print(clusters[[i]])}

# Initialize an empty list to store accuracy results
accuracy_results <- data.frame(Model = character(),
                               Accuracy = numeric(),
                               stringsAsFactors = FALSE)

# Loop through the models and calculate accuracy for each one
for (i in seq_along(logit_models_lusc)) {
  
  # Get the current model
  model <- logit_models_lusc[[i]]
  
  # Get the corresponding test data (assuming target column is in the test data)
  test_data <- test_data_lusc[[i]]
  
  # Make predictions (predicting on the test data)
  predictions <- predict(model, newdata = test_data, type = "response")
  
  # Convert predictions to binary (0 or 1) based on a threshold (typically 0.5)
  predicted_class <- ifelse(predictions > 0.5, 1, 0)
  
  # Calculate accuracy: compare predicted class to actual target in test data
  actual_class <- test_data$target  # Assuming target is the true class in the test data
  
  accuracy <- mean(predicted_class == actual_class)
  
  # Append the results to the accuracy results table
  accuracy_results <- rbind(accuracy_results, data.frame(Model = paste("Cluster", i),
                                                        Accuracy = accuracy))
}

# Print the summary table of accuracy results
print(accuracy_results)

accuracy_results <- accuracy_results[order(accuracy_results$Accuracy, decreasing = TRUE), ] %>% filter(Accuracy >= 0.90)
accuracy_results

clusters <- as.list(names(tcga.lusc))
model <- as.list(c(24,6,7,30,19,31,15,1,18,22,4,11))
for (i in model){
    print(clusters[[i]])}

kfold <- createFolds(labels, k = 10)

install.packages("randomForest")
library("randomForest")

# Create a list to store the random forest models
rf_models <- list()

# Assuming the target variable is 'target' and it is present in each training dataframe.
# Loop through the list of training dataframes
for (i in seq_along(train_data_luad)) {
  
  # Get the current training data frame
  train_df <- train_data_luad[[i]]
  
  # Extract the target column (assuming the target column is named 'target')
  target <- train_df$target
  
  # Remove the target column from the feature data
  features <- train_df[, !colnames(train_df) %in% "target"]
  
  # Fit the random forest model using the features and target variable
  rf_model <- randomForest(features, target, ntree = 20)  # You can adjust ntree and other parameters as needed
  
  # Store the model in the list
  rf_models[[i]] <- rf_model
}

# Create a list to store accuracy results
accuracy_results <- data.frame(Model = character(),
                               Accuracy = numeric(),
                               stringsAsFactors = FALSE)

# Assuming test_data_luad contains your test datasets
for (i in seq_along(rf_models)) {
  
  # Get the current random forest model
  rf_model <- rf_models[[i]]
  
  # Get the corresponding test data (assuming test data has the same structure as training data)
  test_df <- test_data_luad[[i]]
  
  # Extract the target column (assuming the target column is named 'target')
  actual_target <- test_df$target
  
  # Remove the target column from the test data to get the features
  test_features <- test_df[, !colnames(test_df) %in% "target"]
  
  # Make predictions on the test data
  predictions <- predict(rf_model, test_features)
  
  # Calculate accuracy
  accuracy <- mean(predictions == actual_target)
  
  # Store the results in the accuracy results table
  accuracy_results <- rbind(accuracy_results, data.frame(Model = paste("Cluster", i),
                                                         Accuracy = accuracy))
}

# Print the summary table of accuracy results
print(accuracy_results)

rf_models[[2]]

install.packages("e1071")
library(e1071)

library(e1071)

# Create a list to store the SVM models
svm_models <- list()

# Assuming the target variable is 'target' and it is present in each training dataframe.
# Loop through the list of training dataframes
for (i in seq_along(train_data_luad)) {
  
  # Get the current training data frame
  train_df <- train_data_luad[[i]]
  
  # Extract the target column (assuming the target column is named 'target')
  target <- train_df$target
  
  # Remove the target column from the feature data
  features <- train_df[, !colnames(train_df) %in% "target"]
  
  # Train the SVM model using the features and target variable
  svm_model <- svm(features, target, kernel = "linear")  # You can adjust the kernel type as needed
  
  # Store the model in the list
  svm_models[[i]] <- svm_model
}

# Now, svm_models contains all the SVM models for each dataset

train_df <- train_data_luad[[1]]
target <- train_df$target
  
  # Remove the target column from the feature data
features <- train_df[, !colnames(train_df) %in% "target"]

features

svm_model <- svm(features, target, kernel = "linear") 

test_df <- test_data_luad[[1]]
  
  # Extract the target column (assuming the target column is named 'target')
  actual_target <- test_df$target
  
  # Remove the target column from the test data to get the features
  test_features <- test_df[, !colnames(test_df) %in% "target"]

predictions <- predict(svm_model, test_features)

accuracy <- mean(predictions == actual_target)
accuracy

# Create a list to store accuracy results
accuracy_results <- data.frame(Model = character(),
                               Accuracy = numeric(),
                               stringsAsFactors = FALSE)

# Assuming test_data_luad contains your test datasets
for (i in seq_along(svm_models)) {
  
  # Get the current SVM model
  svm_model <- svm_models[[i]]
  
  # Get the corresponding test data (assuming test data has the same structure as training data)
  test_df <- test_data_luad[[i]]
  
  # Extract the target column (assuming the target column is named 'target')
  actual_target <- test_df$target
  
  # Remove the target column from the test data to get the features
  test_features <- test_df[, !colnames(test_df) %in% "target"]
  
  # Make predictions on the test data
  predictions <- predict(svm_model, test_features)
  
  # Calculate accuracy
  accuracy <- mean(predictions == actual_target)
  
  # Store the results in the accuracy results table
  accuracy_results <- rbind(accuracy_results, data.frame(Model = paste("Cluster", i),
                                                         Accuracy = accuracy))
}

# Print the summary table of accuracy results
print(accuracy_results)

# Create a list to store the PCA results and models
pca_models <- list()

# Assuming 'train_data_luad' is a list of training data frames and each has a 'target' column
for (i in seq_along(train_data_luad)) {
  
  # Get the current training data frame
  train_df <- train_data_luad[[i]]
  
  # Extract the target variable (assuming it's named 'target')
  target <- train_df$target
  
  # Remove the target variable from the training data to keep only the features
  features <- train_df[, !colnames(train_df) %in% "target"]
  
  # Apply PCA on the features
  pca_result <- prcomp(features, center = TRUE, scale. = TRUE)
  
  # Get the principal components (scores) for each sample
  pca_scores <- pca_result$x  # These are the transformed features
  
  # Combine the PCA scores with the target variable to prepare the dataset
  pca_train_df <- cbind(pca_scores, target = target)
  
  # Now you can train a machine learning model (e.g., SVM) using the PCA scores
  # Assuming you're using an SVM model here, but it can be any machine learning algorithm
  svm_model <- svm(target ~ ., data = pca_train_df, kernel = "linear")  # Example using SVM
  
  # Store the trained model in the list
  pca_models[[i]] <- svm_model
}


pca_test_results <- list()

for (i in seq_along(test_data_luad)) {
  
  # Get the current test data
  test_df <- test_data_luad[[i]]
  
  # Extract the target variable
  actual_target <- test_df$target
  
  # Remove the target variable to get the features
  test_features <- test_df[, !colnames(test_df) %in% "target"]
  
  # Apply PCA to the test features (using the same PCA transformation as for the training set)
  pca_scores_test <- predict(pca_result, newdata = test_features)
  
  # Use the trained model to make predictions
  predictions <- predict(pca_models[[i]], newdata = pca_scores_test)
  
  # Calculate the accuracy or other evaluation metrics as needed
  accuracy <- mean(predictions == actual_target)
  
  # Store the results
  pca_test_results[[i]] <- list(accuracy = accuracy, predictions = predictions)
}

# You can print or analyze the test results
print(pca_test_results)
