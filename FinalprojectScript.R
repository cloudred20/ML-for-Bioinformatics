## Step 1: Downloading and Preparing TCGA miRNA Expression Data ----

# Create a directory for storing the downloaded data
dir.create("./data")

## 1. Load Required Packages ----

# Load essential libraries for data handling and analysis
library(tidyverse)
library(TCGAbiolinks)
library(DESeq2)

## 2. Download miRNA Expression Data ----

# Define the list of TCGA cancer projects to download
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

## 3. Extract Tumor Purity Estimates ----

# Retrieve tumor purity estimates for the queried miRNA datasets
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

## 5. Download miRNA Datasets ----

# Download the filtered miRNA datasets from the GDC API
lapply(query.miRNA.download, function(x) {
  GDCdownload(x, directory = ".")
})

# Prepare and save the downloaded datasets as summarized experiments
for (i in cancer.list) {
  GDCprepare(query.miRNA.download[[i]], 
             directory = ".", 
             save = TRUE, 
             summarizedExperiment = TRUE, 
             remove.files.prepared = TRUE)
}

# Move all files containing "TCGA" in their name to the "data" directory
file.rename(from = list.files(pattern = "TCGA"), 
            to = file.path("data", list.files(pattern = "TCGA")))

# Load .RData files from the "data/" directory into a list of environments
dat.mir <- lapply(list.files("data/", pattern = "RData", full.names = TRUE), 
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

# Count the total number of miRNAs present in each dataset before pre-filtering
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

## 4. Normalization of MiRNA expression matrix ----

# Prepare the countData 
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
  write.csv(dat.mir[[i]], file = paste0("data/dat_mir_", i, ".csv"), row.names = TRUE)
}

# Save condition (assuming condition is a list of data frames)
for (i in names(condition)) {  # Loop over names of 'condition' instead of 'dat.mir'
  write.csv(condition[[i]], file = paste0("data/condition_", i, ".csv"), row.names = FALSE)
}

## 6. Genomic Annotation of Clustered MiRNAs ----

# Download the human genome GFF3 file from miRBase containing miRNA annotations
# Use curl to download the file
system("curl -o data/hsa.gff3 https://www.mirbase.org/download/hsa.gff3")

# Read the GFF3 file into a data frame
hsa <- read.table("data/hsa.gff3")

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
write.table(hsa, "data/hsa_cluster.txt", sep = "\t")

## . Print version information about R, the OS and attached or loaded packages ----
sessionInfo()
