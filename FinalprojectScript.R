## Step1: Downloading and preparing TCGA miRNA expression data

## We first create a directory for downloading the data
dir.create("./data")

## This script is organized into 7 sub-sections.
## 1. Load packages ----

library(tidyverse); library(TCGAbiolinks); library(sqldf); library(UCSCXenaTools)
library(ExpressionNormalizationWorkflow); library(DESeq2)

## 2. Download MiRNA expression data ----

# list of different cancer data-sets available at TCGA
cancer.list <- getGDCprojects() %>% 
  filter(grepl("TCGA", id)) %>%
  pull(id) %>% sort(.)

# querying GDC API for cancer data-sets where miRNA expression data is available
query.miRNA <- NULL

for (i in cancer.list) {                          
  query.miRNA[[i]] <- GDCquery(project = i,
                               legacy = T,
                               data.category = "Gene expression",
                               data.type = "miRNA gene quantification",
                               experimental.strategy = "miRNA-Seq",
                               file.type = "hg19.mirbase20.mirna.quantification",
                               sample.type = c("Primary Tumor",
                                               "Solid Tissue Normal"))
}

# number of tumor and normal samples in miRNA expression data of each cancer (before filtering)
miRNA_samples <- lapply(query.miRNA, function(x){
  table(x[[1]][[1]]$sample_type)
})

miRNA_samples <- data.frame(names(miRNA_samples), 
                       data.table::rbindlist(lapply(miRNA_samples, as.data.frame.list), fill = T))

colnames(miRNA_samples)[1] <- "cancer.list"

# miRNA_samples has only 32 cancers instead of 33. 
# TCGA-LAML is a blood cancer so we've removed that from our list. 
# We're analyzing only primary tumor samples. 
# Notice some cancers do not have normal samples data.

miRNA_samples

# extracting tumor purity estimates for queried miRNA expression data
miRNA_tumor_purity <- NULL

for (i in miRNA_samples[,1]) {
  miRNA_tumor_purity[[i]] <- TCGAtumor_purity(getResults(query.miRNA[[i]])$cases,
                                              estimate = 0,
                                              absolute = 0,
                                              lump = 0,
                                              ihc = 0,
                                              cpe = 0.6)
}

# removing cancer types with missing information about tumor purity estimates from the list
names(which(lapply(sapply(miRNA_tumor_purity, "[[", "pure_barcodes"), length) > 0))

# removing cancer types with less than 30 normal control samples from the list
miRNA_samples[which(miRNA_samples[,"Solid.Tissue.Normal"] >= 30), "cancer.list"]

cancer.list <- miRNA_samples %>% 
  filter(cancer.list %in% 
           intersect(names(which(lapply(sapply(miRNA_tumor_purity, "[[", "pure_barcodes"), length) > 0)),
                     miRNA_samples[which(miRNA_samples[,"Solid.Tissue.Normal"] >= 30), "cancer.list"])) %>% 
  pull(cancer.list) 

cancer.list

# querying GDC API for cancer data-sets where
# miRNA expression data is available for at least 30 normal controls
# tumor purity estimates information is available
query.miRNA.download <- NULL

for (i in cancer.list) {                          
  query.miRNA.download[[i]] <- GDCquery(project = i,
                                        legacy = T,
                                        data.category = "Gene expression",
                                        data.type = "miRNA gene quantification",
                                        experimental.strategy = "miRNA-Seq",
                                        file.type = "hg19.mirbase20.mirna.quantification",
                                        sample.type = c("Primary Tumor",
                                                        "Solid Tissue Normal"),
                                        barcode = c(miRNA_tumor_purity[[i]]$pure_barcodes,
                                                    miRNA_tumor_purity[[i]]$filtered))
}

# downloading the queried miRNA data-sets from GDC-API 
lapply(query.miRNA.download, function(x){
  GDCdownload(x, directory = ".")
})

# reading the downloaded data and preparing it into an R object

dir.create("data/MiRNA_files")

for (i in cancer.list) {
  GDCprepare(query.miRNA.download[[i]], 
             directory = ".",
             save = T,
             save.filename = paste("data/MiRNA_files/", i,"_miRNA_hg19.rda", sep = ""),
             summarizedExperiment = T,
             remove.files.prepared = F)
}

dat.mir <- lapply(list.files("data/MiRNA_files/", pattern = "rda", full.names = T), function(x){
  new.env <- new.env()
  load(x, envir = new.env)
  new.env
})

dat.mir <- lapply(dat.mir, "[[", "data")
names(dat.mir) <- cancer.list

for (i in cancer.list) {
  dat.mir[[i]] <- subset(dat.mir[[i]],
                         select = c("miRNA_ID",colnames(dat.mir[[i]])[grep("count", colnames(dat.mir[[i]]))]))
  colnames(dat.mir[[i]]) <- gsub("read_count_","", colnames(dat.mir[[i]]))
  
}

# number of tumor and normal samples in miRNA expression data of each cancer (after filtering)
miRNA_samples <- lapply(query.miRNA.download, function(x){
  table(x[[1]][[1]]$sample_type)
})

miRNA_samples <- data.frame(cancer.list, 
                            data.table::rbindlist(lapply(miRNA_samples, as.data.frame.list), fill = T))

miRNA_samples

# no. of total miRNAs present before pre-filtering
as.data.frame(lapply(dat.mir, function(x){
  dim(x)[1]
}))

# filter miRNAs with low read count

row.sum <- lapply(dat.mir, function(x){
  which((rowSums(x[,-1]) >=10) == T) # the first column store miRNA_ID
})

for (i in cancer.list) {
  dat.mir[[i]] <- dat.mir[[i]][row.sum[[i]],]
}

# no. of total miRNAs present after pre-filtering

as.data.frame(lapply(dat.mir, function(x){
  dim(x)[1]
}))

# Table.S1: List of cancers analyzed along with number of patient samples and miRNAs post filtering

Table.S1 <- miRNA_samples %>% mutate(MiRNA.Count = 
                                       as.numeric(lapply(dat.mir, function(x){
                                         dim(x)[1]
                                       })))

colnames(Table.S1) <- c("Cancer", "Tumor(n)", "Normal(n)", "MiRNA(n)")
Table.S1

# saving miRNA expression data for later use
for (i in cancer.list) {
  write.table(dat.mir[[i]], file = paste("data/MiRNA_files/", i, "_miRNA_expression_matrix.txt", sep = ""),
              quote = F, sep = "\t")
}

## 3. Download covariate files ----

# technical covariates - extracting information from TCGA barcodes

covrts.tech <- NULL

for (i in cancer.list) {
  covrts.tech[[i]] <- getResults(query.miRNA.download[[i]])[,c("cases", "platform", "sample_type")]
  covrts.tech[[i]] <- covrts.tech[[i]] %>% separate(cases, into = c("project","tss","participant",
                                                                    "sample", "portion", "plate", 
                                                                    "center"), remove = F)
  covrts.tech[[i]] <- covrts.tech[[i]] %>% select(cases, everything())
  covrts.tech[[i]] <- covrts.tech[[i]] %>% select(-participant)
  
}

# saving the technical covariate data for later use

dir.create("data/Technical_covariates_files")

for (i in cancer.list) {
  write.table(covrts.tech[i], file = paste("data/Technical_covariates_files/", 
                                           i, "_technical_covariates", sep = ""))
}

# biological covariates - downloading information from UCSC Xena

dir.create("data/Biological_covariates_files")

data(XenaData)

XenaGenerate(subset = XenaHostNames=="tcgaHub") %>% 
  XenaFilter(filterDatasets = "clinical") %>% 
  XenaFilter(filterDatasets = paste(gsub("TCGA-", "", cancer.list), collapse = "|")) -> df_todo

XenaQuery(df_todo) %>%
  XenaDownload(destdir = "data/Biological_covariates_files/")

system("mv ./data/Biological_covariates_files/*/* ./data/Biological_covariates_files/.")
system("rm -r ./data/Biological_covariates_files/*.sampleMap")

covrts.biol <- lapply(list.files(path = "data/Biological_covariates_files/", full.names = T), 
                      read.delim)

names(covrts.biol) <- cancer.list

# selecting columns with relevant biological information

for (i in cancer.list) {
  covrts.biol[[i]] <- covrts.biol[[i]] %>% 
    filter(sampleID %in% substr(colnames(dat.mir[[i]])[-1], start = 1, stop = 15))
  covrts.biol[[i]] <- covrts.biol[[i]] %>% 
    arrange(match(sampleID, substr(colnames(dat.mir[[i]])[-1], start = 1, stop = 15)))
  covrts.biol[[i]] <- covrts.biol[[i]] %>% 
    select("age_at_initial_pathologic_diagnosis", "gender", "initial_weight",
           "year_of_initial_pathologic_diagnosis", "histological_type",
           "new_tumor_event_after_initial_treatment", 
           contains(c("pathologic_M", "pathologic_N", "pathologic_T", "pathologic_stage")), 
           "sample_type")
  rownames(covrts.biol[[i]]) <- colnames(dat.mir[[i]])[-1]
}

for (i in cancer.list[-c(8,10)]) {
  colnames(covrts.biol[[i]]) <- c("age", "gender", "weight", "year", "histology_type",
                                  "new_tumor", "metastasis", "lymph.node", "tumor_size", "stage",
                                  "sample_type")
}

colnames(covrts.biol[[8]]) <- c("age", "gender", "weight", "year", "histology_type",
                                "new_tumor", "lymph.node", "tumor_size", "sample_type")

colnames(covrts.biol[[10]]) <- c("age", "gender", "weight", "year", "histology_type",
                                "new_tumor", "sample_type")

# saving the biological covariate data for later use

for (i in cancer.list) {
  covrts.biol[[i]] <- covrts.biol[[i]][, which(colMeans(!is.na(covrts.biol[[i]])) > 0.5)]
  write.table(covrts.biol[i], file = paste("data/Biological_covariates_files/", 
                                           i, "_biological_covariates", sep = ""), sep = "\t")
  
}

system("rm ./data/Biological_covariates_files/*_clinicalMatrix")

## 4. Normalization of MiRNA expression matrix ----

# prepare countData and colData for DESeqDataSetFromMatrix function

dat_mir <- lapply(
  list.files(path = "data/MiRNA_files/", 
             pattern = "miRNA_expression_matrix.txt", 
             full.names = T), 
  read.table)

dat_mir <- lapply(dat_mir, function(x){
  rownames(x) <- x[,1]
  x <- x[,-1]
})

condition <- lapply(dat_mir, function(x){
  data.frame(colnames(x), 
             substr(colnames(x), start = 14, stop = 15))
})

condition <- lapply(condition, function(x){
  i <- setNames(object = x, nm = c("cases", "sample_type"))
  i <- mutate(i, sample_type = ifelse(sample_type == "01", "Tumor", "Normal"))
})

# create DESeqDataSet object

dds <- NULL

for (i in 1:length(dat_mir)) {
  dds[[i]] <- DESeqDataSetFromMatrix(countData = dat_mir[[i]],
                                     colData = condition[[i]],
                                     design = ~sample_type)
}

# normalization for sequencing depth and library composition

vsd <- NULL

for (i in 1:length(dds)) {
 vsd[[i]] <- varianceStabilizingTransformation(dds[[i]], blind=F) 
}

## 5. Batch effect estimation ----

## technical covariates 

covrts_tech <- lapply(
  list.files(path = "data/Technical_covariates_files/", 
             pattern = "technical_covariates", full.names = T),
  read.table)
  
covrts_tech <- lapply(covrts_tech, function(x){
  rownames(x) <- gsub("-", ".", x[,1])
  x <- x[,-1]
  colnames(x) <- substring(colnames(x), 11)
  x[, colnames(x)[which(x %>% 
                          summarise_all(funs(n_distinct(., na.rm = T))) > 1)]]
})

# create expSetobj that store the expression values and the technical covariate information
inpData_tech_BRCA <- expSetobj(expression = assay(vsd[[1]]), covariates = covrts_tech[[1]])
inpData_tech_HNSC <- expSetobj(expression = assay(vsd[[2]]), covariates = covrts_tech[[2]])
inpData_tech_KIRC <- expSetobj(expression = assay(vsd[[3]]), covariates = covrts_tech[[3]])
inpData_tech_KIRP <- expSetobj(expression = assay(vsd[[4]]), covariates = covrts_tech[[4]])
inpData_tech_LIHC <- expSetobj(expression = assay(vsd[[5]]), covariates = covrts_tech[[5]])
inpData_tech_LUAD <- expSetobj(expression = assay(vsd[[6]]), covariates = covrts_tech[[6]])
inpData_tech_LUSC <- expSetobj(expression = assay(vsd[[7]]), covariates = covrts_tech[[7]])
inpData_tech_PRAD <- expSetobj(expression = assay(vsd[[8]]), covariates = covrts_tech[[8]])
inpData_tech_THCA <- expSetobj(expression = assay(vsd[[9]]), covariates = covrts_tech[[9]])
inpData_tech_UCEC <- expSetobj(expression = assay(vsd[[10]]),covariates = covrts_tech[[10]])

# Principal variance Component Analysis
pvca_tech_BRCA <- pvcAnaly(exp_datObj = inpData_tech_BRCA,
                           pct_threshold = 0.60,
                           batch_factors = colnames(covrts_tech[[1]]))
pvca_tech_HNSC <- pvcAnaly(exp_datObj = inpData_tech_HNSC,
                           pct_threshold = 0.60,
                           batch_factors = colnames(covrts_tech[[2]]))
pvca_tech_KIRC <- pvcAnaly(exp_datObj = inpData_tech_KIRC,
                           pct_threshold = 0.60,
                           batch_factors = colnames(covrts_tech[[3]]))
pvca_tech_KIRP <- pvcAnaly(exp_datObj = inpData_tech_KIRP,
                           pct_threshold = 0.60,
                           batch_factors = colnames(covrts_tech[[4]]))
pvca_tech_LIHC <- pvcAnaly(exp_datObj = inpData_tech_LIHC,
                           pct_threshold = 0.60,
                           batch_factors = colnames(covrts_tech[[5]]))
pvca_tech_LUAD <- pvcAnaly(exp_datObj = inpData_tech_LUAD,
                           pct_threshold = 0.60,
                           batch_factors = colnames(covrts_tech[[6]]))
pvca_tech_LUSC <- pvcAnaly(exp_datObj = inpData_tech_LUSC,
                           pct_threshold = 0.60,
                           batch_factors = colnames(covrts_tech[[7]]))
pvca_tech_PRAD <- pvcAnaly(exp_datObj = inpData_tech_PRAD,
                           pct_threshold = 0.60,
                           batch_factors = colnames(covrts_tech[[8]]))
pvca_tech_THCA <- pvcAnaly(exp_datObj = inpData_tech_THCA,
                           pct_threshold = 0.60,
                           batch_factors = colnames(covrts_tech[[9]])) 
pvca_tech_UCEC <- pvcAnaly(exp_datObj = inpData_tech_UCEC,
                           pct_threshold = 0.60,
                           batch_factors = colnames(covrts_tech[[10]]))

# Figure.S1: Principal variant component analysis to estimate the proportion of variation in miRNA 
# gene expression contributed by technical variables. The maximum variation in expression is due to 
# “residual” which is due to inter individual differences followed by sample (normal or tumor) for 
# each cancer except for LUAD and PRAD, where the miRNA counts were modeled for batch effect 
# “platform” and “plate” respectively. 

pvca_tech <- list(pvca_tech_BRCA, pvca_tech_HNSC, pvca_tech_KIRC, pvca_tech_KIRP, pvca_tech_LIHC, 
                  pvca_tech_LUAD, pvca_tech_LUSC, pvca_tech_PRAD, pvca_tech_THCA, pvca_tech_UCEC)

pvca_tech <- lapply(pvca_tech, function(x){
  i <- data.frame(x[[2]], round(sapply(x[[1]], "[[", 1)*100, 2))
  i <- setNames(i, nm = c("variables", "variation"))
  i <- i %>% mutate(variables = recode(variables, "resid" = "residual"))
})

names(pvca_tech) <- cancer.list

p <- NULL

for (i in cancer.list) {
  p[[i]] <- ggplot(pvca_tech[[i]], aes(reorder(variables, variation), 
                                       variation, 
                                       label = variation)) +
    geom_bar(stat = "identity", fill = "#F8766D", color = "black") + 
    geom_text(nudge_y = 3, size = 4) +
    theme_bw() + 
    labs(x ="", y = "Weighted Average Percentage Variance", title = gsub("TCGA-", "", i)) + 
    theme(plot.title = element_text(face = "bold", colour = "black", size = 14, hjust = 0.5),
          axis.title = element_text(face = "bold", colour = "black", size = 14),
          axis.text.x = element_text(face = "bold", colour = "black", size = 14, angle = 90, 
                                     hjust = 1, vjust = 0.25),
          axis.text.y = element_text(face = "bold", colour = "black", size = 14))
}

lapply(names(p), function(x) {
  ggsave(filename=paste("results/Step1/covrts_tech_",x,".png",sep=""), 
         plot=p[[x]])
})

## biological covariates 

covrts_biol <- lapply(
  list.files(path = "data/Biological_covariates_files/", 
             pattern = "biological_covariates", full.names = T),
  read.table)

covrts_biol <- lapply(covrts_biol, function(x){
  rownames(x) <- gsub("-", ".", rownames(x))
  colnames(x) <- substring(colnames(x), 11)
  i <- x %>% mutate_all(list(~na_if(.,"")))
  i <- i[, colnames(i)[which(i %>% 
                               summarise_all(funs(n_distinct(., na.rm = T))) > 1)]]
})

# create expSetobj that store the expression values and the biological covariate information
inpData_biol_BRCA <- expSetobj(expression = assay(vsd[[1]]), covariates = covrts_biol[[1]])
inpData_biol_HNSC <- expSetobj(expression = assay(vsd[[2]]), covariates = covrts_biol[[2]])
inpData_biol_KIRC <- expSetobj(expression = assay(vsd[[3]]), covariates = covrts_biol[[3]])
inpData_biol_KIRP <- expSetobj(expression = assay(vsd[[4]]), covariates = covrts_biol[[4]])
inpData_biol_LIHC <- expSetobj(expression = assay(vsd[[5]]), covariates = covrts_biol[[5]])
inpData_biol_LUAD <- expSetobj(expression = assay(vsd[[6]]), covariates = covrts_biol[[6]])
inpData_biol_LUSC <- expSetobj(expression = assay(vsd[[7]]), covariates = covrts_biol[[7]])
inpData_biol_PRAD <- expSetobj(expression = assay(vsd[[8]]), covariates = covrts_biol[[8]])
inpData_biol_THCA <- expSetobj(expression = assay(vsd[[9]]), covariates = covrts_biol[[9]])
inpData_biol_UCEC <- expSetobj(expression = assay(vsd[[10]]),covariates = covrts_biol[[10]])

# Principal variance Component Analysis
pvca_biol_BRCA <- pvcAnaly(exp_datObj = inpData_biol_BRCA,
                           pct_threshold = 0.60,
                           batch_factors = colnames(covrts_biol[[1]]))
pvca_biol_HNSC <- pvcAnaly(exp_datObj = inpData_biol_HNSC,
                           pct_threshold = 0.60,
                           batch_factors = colnames(covrts_biol[[2]]))
pvca_biol_KIRC <- pvcAnaly(exp_datObj = inpData_biol_KIRC,
                           pct_threshold = 0.60,
                           batch_factors = colnames(covrts_biol[[3]]))
pvca_biol_KIRP <- pvcAnaly(exp_datObj = inpData_biol_KIRP,
                           pct_threshold = 0.60,
                           batch_factors = colnames(covrts_biol[[4]]))
pvca_biol_LIHC <- pvcAnaly(exp_datObj = inpData_biol_LIHC,
                           pct_threshold = 0.60,
                           batch_factors = colnames(covrts_biol[[5]]))
pvca_biol_LUAD <- pvcAnaly(exp_datObj = inpData_biol_LUAD,
                           pct_threshold = 0.60,
                           batch_factors = colnames(covrts_biol[[6]]))
pvca_biol_LUSC <- pvcAnaly(exp_datObj = inpData_biol_LUSC,
                           pct_threshold = 0.60,
                           batch_factors = colnames(covrts_biol[[7]]))
pvca_biol_PRAD <- pvcAnaly(exp_datObj = inpData_biol_PRAD,
                           pct_threshold = 0.60,
                           batch_factors = colnames(covrts_biol[[8]]))
pvca_biol_THCA <- pvcAnaly(exp_datObj = inpData_biol_THCA,
                           pct_threshold = 0.60,
                           batch_factors = colnames(covrts_biol[[9]])) 
pvca_biol_UCEC <- pvcAnaly(exp_datObj = inpData_biol_UCEC,
                           pct_threshold = 0.60,
                           batch_factors = colnames(covrts_biol[[10]]))

# Figure S2: Principal variant component analysis to estimate the proportion of variation in miRNA 
# gene expression contributed by biological variables. The maximum variation in expression is due to 
# “residual” which is due to inter individual differences followed by sample (normal or tumor)

pvca_biol <- list(pvca_biol_BRCA, pvca_biol_HNSC, pvca_biol_KIRC, pvca_biol_KIRP, pvca_biol_LIHC, 
                  pvca_biol_LUAD, pvca_biol_LUSC, pvca_biol_PRAD, pvca_biol_THCA, pvca_biol_UCEC)

pvca_biol <- lapply(pvca_biol, function(x){
  i <- data.frame(x[[2]], round(sapply(x[[1]], "[[", 1)*100, 2))
  i <- setNames(i, nm = c("variables", "variation"))
  i <- i %>% mutate(variables = recode(variables, "resid" = "residual"))
})

names(pvca_biol) <- cancer.list

p <- NULL

for (i in cancer.list) {
  p[[i]] <- ggplot(pvca_biol[[i]], aes(reorder(variables, variation), 
                                       variation, 
                                       label = variation)) +
    geom_bar(stat = "identity", fill = "#F8766D", color = "black") + 
    geom_text(nudge_y = 3, size = 4) +
    theme_bw() + 
    labs(x ="", y = "Weighted Average Percentage Variance", title = gsub("TCGA-", "", i)) + 
    theme(plot.title = element_text(face = "bold", colour = "black", size = 14, hjust = 0.5),
          axis.title = element_text(face = "bold", colour = "black", size = 14),
          axis.text.x = element_text(face = "bold", colour = "black", size = 14, angle = 90, 
                                     hjust = 1, vjust = 0.25),
          axis.text.y = element_text(face = "bold", colour = "black", size = 14))
}

lapply(names(p), function(x) {
  ggsave(filename=paste("results/Step1/covrts_biol_",x,".png",sep=""), 
         plot=p[[x]])
})

## 6. Genomic Annotation of clustered MiRNAs ----

# obtaining a list of clustered miRNAs - same chromosome, same strand and
# inter-miRNA distance of 10Kb

download.file(url = "https://www.mirbase.org/ftp/18/genomes/hsa.gff3", 
              destfile = "data/hsa.gff3")

hsa <- read.table("data/hsa.gff3")

hsa <- hsa %>% filter(V3 == "miRNA_primary_transcript") %>%
  select(V1, V4, V5, V7, V9)

hsa <- hsa %>% mutate(V9 = gsub(".*Name=", "", V9))

colnames(hsa) <- c("Chr", "Start", "End", "Strand", "MiRNA")

hsa <- hsa %>% group_by(Chr) %>%
  mutate(difference_1= Start-lag(End, default= dplyr::first(Start)),
         difference_2= abs(End-lead(Start, default= dplyr::first(End)))) %>% 
  mutate(flag_1 = ifelse(difference_1 < 10000, T, F),
         flag_2 = ifelse(difference_2 < 10000, T, F)) 


hsa <- hsa %>% filter(!(difference_1 == 0 & flag_1 == T & flag_2 == F))

hsa <- rbind(hsa %>% filter(flag_1 == F & lead(flag_1) == T),
             hsa %>% filter(flag_1 == T))

hsa <- hsa %>% arrange(Chr, Start)

hsa <- sqldf::sqldf("SELECT hs1.*, SUM(flag)
OVER (PARTITION BY Chr ORDER BY Start) as cluster
FROM ( SELECT hs1.*, ((Start - LAG(Start, 1, Start)
OVER (PARTITION BY Chr ORDER BY Start) > 10000) ) as flag
FROM hsa hs1) hs1;)")

hsa <- hsa%>%
  mutate(cluster = cluster + 1,
         cluster = paste( "chr", Chr , "_" , cluster, sep="")) %>%
  select(c(Chr, Start, End, Strand, MiRNA, cluster))

write.table(hsa,"data/hsa_cluster.txt", sep = "\t")

## 7. Print version information about R, the OS and attached or loaded packages ----
sessionInfo()
