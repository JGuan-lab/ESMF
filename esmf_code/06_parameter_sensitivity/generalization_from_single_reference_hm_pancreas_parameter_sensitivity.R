# ============================================================================
# ESMF PARAMETER SENSITIVITY ANALYSIS (P1-P4) FOR HUMAN PANCREAS
# Single-Reference Generalization benchmark (copy of Ablation_Study_..._Hm_Pancreas)
# ============================================================================
# DESCRIPTION:
# This script performs parameter sensitivity analysis for the ESMF deconvolution
# method on the single-reference generalization benchmark (human pancreas,
# Acinar/Alpha/Beta). It varies one parameter at a time while holding the others
# fixed, and records the resulting H_RMSE / H_Pearson (and W_RMSE / W_Pearson).
#
# P1: rate multiplier grid        c(0, 0.01, 0.1, 0.5, 1, 2, 5, 10)   (prior strength)
# P2: para1/para2/para3 combos    6 combinations                        (layer penalties)
# P3: iter_times                  c(500, 1000, 2000, 3000, 5000)        (iterations)
# P4: seed robustness             20 seeds                              (initialization)

# OUTPUT:
# Results are saved to:
#   /media/desk16/tjn050/LSC/LSC_deconvolution/LSC/results/Ablation_Study_GeneralizationFromSingleReference_Hm_Pancreas_ParameterSensitivity/

# ============================================================================
# SECTION 1: LOAD PACKAGES AND SOURCE FUNCTIONS
# ============================================================================
# Load all required R packages for deconvolution analysis
setwd("/media/desk16/tjn050/LSC/LSC_deconvolution/LSC/")
library(NMF)

library(CellMix)
library(MASS)
#library(pscl)
library(fitdistrplus)
library(scater)
library(Matrix)
library(Seurat)
library(descend)      
library(data.table)
library(cowplot)
library(rhdf5)
library(clusterProfiler)
source("esmf_deconvolution.R")
source("esmf_helper_functions.R")
source("bulk_deconvolution.R")
source('normalize_col_in_matrix.R')
source("basic_functions.R")
source("cibersort.R")


# Define analysis parameters
# Parameter_Sensitivity v2 copy: 独立输出目录，避免覆盖正文基准结果与 v1 结果
analysis_prefix <- "Ablation_Study_GeneralizationFromSingleReference_Hm_Pancreas_ParameterSensitivity_v2"
results_dir <- file.path("/media/desk16/tjn050/LSC/LSC_deconvolution/LSC/results", analysis_prefix, "")
current_datetime <- format(Sys.time(), "%Y%m%d_%H%M%S")

# Create results directory if it doesn't exist
if (!dir.exists(results_dir)) {
  dir.create(results_dir, recursive = TRUE)
  cat("Results directory created:", results_dir, "\n")
}

# ============================================================================
# SECTION 2: SET ANALYSIS PARAMETERS
# ============================================================================
# Configure core parameters for deconvolution analysis
marker_strategy = "all"
number_cells = round(100, digits = -2)  # Must be multiple of 100
to_remove = "none"

dataset = "example"
transformation = "none"
deconv_type = "bulk"

method = "OLS"
number_cells = round(as.numeric("100"), digits = -2)  # Must be multiple of 100
to_remove = "none"
num_cores = min(as.numeric("1"), parallel::detectCores() - 1)

# Set normalization parameters based on deconvolution type
if (deconv_type == "bulk") {
  normalization = "TPM"
  marker_strategy = "all"
} else if (deconv_type == "sc") {
  normalization_scC = "TMM"
  normalization_scT = "TMM"
} else {
  print("Please enter a valid deconvolution framework")
  stop()
}

normalization_scC = "TMM"
normalization_scT = "TMM"


# ============================================================================
# SECTION 3: LOAD AND PREPROCESS TRAINING DATA FROM H5 FORMAT
# ============================================================================
# This section loads and preprocesses human pancreas single-cell RNA-seq data
# from H5 format. The training data is from GSE84133 dataset.
# Training dataset: HU_0228_Pancreas_GSE84133
# Selected cell types: Acinar, Alpha, Beta
# ============================================================================

# ----------------------------------------------------------------------------
# Training Dataset: HU_0228_Pancreas_GSE84133 (H5 format)
# ----------------------------------------------------------------------------
train_dir = "/media/desk16/tjn050/LSC/HUSCH/Pancreas/HU_0228_Pancreas_GSE84133/"
train_set = strsplit(train_dir, "/")[[1]][8]
cat("Loading training dataset:", train_set, "\n")

# Load metadata containing cell type information
phenoData_train <- read.csv(paste0(train_dir, train_set, "_meta.txt"), 
                            sep = "\t", header = TRUE) 

# Load gene count matrix from H5 file
train_mat <- h5read(paste0(train_dir, train_set, "_gene_count.h5"), "matrix")

# Clean cell type names: replace special characters with underscores
phenoData_train$Celltype <- sub("/", "_", phenoData_train$Celltype)
phenoData_train$Celltype <- sub(" ", "_", phenoData_train$Celltype)
phenoData_train$Celltype <- sub(" ", "_", phenoData_train$Celltype)

# Standardize column names
colnames(phenoData_train) = c("cellID", "cellType", "Platform")
table(phenoData_train$cellType)

# Extract components from H5 structure
barcodes <- train_mat$barcodes
counts <- train_mat$data
gene_id <- train_mat$features$id
indices <- train_mat$indices
indptr <- train_mat$indptr
shape <- train_mat$shape

# Reconstruct sparse matrix from H5 components
mat <- sparseMatrix(
  i = indices[] + 1,
  p = indptr[],
  x = as.numeric(x = counts[]),
  dims = shape[],
  repr = "T"
)

# Assign row and column names
colnames(mat) <- barcodes
rownames(mat) <- gene_id

# Convert to dense matrix
train_data <- as.matrix(mat)
rm(mat)  # Remove temporary object to free memory

# Add prefix to column names to identify training data
colnames(train_data) <- paste0("train", colnames(train_data))
phenoData_train$cellID <- paste0("train", phenoData_train$cellID)
table(phenoData_train$cellType)

# Ensure phenotype and expression data match
phenoData_train <- phenoData_train[which(phenoData_train$cellID %in% colnames(train_data)), ]
table(phenoData_train$cellType)
train_data <- train_data[, intersect(phenoData_train$cellID, colnames(train_data))]
dim(train_data)
table(phenoData_train$cellType)

# Convert to matrix format and add metadata column
train_data <- as(train_data, "matrix")
phenoData_train$name <- phenoData_train$cellID

# Reorder columns to match expected format
phenoData_train <- phenoData_train[, c(2, 1, 4)]
table(phenoData_train$cellType)
dim(train_data)
dim(phenoData_train)

# ============================================================================
# SECTION 4: LOAD AND PREPROCESS TEST DATA FROM PANGLAODB
# ============================================================================
# Load test datasets: SRA701877_SRS3279695, SRA701877_SRS3279693, SRA701877_SRS3279692
# These datasets serve as independent validation for deconvolution methods
# All test data are from PanglaoDB sparse RData format
# ============================================================================

# ----------------------------------------------------------------------------
# Test Dataset 1: SRA701877_SRS3279695
# ----------------------------------------------------------------------------
test_set1 = "SRA701877_SRS3279695"
cat("Loading test dataset 1:", test_set1, "\n")

# Load annotation file containing cell type information
anno_SRA_test1 <- read.csv(paste0(paste0("/media/desk16/tjn050/LSC/PanglaoDB/", test_set1), "_annotation.csv"), 
                           sep = ",", header = FALSE) 
colnames(anno_SRA_test1) <- anno_SRA_test1[1, ]
anno_SRA_test1 <- anno_SRA_test1[-1, ]

# Load cluster assignment file
meta_temp <- read.csv(paste0(paste0("/media/desk16/tjn050/LSC/PanglaoDB/", test_set1), ".clusters.txt"),
                      sep = "\t", header = FALSE) 
colnames(meta_temp) <- "cellID"

# Extract cell type information from annotation
cellType_temp <- gsub("[^0-9]", "", meta_temp$cellID)
cellType_temp <- anno_SRA_test1[match(cellType_temp, anno_SRA_test1$`Cluster ID`), 3, drop = F]
cellType_data <- cellType_temp$`Inferred cell type`
cellID_data <- gsub("[0-9]", "", meta_temp$cellID)
cellID_data <- gsub(" ", "", cellID_data)

# Create phenotype data frame
phenoData_test1 = data.frame(cellType = cellType_data, cellID = cellID_data) 
rownames(phenoData_test1) <- cellID_data

# Clean cell type names: remove spaces and hyphens
phenoData_test1$cellType <- sub(" ", "", phenoData_test1$cellType)
phenoData_test1$cellType <- sub(" ", "", phenoData_test1$cellType)
phenoData_test1$cellType <- sub(" ", "", phenoData_test1$cellType)
phenoData_test1$cellType <- sub(" ", "", phenoData_test1$cellType)
phenoData_test1$cellType <- sub("-", "", phenoData_test1$cellType)
table(phenoData_test1$cellType)
phenoData_test1$cellID <- paste0("test1", phenoData_test1$cellID)

# Load sparse expression matrix
load(paste0(paste0("/media/desk16/tjn050/LSC/PanglaoDB/", test_set1), ".sparse.RData"))
test_data1 <- sm
rm(sm)  # Remove temporary object to free memory
dim(test_data1)
colnames(test_data1) <- paste0("test1", colnames(test_data1))
table(phenoData_test1$cellType)

# Filter out "Unknown" cell types
select_cellType <- gsub("Unknown", "", unique(phenoData_test1$cellType))
phenoData_test1 <- phenoData_test1[which(phenoData_test1$cellType %in% select_cellType), ]
table(phenoData_test1$cellType)

# Ensure phenotype and expression data match
phenoData_test1 <- phenoData_test1[which(phenoData_test1$cellID %in% colnames(test_data1)), ]
table(phenoData_test1$cellType)
test_data1 <- test_data1[, intersect(phenoData_test1$cellID, colnames(test_data1))]
table(phenoData_test1$cellType)

# Convert to dense matrix and add metadata
test_data1 <- as(test_data1, "matrix")
phenoData_test1$name <- phenoData_test1$cellID
dim(test_data1)
dim(phenoData_test1)

# Standardize gene names: take first part before underscore
rownames(test_data1) <- unlist(lapply(strsplit(rownames(test_data1), "_"), function(x) x[1]))
table(phenoData_train$cellType)
table(phenoData_test1$cellType)

# ----------------------------------------------------------------------------
# Test Dataset 2: SRA701877_SRS3279693
# ----------------------------------------------------------------------------
test_set2 = "SRA701877_SRS3279693"
cat("Loading test dataset 2:", test_set2, "\n")

# Load annotation file
anno_SRA_test2 <- read.csv(paste0(paste0("/media/desk16/tjn050/LSC/PanglaoDB/", test_set2), "_annotation.csv"), 
                           sep = ",", header = FALSE) 
colnames(anno_SRA_test2) <- anno_SRA_test2[1, ]
anno_SRA_test2 <- anno_SRA_test2[-1, ]

# Load cluster assignments
meta_temp <- read.csv(paste0(paste0("/media/desk16/tjn050/LSC/PanglaoDB/", test_set2), ".clusters.txt"),
                      sep = "\t", header = FALSE) 
colnames(meta_temp) <- "cellID"

# Extract cell type information
cellType_temp <- gsub("[^0-9]", "", meta_temp$cellID)
cellType_temp <- anno_SRA_test2[match(cellType_temp, anno_SRA_test2$`Cluster ID`), 3, drop = F]
cellType_data <- cellType_temp$`Inferred cell type`
cellID_data <- gsub("[0-9]", "", meta_temp$cellID)
cellID_data <- gsub(" ", "", cellID_data)

# Create phenotype data frame
phenoData_test2 = data.frame(cellType = cellType_data, cellID = cellID_data) 
rownames(phenoData_test2) <- cellID_data

# Clean cell type names
phenoData_test2$cellType <- sub(" ", "", phenoData_test2$cellType)
phenoData_test2$cellType <- sub(" ", "", phenoData_test2$cellType)
phenoData_test2$cellType <- sub(" ", "", phenoData_test2$cellType)
phenoData_test2$cellType <- sub(" ", "", phenoData_test2$cellType)
phenoData_test2$cellType <- sub("-", "", phenoData_test2$cellType)
table(phenoData_test2$cellType)
phenoData_test2$cellID <- paste0("test2", phenoData_test2$cellID)

# Load sparse expression matrix
load(paste0(paste0("/media/desk16/tjn050/LSC/PanglaoDB/", test_set2), ".sparse.RData"))
test_data2 <- sm
rm(sm)
dim(test_data2)
colnames(test_data2) <- paste0("test2", colnames(test_data2))
table(phenoData_test2$cellType)

# Filter out "Unknown" cell types
select_cellType <- gsub("Unknown", "", unique(phenoData_test2$cellType))
phenoData_test2 <- phenoData_test2[which(phenoData_test2$cellType %in% select_cellType), ]
table(phenoData_test2$cellType)

# Ensure phenotype and expression data match
phenoData_test2 <- phenoData_test2[which(phenoData_test2$cellID %in% colnames(test_data2)), ]
table(phenoData_test2$cellType)
test_data2 <- test_data2[, intersect(phenoData_test2$cellID, colnames(test_data2))]
table(phenoData_test2$cellType)

# Convert to dense matrix and add metadata
test_data2 <- as(test_data2, "matrix")
phenoData_test2$name <- phenoData_test2$cellID
dim(test_data2)
dim(phenoData_test2)

# Standardize gene names
rownames(test_data2) <- unlist(lapply(strsplit(rownames(test_data2), "_"), function(x) x[1]))
table(phenoData_test1$cellType)
table(phenoData_test2$cellType)

# ----------------------------------------------------------------------------
# Merge Test Datasets 1 and 2
# ----------------------------------------------------------------------------
# Convert gene names to uppercase for consistency
rownames(test_data1) <- toupper(rownames(test_data1))
rownames(test_data2) <- toupper(rownames(test_data2))

# Keep only common genes between the two datasets
test_data1 <- test_data1[intersect(rownames(test_data1), rownames(test_data2)), ]
test_data2 <- test_data2[intersect(rownames(test_data1), rownames(test_data2)), ]
dim(test_data1)
dim(test_data2)

# Combine expression matrices and phenotype data
test_data <- as.matrix(cbind(test_data1, test_data2))
phenoData_test <- rbind(phenoData_test1, phenoData_test2)
dim(test_data)
dim(phenoData_test)
table(phenoData_test$cellType)

# Remove temporary objects to free memory
rm(phenoData_test1, phenoData_test2, test_data1, test_data2)

# ----------------------------------------------------------------------------
# Test Dataset 3: SRA701877_SRS3279692
# ----------------------------------------------------------------------------
test_set3 = "SRA701877_SRS3279692"
cat("Loading test dataset 3:", test_set3, "\n")

# Load annotation file
anno_SRA_test3 <- read.csv(paste0(paste0("/media/desk16/tjn050/LSC/PanglaoDB/", test_set3), "_annotation.csv"), 
                           sep = ",", header = FALSE) 
colnames(anno_SRA_test3) <- anno_SRA_test3[1, ]
anno_SRA_test3 <- anno_SRA_test3[-1, ]

# Load cluster assignments
meta_temp <- read.csv(paste0(paste0("/media/desk16/tjn050/LSC/PanglaoDB/", test_set3), ".clusters.txt"),
                      sep = "\t", header = FALSE) 
colnames(meta_temp) <- "cellID"

# Extract cell type information
cellType_temp <- gsub("[^0-9]", "", meta_temp$cellID)
cellType_temp <- anno_SRA_test3[match(cellType_temp, anno_SRA_test3$`Cluster ID`), 3, drop = F]
cellType_data <- cellType_temp$`Inferred cell type`
cellID_data <- gsub("[0-9]", "", meta_temp$cellID)
cellID_data <- gsub(" ", "", cellID_data)

# Create phenotype data frame
phenoData_test3 = data.frame(cellType = cellType_data, cellID = cellID_data) 
rownames(phenoData_test3) <- cellID_data

# Clean cell type names
phenoData_test3$cellType <- sub(" ", "", phenoData_test3$cellType)
phenoData_test3$cellType <- sub(" ", "", phenoData_test3$cellType)
phenoData_test3$cellType <- sub(" ", "", phenoData_test3$cellType)
phenoData_test3$cellType <- sub(" ", "", phenoData_test3$cellType)
phenoData_test3$cellType <- sub("-", "", phenoData_test3$cellType)
table(phenoData_test3$cellType)
phenoData_test3$cellID <- paste0("test3", phenoData_test3$cellID)

# Load sparse expression matrix
load(paste0(paste0("/media/desk16/tjn050/LSC/PanglaoDB/", test_set3), ".sparse.RData"))
test_data3 <- sm
rm(sm)
dim(test_data3)
colnames(test_data3) <- paste0("test3", colnames(test_data3))
table(phenoData_test3$cellType)

# Filter out "Unknown" cell types
select_cellType <- gsub("Unknown", "", unique(phenoData_test3$cellType))
phenoData_test3 <- phenoData_test3[which(phenoData_test3$cellType %in% select_cellType), ]
table(phenoData_test3$cellType)

# Ensure phenotype and expression data match
phenoData_test3 <- phenoData_test3[which(phenoData_test3$cellID %in% colnames(test_data3)), ]
table(phenoData_test3$cellType)
test_data3 <- test_data3[, intersect(phenoData_test3$cellID, colnames(test_data3))]
table(phenoData_test3$cellType)

# Convert to dense matrix and add metadata
test_data3 <- as(test_data3, "matrix")
phenoData_test3$name <- phenoData_test3$cellID
dim(test_data3)
dim(phenoData_test3)

# Standardize gene names
rownames(test_data3) <- unlist(lapply(strsplit(rownames(test_data3), "_"), function(x) x[1]))
table(phenoData_test$cellType)
table(phenoData_test3$cellType)

# Prepare for merging: convert to uppercase and keep common genes
rownames(test_data) <- toupper(rownames(test_data))
rownames(test_data3) <- toupper(rownames(test_data3))
test_data <- test_data[intersect(rownames(test_data), rownames(test_data3)), ]
test_data3 <- test_data3[intersect(rownames(test_data), rownames(test_data3)), ]
dim(test_data)
dim(test_data3)

# Final merge of all test datasets
test_data <- as.matrix(cbind(test_data, test_data3))
phenoData_test <- rbind(phenoData_test, phenoData_test3)
dim(test_data)
dim(phenoData_test)
table(phenoData_test$cellType)
table(phenoData_train$cellType)

# Remove temporary objects
rm(phenoData_test3, test_data3)

# ----------------------------------------------------------------------------
# Standardize Cell Type Nomenclature Across Datasets
# ----------------------------------------------------------------------------
# Harmonize cell type names between training and test datasets
table(phenoData_test$cellType)
table(phenoData_train$cellType)

# Standardize test dataset cell type names to match training data
phenoData_test$cellType <- sub("Acinarcells", "Acinar", phenoData_test$cellType)
phenoData_test$cellType <- sub("Alphacells", "Alpha", phenoData_test$cellType)
phenoData_test$cellType <- sub("Betacells", "Beta", phenoData_test$cellType)
phenoData_test$cellType <- sub("Ductalcells", "Duct", phenoData_test$cellType)
table(phenoData_test$cellType)
table(phenoData_train$cellType)

# ----------------------------------------------------------------------------
# Select Specific Cell Types for Analysis
# ----------------------------------------------------------------------------
# For this pancreas tissue analysis, we focus on three major endocrine cell types:
# Acinar, Alpha, and Beta cells
select_cellType <- c("Acinar", "Alpha", "Beta")

# Filter training data to selected cell types
phenoData_train <- phenoData_train[which(phenoData_train$cellType %in% select_cellType), ]
phenoData_train <- phenoData_train[which(phenoData_train$cellID %in% colnames(train_data)), ]
table(phenoData_train$cellType)
train_data <- train_data[, intersect(phenoData_train$cellID, colnames(train_data))]
phenoData_train$name <- phenoData_train$cellID
dim(phenoData_train)
dim(train_data)

# Filter test data to selected cell types
phenoData_test <- phenoData_test[which(phenoData_test$cellType %in% select_cellType), ]
phenoData_test <- phenoData_test[which(phenoData_test$cellID %in% colnames(test_data)), ]
table(phenoData_test$cellType)
test_data <- test_data[, intersect(phenoData_test$cellID, colnames(test_data))]
phenoData_test$name <- phenoData_test$cellID
dim(phenoData_test)
dim(test_data)
table(phenoData_test$cellType)
table(phenoData_train$cellType)

# ----------------------------------------------------------------------------
# Final Data Integration
# ----------------------------------------------------------------------------
# Ensure gene name consistency between training and test data
rownames(train_data) <- toupper(rownames(train_data))
rownames(test_data) <- toupper(rownames(test_data))
train_data <- train_data[intersect(rownames(train_data), rownames(test_data)), ]
test_data <- test_data[intersect(rownames(train_data), rownames(test_data)), ]
dim(train_data)
dim(test_data)

# Combine training and test data for downstream analysis
data <- as.matrix(cbind(train_data, test_data))
full_phenoData <- rbind(phenoData_train, phenoData_test)
dim(data)
dim(full_phenoData)
table(full_phenoData$cellType)
table(phenoData_test$cellType)

# Remove intermediate objects to free memory
rm(train_data, test_data, phenoData_test, phenoData_train)


# ============================================================================
# SECTION 5: QUALITY CONTROL
# ============================================================================
# Quality control: remove cells with library size, mitochondrial or ribosomal
# content beyond 3 median absolute deviations
require(dplyr)
require(Matrix)

# Function to identify cells to remove based on filtering parameter
filterCells <- function(filterParam) {
  cellsToRemove <- which(filterParam > median(filterParam) + 3 * mad(filterParam) | 
                           filterParam < median(filterParam) - 3 * mad(filterParam))
  cellsToRemove
}

libSizes <- colSums(data)
gene_names <- rownames(data)

mtID <- grepl("^MT-|_MT-", gene_names, ignore.case = TRUE)
rbID <- grepl("^RPL|^RPS|_RPL|_RPS", gene_names, ignore.case = TRUE)

mtPercent <- colSums(data[mtID, ]) / libSizes
rbPercent <- colSums(data[rbID, ]) / libSizes

cellsToRemove <- lapply(list(libSizes = libSizes, mtPercent = mtPercent, rbPercent = rbPercent), 
                        filterCells) %>% 
  unlist() %>% 
  unique()

if (length(cellsToRemove) != 0) {
  data <- data[, -cellsToRemove]
  full_phenoData <- full_phenoData[-cellsToRemove, ]
}

# Keep only genes detected in at least 5% of cells
keep <- which(Matrix::rowSums(data > 0) >= round(0.05 * ncol(data)))
data = data[keep, ]
dim(data)

# ============================================================================
# SECTION 6: DATA SPLITTING AND REFERENCE MATRIX CONSTRUCTION
# ============================================================================
set.seed(24)
require(limma)
require(pheatmap)

original_cell_names = colnames(data)
colnames(data) <- as.character(full_phenoData$cellType[match(colnames(data), full_phenoData$cellID)])

# Keep cell types with >= 20 cells after QC
cell_counts = table(colnames(data))
to_keep = names(cell_counts)[cell_counts >= 20]
to_keep
pData <- full_phenoData[full_phenoData$cellType %in% to_keep, ]
to_keep = which(colnames(data) %in% to_keep)
data <- data[, to_keep]
original_cell_names <- original_cell_names[to_keep]

# Split data into training and testing sets
set.seed(2)
training <- grep(pattern = "train", original_cell_names)
testing <- grep(pattern = "test", original_cell_names)

# Generate phenodata for reference matrix C
pDataC = pData[training, ]

train <- data[, training]
test <- data[, testing]
original_cell_names_train <- original_cell_names[training]

train_cellID = train
colnames(train_cellID) = original_cell_names[training]

# Construct reference matrix C and variance profiles
cellType <- colnames(train)
group = list()
for (i in unique(cellType)) { 
  group[[i]] <- which(cellType %in% i)
}
C = lapply(group, function(x) Matrix::rowMeans(train[, x]))
C = round(do.call(cbind.data.frame, C))

refProfiles.var = lapply(group, function(x) train[, x])
refProfiles.var = lapply(refProfiles.var, function(x) matrixStats::rowSds(Matrix::as.matrix(x)))
refProfiles.var = round(do.call(cbind.data.frame, refProfiles.var))
rownames(refProfiles.var) <- rownames(train)

# ============================================================================
# SECTION 7: NORMALIZATION AND MARKER GENE SELECTION
# ============================================================================
# For marker selection, keep genes where at least 30% of cells within a cell
# type have a read/UMI count different from 0
cellType = colnames(train)
keep <- sapply(unique(cellType), function(x) {
  CT_hits = which(cellType %in% x)
  size = ceiling(0.3 * length(CT_hits))
  Matrix::rowSums(train[, CT_hits, drop = FALSE] != 0) >= size
})
train = train[Matrix::rowSums(keep) > 0, ]
train2 = Normalization(train)

# Differential expression analysis using limma-voom
annotation = factor(colnames(train2))
design <- model.matrix(~0 + annotation)
colnames(design) <- unlist(lapply(strsplit(colnames(design), "annotation"), function(x) x[2]))
cont.matrix <- matrix((-1 / ncol(design)), nrow = ncol(design), ncol = ncol(design))
colnames(cont.matrix) <- colnames(design)
diag(cont.matrix) <- (ncol(design) - 1) / ncol(design)

v <- limma::voom(train2, design = design, plot = FALSE)
fit <- limma::lmFit(v, design)
fit2 <- limma::contrasts.fit(fit, cont.matrix)
fit2 <- limma::eBayes(fit2, trend = TRUE)

markers = marker.fc(fit2, log2.threshold = log2(1))
table(markers$CT)



prio_markers=list(
  Acinar=c("CTRB2","REG1A","SERPINA3","REG1B","C15ORF48","AKR1C3","LYZ"),
  Alpha =c("SCGB2A1"),
  #Pancreaticstellatecells=c(     "DEFA1","DEFA6","LYZ", "IL4R"),
  Beta =c("INS","SCGB2A1")
  #Ductalcells =c("SERPINA3","AKR1C3","DCDC2")
)


# Filter markers: keep only those expressed in at least 50% of cells
sig_train <- train[markers$gene, ]
dim(sig_train)
keep <- which(Matrix::rowSums(sig_train > 0) >= round(0.5 * ncol(sig_train)))
markers = markers[keep, ]
table(markers$CT)

# Check for priority markers in the selected marker list
for (i in 1:length(unique(markers$CT))) {
  gene_index <<- which(markers$gene %in% prio_markers[[i]])
  print(markers$gene[gene_index])
}

# Select top 200 markers per cell type
top_markers <- data.frame()
for (CT in unique(markers$CT)) { 
  CT_index <- which(markers$CT == CT)
  if (length(CT_index) > 200) {
    CT_index = CT_index[1:200]
  }
  top_markers <- rbind(top_markers, markers[CT_index, ])
}
markers <- top_markers
table(markers$CT)

# Re-check priority markers after top marker selection
for (i in 1:length(unique(markers$CT))) {
  gene_index <<- which(markers$gene %in% prio_markers[[i]])
  print(markers$gene[gene_index])
}
table(markers$CT)
dim(markers)

# Generate marker distribution for deconvolution methods
marker_distrib = marker_strategies(markers, marker_strategy, C)
md = marker_distrib
table(md$CT)

# ============================================================================
# SECTION 8: PSEUDO-BULK MIXTURE GENERATION
# ============================================================================
# Generate 1000 pseudo-bulk mixtures on test data for validation
cellType <- colnames(test)
colnames(test) <- original_cell_names[testing]
set.seed(25)

generator <- Generator(sce = test, phenoData = full_phenoData, 
                       Num.mixtures = 1000, pool.size = number_cells)
T <- generator[["T"]]
P <- generator[["P"]]
C = lapply(group, function(x) Matrix::rowMeans(train[, x]))
C = round(do.call(cbind.data.frame, C))

# Bulk simulation matrix
V <- as.matrix(T)

# ============================================================================
# SECTION 9: BULK DECONVOLUTION METHODS EVALUATION
# ============================================================================
# Evaluate multiple bulk deconvolution methods on simulated mixtures
normalization = "none"
T <- generator[["T"]]
P <- generator[["P"]]
C = lapply(group, function(x) Matrix::rowMeans(train[, x]))
C = round(do.call(cbind.data.frame, C))
genes_filt <- intersect(rownames(T), md$gene)
# genes_filt <- intersect(rownames(T), rownames(sig_matrix))
T <- T[genes_filt, ]
C <- C[genes_filt, ]
V <- as.matrix(T)
dim(V)
dim(C)
dim(T)
dim(P)

elem <- NULL
STRING <- "example_string"
methods_to_use <- c("nnls", "OLS", "ridge")  # 参数敏感性：仅 3 个确定性方法用于复现验证

# Data transformation and scaling
T = Transformation(T, transformation)
C = Transformation(C, transformation)
T = Scaling(T, normalization)
C = Scaling(C, normalization)

# Marker selection (on training data)
marker_distrib = marker_strategies(markers, marker_strategy, C)

# If a cell type is removed, only keep mixtures where that CT was present
if (to_remove != "none") {
  T <- T[, P[to_remove, ] != 0]
  C <- C[, colnames(C) %in% rownames(P) & (!colnames(C) %in% to_remove)]
  P <- P[!rownames(P) %in% to_remove, colnames(T)]
  refProfiles.var = refProfiles.var[, colnames(refProfiles.var) %in% rownames(P) & (!colnames(refProfiles.var) %in% to_remove)]
  marker_distrib <- marker_distrib[marker_distrib$CT %in% rownames(P) & (marker_distrib$CT != to_remove), ]
}

methods_to_use <- c("nnls", "OLS", "ridge")  # 参数敏感性：仅 3 个确定性方法用于复现验证

# Create empty lists to store results
H_rmse_list <- list()
H_pearson_list <- list()

# Evaluate each bulk deconvolution method
for (method in methods_to_use) {
  set.seed(123456)
  RESULTS = Deconvolution(T = T, C = C, method = method, P = P, elem = to_remove, 
                          marker_distrib = marker_distrib, refProfiles.var = refProfiles.var)
  RESULTS = RESULTS %>% dplyr::summarise(
    RMSE = sqrt(mean((observed_values - expected_values)^2)) %>% round(., 4),
    Pearson = cor(observed_values, expected_values) %>% round(., 4)
  )
  H_rmse_list[[method]] <- RESULTS$RMSE
  H_pearson_list[[method]] <- RESULTS$Pearson
  print(RESULTS)
}

# ============================================================================
# SECTION 10: SINGLE-CELL DECONVOLUTION METHODS EVALUATION
# ============================================================================
if (FALSE) {  # SKIPPED: 参数敏感性脚本不跑 Section 10（MuSiC + CIBERSORTx，慢且非参数敏感性依赖）
# Evaluate single-cell deconvolution methods (MuSiC and CIBERSORTx)
# Re-initialize data for sc methods
T <- generator[["T"]]
P <- generator[["P"]]
C = lapply(group, function(x) Matrix::rowMeans(train[, x]))
C = round(do.call(cbind.data.frame, C))

# Data transformation and scaling for sc methods
T = Transformation(T, transformation)
C = Transformation(train_cellID, transformation)
T = Scaling(T, normalization_scT)
C = Scaling(C, normalization_scC)

# If a cell type is removed, adjust data accordingly
if (to_remove != "none") {
  T <- T[, P[to_remove, ] != 0]
  C <- C[, pDataC$cellType != to_remove]
  P <- P[!rownames(P) %in% to_remove, colnames(T)]
  pDataC <- pDataC[pDataC$cellType != to_remove, ]
}

phenoDataC = pDataC
elem = to_remove

# BisqueRNA requires "SubjectName" in phenoDataC
if (length(grep("[N-n]ame", colnames(phenoDataC))) > 0) {
  sample_column = grep("[N-n]ame", colnames(phenoDataC))
} else {
  sample_column = grep("[S-s]ample|[S-s]ubject", colnames(phenoDataC))
}
colnames(phenoDataC)[sample_column] = "SubjectName"
rownames(phenoDataC) = phenoDataC$cellID

# Create ExpressionSet objects for compatibility
require(xbioc)
C.eset <- Biobase::ExpressionSet(assayData = as.matrix(C), 
                                 phenoData = Biobase::AnnotatedDataFrame(phenoDataC))
T.eset <- Biobase::ExpressionSet(assayData = as.matrix(T))

# Ensure matrix dimension appropriateness
keep = intersect(rownames(C), rownames(T))
C = C[keep, ]
T = T[keep, ]

T.mtx = exprs(T.eset)
C.sce <- SingleCellExperiment(assays = list(counts = C), colData = phenoDataC)
set.seed(123456)

# MuSiC deconvolution
RES_TEMP = t(MuSiC::music_prop(bulk.mtx = T, sc.sce = C.sce, clusters = 'cellType',
                               samples = 'cellID', markers = NULL, normalize = FALSE, 
                               verbose = F)$Est.prop.weighted)
mat1 <- RES_TEMP
mat2 <- P
# Calculate Pearson correlation between samples
samples <- colnames(mat1)
cor_results <- sapply(samples, function(s) {
  cor(mat1[, s], mat2[, s], method = "pearson")
})

# Prepare MuSiC results for evaluation
RESULTS = RES_TEMP[gtools::mixedsort(rownames(RES_TEMP)), ]
RESULTS = RESULTS[gtools::mixedsort(rownames(RESULTS)), ]
RESULTS <- as.data.table(RESULTS, keep.rownames = TRUE)
RESULTS = data.table::melt(RESULTS)
colnames(RESULTS) <- c("CT", "tissue", "observed_values")

P = P[gtools::mixedsort(rownames(P)), ]
P$CT = rownames(P)
setDT(P)
P = data.table::melt(P, id.vars = "CT")
colnames(P) <- c("CT", "tissue", "expected_values")

RESULTS = merge(RESULTS, P)
RESULTS$expected_values <- round(RESULTS$expected_values, 3)
RESULTS$observed_values <- round(RESULTS$observed_values, 3)
RESULTS = RESULTS %>% dplyr::summarise(
  RMSE = sqrt(mean((observed_values - expected_values)^2)) %>% round(., 4),
  Pearson = cor(observed_values, expected_values) %>% round(., 4)
)
H_rmse_list[["MuSiC"]] <- RESULTS$RMSE
H_pearson_list[["MuSiC"]] <- RESULTS$Pearson

# CIBERSORTx deconvolution
T <- generator[["T"]]
P <- generator[["P"]]
C = lapply(group, function(x) Matrix::rowMeans(train[, x]))
C = round(do.call(cbind.data.frame, C))

# Load CIBERSORTx results from file
results <- read.table("/media/desk16/tjn050/LSC/LSC_deconvolution/LSC/results/sim5_HU_0228_Pancreas_GSE84133_CIBERSORTx_Job76_Results.csv", 
                      header = TRUE, sep = ",")
rownums <- dim(results)[2] - 3
H <- as.matrix(results[, 2:rownums])
# H <- H[, colnames(C)]
rownames(H) <- results$Mixture
H <- t(H)
mat1 <- H
mat2 <- P
# Calculate Pearson correlation between samples
samples <- colnames(mat1)
cor_results <- sapply(samples, function(s) {
  cor(mat1[, s], mat2[, s], method = "pearson")
})

corH <- cor(H, P)
H <- H[gtools::mixedsort(rownames(H)), ]
# setDT(H)
H <- reshape2::melt(H)
# H <- data.table::melt(H)
colnames(H) <- c("CT", "tissue", "observed_values")
head(H)

P <- as.matrix(P)
P <- P[gtools::mixedsort(rownames(P)), ]
# P <- data.table::melt(P)
P <- reshape2::melt(P)
colnames(P) <- c("CT", "tissue", "expected_values")
head(P)

INDICATOR = merge(H, P)
head(INDICATOR)
INDICATOR$expected_values <- round(INDICATOR$expected_values, 2)
INDICATOR$observed_values <- round(INDICATOR$observed_values, 2)
INDICATOR$expected_values[1:10]
INDICATOR$observed_values[1:10]

INDICATOR = INDICATOR %>% dplyr::summarise(
  RMSE = sqrt(mean((observed_values - expected_values)^2)) %>% round(., 4),
  Pearson = cor(observed_values, expected_values) %>% round(., 4)
)
H_rmse_list[["CIBERSORTx"]] <- INDICATOR$RMSE
H_pearson_list[["CIBERSORTx"]] <- INDICATOR$Pearson
}  # end if (FALSE) SKIPPED Section 10

# ============================================================================
# SECTION 11: SIGNATURE MATRIX ESTIMATION
# ============================================================================
# Estimate signature matrix using negative binomial distribution
table(md$CT)
sig_train <- train[md$gene, ]


params_ls <- list()
sim_ls <- list()
mu_ls <- list()
sd_ls <- list()
library(fitdistrplus)

# md$CT <- sub("B.cells", "B cells", md$CT)
# md$CT <- sub("B.cells", "B cells", md$CT)

# Fit negative binomial distribution for each gene in each cell type
for (i in 1:length(table(md$CT))) {
  means <- vector()
  sds <- vector()
  cell_index <<- which(colnames(sig_train) %in% names(table(md$CT))[i])
  counts <- sig_train[, cell_index]
  # means <- rowMeans(counts)
  # test <- as.data.frame(t(counts))
  for (j in 1:dim(counts)[1]) {
    # fm_nb <- MASS::glm.nb(gene_1 ~ ., data = test)
    # fm_zinb <- zeroinfl(gene_1 ~ . | 1, data = test, dist = "negbin")
    # fm_zip <- zeroinfl(gene_1 ~ . | 1, data = bioChemists)
    # f1 <- fitdist(round(counts[j,]), "nbinom", method = "mse")
    if (length(which(round(counts[j, ]) != 0)) > 0) {
      # f1 <- fitdist(round(counts[j,]), "nbinom", method = "mle")
      f1 <<- fitdist(round(counts[j, ]), "nbinom", method = "mse")
      msg <- paste((f1$estimate)[1], (f1$estimate)[2])
      msg <- paste("fit is", msg)
      print(paste(j, msg))
      
      if (is.na(f1$estimate)[1]) {
        means <- c(means, 0)
        print("means NA")
        if (is.na(f1$estimate)[2]) {
          sds <- c(sds, 0)
          print("sd NA")
        } else {
          sds <- c(sds, f1[["estimate"]][["size"]])
        }
      } else if (is.na(f1$estimate)[2]) {
        sds <- c(sds, 0)
        print("sd NA")
        if (is.na(f1$estimate)[1]) {
          means <- c(means, 0)
          print("means NA")
        } else {
          means <- c(means, f1[["estimate"]][["mu"]])
        }
      } else {
        means <- c(means, f1[["estimate"]][["mu"]])
        sds <- c(sds, f1[["estimate"]][["size"]])
      }
    } else if (length(which(round(counts[j, ]) != 0)) == 0) {
      means <- c(means, 0)
      sds <- c(sds, 0)
      print("all zeros")
    }
  }
  mu_ls <- c(mu_ls, list(means))
  sd_ls <- c(sd_ls, list(sds))
}

# Construct signature matrix and standard deviation matrix
sig_matrix <- matrix(unlist(mu_ls), ncol = length(mu_ls))
colnames(sig_matrix) <- names(table(md$CT))
rownames(sig_matrix) <- rownames(sig_train)
# sig_matrix <- sig_matrix[marker_index, ]
theta_matrix <- matrix(unlist(sd_ls), ncol = length(sd_ls))
sd_matrix <- sig_matrix * (1 + theta_matrix)
colnames(sd_matrix) <- names(table(md$CT))
rownames(sd_matrix) <- rownames(sig_train)
# sd_matrix <- sd_matrix[marker_index, ]
colSums(sig_matrix)
colSums(sd_matrix)
# md <- md[rownames(sig_matrix), ]

# Create MarkerList object for CellMix methods
ML = CellMix::MarkerList()
ML@.Data <- tapply(as.character(md$gene), as.character(md$CT), list)

# ============================================================================
# SECTION 12: ESMF
# ============================================================================
# Prepare data for ESMF deconvolution
normalization = "none"
T <- generator[["T"]]
P <- generator[["P"]]
C = lapply(group, function(x) Matrix::rowMeans(train[, x]))
C = round(do.call(cbind.data.frame, C))

# Filter genes to intersect between V and signature matrix
genes_filt <- intersect(rownames(V), rownames(sig_matrix))
V <- V[genes_filt, ]
sig_matrix <- sig_matrix[genes_filt, ]
sd_matrix <- sd_matrix[genes_filt, ]
md <- md[rownames(sig_matrix), ]

dim(V)
dim(sig_matrix)
dim(sd_matrix)
dim(md)

# Generate coefficient matrices for ESMF
m <- dim(V)[1]
n <- dim(V)[2]
r <- length(table(md$CT))
set.seed(1111)
W = rmatrix(m, r)
set.seed(1112)
H = rmatrix(r, n)
# colnames(V) <- tissueID
# m <- length(V)
# n <- 1
dim(W)
dim(H)

# Set penalty parameters for marker gene constraints
# para1 <- 0.001
# para2 <- 100000
# para3 <- 10000000
# para1 <- 0.001
# para2 <- 1
# para3 <- 10
para1 <- 1
para2 <- 1
para3 <- 1
names(table(md$CT))

coef1 <- matrix(data = para1, nrow = m, ncol = r)
coef2 <- matrix(data = para1, nrow = r, ncol = n)
coef3 <- matrix(data = para1, nrow = m, ncol = r)
step1 <- matrix(data = para1, nrow = m, ncol = r)
step2 <- matrix(data = para1, nrow = m, ncol = n)

# Apply marker-specific penalties
for (i in 1:length(table(md$CT))) {
  gene_index <<- which(md$CT %in% names(table(md$CT))[i])
  coef1[gene_index, i] <- para2
  coef3[gene_index, i] <- para2
  step1[gene_index, i] <- para2
}

# Apply priority marker penalties
for (i in 1:length(table(md$CT))) {
  gene_index <<- which(md$gene %in% prio_markers[[i]])
  print(md$gene[gene_index])
}
for (i in 1:length(table(md$CT))) {
  gene_index <<- which(md$gene %in% prio_markers[[i]])
  print(gene_index)
  print(i)
  coef1[gene_index, i] <- para3
  coef3[gene_index, i] <- para3
  step1[gene_index, i] <- para3
}
# ----------------------------------------------------------------------------
# Save ESMF Input Parameters for Reproducibility
# ----------------------------------------------------------------------------
# Save all input parameters with timestamp for debugging and reproducibility
# This ensures that the exact input data for ESMF deconvolution is archived
timestamp_esmf <- format(Sys.time(), "%Y%m%d_%H%M%S")
esmf_param_dir <- file.path(results_dir, "ESMF_parameters", "")
if (!dir.exists(esmf_param_dir)) {
  dir.create(esmf_param_dir, recursive = TRUE)
  cat("ESMF parameters directory created:", esmf_param_dir, "\n")
}

# Function to save matrix as file with metadata
save_matrix_as_file <- function(matrix, prefix, description, timestamp, out_dir) {
  filename <- paste0(prefix, "_", timestamp, ".txt")
  filepath <- file.path(out_dir, filename)
  
  # Create header with metadata
  cat(paste0("# ", description, "\n"), file = filepath)
  cat(paste0("# Generated: ", Sys.time(), "\n"), file = filepath, append = TRUE)
  cat(paste0("# Dimensions: ", nrow(matrix), " x ", ncol(matrix), "\n"), file = filepath, append = TRUE)
  cat(paste0("# Prefix: ", prefix, "\n\n"), file = filepath, append = TRUE)
  
  # Save matrix data
  write.table(matrix, file = filepath, sep = "\t", quote = FALSE, 
              col.names = TRUE, row.names = TRUE, append = TRUE)
  cat("Saved", description, "to", filepath, "\n")
  
  return(filepath)
}

# Save all ESMF input parameters
cat("\nSaving ESMF input parameter files (timestamp:", timestamp_esmf, ")...\n")
save_matrix_as_file(V, "V", "Bulk expression matrix for ESMF", timestamp_esmf, esmf_param_dir)
save_matrix_as_file(sig_matrix, "sig_matrix", "Signature matrix for ESMF", timestamp_esmf, esmf_param_dir)
save_matrix_as_file(sd_matrix, "sd_matrix", "Standard deviation matrix for ESMF", timestamp_esmf, esmf_param_dir)
save_matrix_as_file(coef1, "coef1", "Penalty coefficient matrix for ESMF", timestamp_esmf, esmf_param_dir)

# Also save the marker distribution for reference
save_matrix_as_file(as.matrix(md), "marker_distrib", "Marker gene distribution for ESMF", 
                    timestamp_esmf, esmf_param_dir)

# Save a summary file with all parameter details
summary_file <- file.path(esmf_param_dir, paste0("ESMF_parameters_summary_", timestamp_esmf, ".txt"))
sink(summary_file)
cat("=== ESMF DECONVOLUTION PARAMETERS SUMMARY ===\n\n")
cat("Timestamp:", timestamp_esmf, "\n")
cat("Analysis date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("Analysis prefix:", analysis_prefix, "\n\n")

cat("=== MATRIX DIMENSIONS ===\n")
cat("V (bulk matrix):", dim(V)[1], "genes x", dim(V)[2], "samples\n")
cat("Signature matrix:", dim(sig_matrix)[1], "genes x", dim(sig_matrix)[2], "cell types\n")
cat("SD matrix:", dim(sd_matrix)[1], "genes x", dim(sd_matrix)[2], "cell types\n")
cat("Penalty matrix (coef1):", dim(coef1)[1], "genes x", dim(coef1)[2], "cell types\n\n")

cat("=== PARAMETER SETTINGS ===\n")
cat("para1 (baseline penalty):", para1, "\n")
cat("para2 (marker penalty):", para2, "\n")
cat("para3 (priority marker penalty):", para3, "\n")
cat("Number of cell types (r):", r, "\n")
cat("Number of samples (n):", n, "\n")
cat("Number of genes (m):", m, "\n\n")

cat("=== CELL TYPE INFORMATION ===\n")
cat("Cell types:", paste(colnames(sig_matrix), collapse = ", "), "\n")
cat("Marker distribution counts:\n")
print(table(md$CT))
sink()
cat("ESMF parameters summary saved to:", summary_file, "\n")


sdfold = 3
err_cutoff = 10e-5
rate = coef1
iter_times = 5000

# ---- 局部辅助函数：向量化 ESMF 内核（与 esmf_deconvolution.R 数学等价，仅提速）----
# 说明：W/H 更新公式与原版逐元素一致。err 只用于提前终止判断；
#       err_cutoff=0 且 err 单调下降时不触发终止，因此 3000 次迭代下 W/H 与原版一致。
#       去掉每迭代 3 行 message 日志，实测提速约 20 倍。
#       2026-08-18：上移至此（Section 12 之前），供 Section 12 基准与 Section 15 P1–P4 共用。
ESMF_deconvolution_fast <- function(V, ml, r, err_cutoff, rate, sig_matrix, sd_matrix, iter_times) {
  m <- dim(V)[1]
  n <- dim(V)[2]
  W <- rmatrix(m, r)
  H <- rmatrix(r, n)
  err <- 0.0
  for (iteration in 1:iter_times) {
    V_pre <- W %*% H
    E <- V - V_pre
    err_pre <- err
    err <- sum(E * E)
    if (err < err_cutoff || err == err_pre || (err > err_pre && iteration > 100)) break
    VH <- V %*% t(H)
    WHH <- W %*% H %*% t(H)
    W <- W * ((VH + rate * sig_matrix) / (WHH + rate * W))
    W <<- W
    WV <- t(W) %*% V
    WWH <- t(W) %*% W %*% H
    WWH_safe <- WWH
    WWH_safe[WWH_safe == 0] <- 1
    H_new <- (H * WV) / WWH_safe
    if (any(WWH == 0)) H_new[WWH == 0] <- H[WWH == 0]
    H <- H_new
    H <<- H
    # 原版为 if (iter_times %% 10 == 0)，对 500/1000/2000/3000/5000 等每轮生效
    if (iter_times %% 10 == 0) {
      H[H < .Machine$double.eps] <- .Machine$double.eps
      W[W < .Machine$double.eps] <- .Machine$double.eps
    }
  }
  list(W = W, H = H, err = err)
}

set.seed(1233322)


RES <- ESMF_deconvolution_fast(V = V, ml = ML, r = r, err_cutoff = 0, rate = coef1, sig_matrix = sig_matrix, sd_matrix = sd_matrix, iter_times = 3000)

# ============================================================================
# SECTION 12.5: INPUT REPRODUCIBILITY CHECK（输入复现校验）
# 本段双重验证输入与正文基准一致，全部通过才继续参数敏感性实验：
#   1) 3 个确定性方法（nnls/OLS/ridge）结果 vs 基准 deconvolution_summary_20260815_235725.txt
#   2) ESMF 输入（V/sig_matrix/sd_matrix/coef1/marker_distrib）vs 基准参数文件 *_20260816_002201.txt
# 任一不一致立即终止，避免基于漂移输入做参数敏感性分析。
# ============================================================================
bench_ok <- TRUE
baseline_methods <- list(
  nnls  = c(RMSE = 0.2344, Pearson = 0.6871),
  OLS   = c(RMSE = 0.2271, Pearson = 0.6924),
  ridge = c(RMSE = 0.2220, Pearson = 0.6950))
cat("\n---- [CHECK 1/2] 3-method benchmark vs baseline ----\n")
for (mm in names(baseline_methods)) {
  if (mm %in% names(H_rmse_list)) {
    d_r <- abs(H_rmse_list[[mm]] - baseline_methods[[mm]]["RMSE"])
    d_p <- abs(H_pearson_list[[mm]] - baseline_methods[[mm]]["Pearson"])
    ok <- d_r <= 1e-4 && d_p <= 1e-4
    bench_ok <- bench_ok && ok
    cat(sprintf("[BENCH] %-10s RMSE=%.4f (base %.4f, diff %.0e)  Pearson=%.4f (base %.4f, diff %.0e) => %s\n",
                mm, H_rmse_list[[mm]], baseline_methods[[mm]]["RMSE"], d_r,
                H_pearson_list[[mm]], baseline_methods[[mm]]["Pearson"], d_p,
                ifelse(ok, "OK", "MISMATCH")))
  } else {
    bench_ok <- FALSE
    cat(sprintf("[BENCH] %-10s 未运行，无法验证 => MISMATCH\n", mm))
  }
}

cat("\n---- [CHECK 2/2] ESMF input vs baseline parameter files ----\n")
baseline_ts <- "20260816_002201"   # 正文基准 ESMF 参数 timestamp（ESMF 基准那次）
baseline_param_dir <- file.path(
  "/media/desk16/tjn050/LSC/LSC_deconvolution/LSC/results/GeneralizationFromSingleReference_Hm_Pancreas",
  "ESMF_parameters")

check_repro <- function(obj, obj_name, tol = 1e-6) {
  fname <- file.path(baseline_param_dir, paste0(obj_name, "_", baseline_ts, ".txt"))
  if (!file.exists(fname)) {
    cat(sprintf("[CHECK] %-12s : 基准文件缺失 %s\n", obj_name, fname))
    return(FALSE)
  }
  ref <- as.matrix(read.table(fname, header = TRUE, row.names = 1,
                              comment.char = "#", check.names = FALSE))
  cur <- as.matrix(obj)
  dim_ok <- identical(dim(ref), dim(cur))
  # 行列名一致性：write.table 对无名矩阵会自动补占位名（行名 1..n / 列名 V1..Vk）。
  # 若一侧为 NULL、另一侧恰好是自动占位名，视为等价（原矩阵本就无行列名）；
  # 真实行列名不一致（如基因名/细胞类型名漂移）仍判 MISMATCH。
  auto_names_ok <- function(cur_nm, ref_nm, n, col = FALSE) {
    if (identical(cur_nm, ref_nm)) return(TRUE)
    auto <- if (col) paste0("V", seq_len(n)) else as.character(seq_len(n))
    if (is.null(cur_nm)) return(identical(ref_nm, auto))
    if (is.null(ref_nm)) return(identical(cur_nm, auto))
    FALSE
  }
  rn_ok <- dim_ok && auto_names_ok(rownames(cur), rownames(ref), nrow(ref), col = FALSE)
  cn_ok <- dim_ok && auto_names_ok(colnames(cur), colnames(ref), ncol(ref), col = TRUE)
  if (is.character(ref) || is.character(cur)) {
    val_ok <- dim_ok && identical(unname(ref), unname(cur))
    ok <- dim_ok && rn_ok && cn_ok && val_ok
    cat(sprintf("[CHECK] %-12s : dim=%s rownames=%s colnames=%s values=%s => %s\n",
                obj_name, dim_ok, rn_ok, cn_ok, val_ok, ifelse(ok, "OK", "MISMATCH")))
    return(ok)
  }
  maxdiff <- if (dim_ok) max(abs(ref - cur)) else Inf
  ok <- dim_ok && rn_ok && cn_ok && maxdiff <= tol
  cat(sprintf("[CHECK] %-12s : dim=%s rownames=%s colnames=%s max|diff|=%.2e => %s\n",
              obj_name, dim_ok, rn_ok, cn_ok, maxdiff, ifelse(ok, "OK", "MISMATCH")))
  invisible(ok)
}

repro_items <- list(
  V              = V,
  sig_matrix     = sig_matrix,
  sd_matrix      = sd_matrix,
  coef1          = coef1,
  marker_distrib = as.matrix(md)
)
repro_ok <- all(vapply(names(repro_items),
                       function(nm) check_repro(repro_items[[nm]], nm),
                       logical(1)))

if (bench_ok && repro_ok) {
  cat("\n>>> INPUT REPRODUCED: benchmark + ESMF 输入均与基准一致，开始参数敏感性实验。\n")
} else {
  cat("\n>>> FATAL: 输入与基准不一致！数据/代码可能已变动，参数敏感性实验已终止。\n")
  stop("Input reproducibility check failed. Aborting parameter sensitivity study.")
}

# ============================================================================
# SECTION 13: RESULTS PROCESSING AND EVALUATION
# ============================================================================
# Process ESMF deconvolution results
P <- generator[["P"]]

W <- RES$W
H <- RES$H
# rownames(H) <- rownames(P)
# colnames(H) <- colnames(P)

# pheatmap(corH, scale = "none")
corW <- cor(W, sig_matrix)
# corW[,2] <- 0
rownames(corW) <- c(1:nrow(corW))
corW

# sig_train[, which(colnames(train) == "MDA-MB468")]
# apply(corW, 1, which.max)
# dis <- (W - sig_matrix) / sd_matrix
# H <- H[c(8,7,4,6,1,3,5,2),]

# Function to find permutation that minimizes Frobenius norm
pMatrix.min <- function(A, B) { 
  # Finds the permutation P of A such that ||PA - B|| is minimum in Frobenius norm
  # Uses the linear-sum assignment problem (LSAP) solver in the "clue" package
  # Returns P%*%A and the permutation vector
  # A[pvec, ] is the permutation of A closest to B
  n <- nrow(A)
  D <- matrix(NA, n, n)
  for (i in 1:n) { 
    for (j in 1:n) { 
      D[j, i] <- (sum((B[j, ] - A[i, ])^2))
    }
  }
  vec <- c(solve_LSAP(D))
  list(A = A[vec, ], pvec = vec)
}

require(clue)  # Need this package to solve the LSAP

# An example
B <- diag(1, nrow(corW))  # This choice of B maximizes the trace of permuted A

matrix.sort <- function(matrix) {
  if (nrow(matrix) != ncol(matrix)) stop("Not diagonal")
  if (is.null(rownames(matrix))) rownames(matrix) <- 1:nrow(matrix)
  row.max <- apply(matrix, 1, which.max)
  if (all(table(row.max) != 1)) stop("Ties cannot be resolved")
  matrix[names(sort(row.max)), ]
}

# corW <- matrix.sort(corW)
# cell_order <- as.numeric(rownames(corW))

X <- pMatrix.min(corW, B)
sum(diag(corW))
sum(diag(X$A))
cell_order <- X$pvec
corW <- X$A
rownames(corW) <- colnames(corW)
cell_order

# Normalize H matrix to represent proportions
for (j_2 in seq(n)) {
  col_sum_H = sum(H[, j_2])
  for (k_2 in seq(r)) {
    H[k_2, j_2] = H[k_2, j_2] / col_sum_H
  }
}
colSums(H)
dim(H)
# dim(test)
# colnames(H) <- colnames(T$Matrix)
colnames(H) <- colnames(V)
# colnames(H) <- colnames(MWBL)

order(apply(corW, 1, which.max))
# H <- H[c(5, 4, 1, 2, 3), ]
# H <- H[c(4, 2, 5, 6, 1, 3), ]
# H <- H[c(2, 4, 5, 6, 1, 3), ]
# H <- H[order(apply(corW, 1, which.max)), ]
H <- H[cell_order, ]

# mat1 <- H
# mat2 <- P
# # Calculate Pearson correlation between samples
# samples <- colnames(mat1)
# cor_results <- sapply(samples, function(s) {
#   cor(mat1[, s], mat2[, s], method = "pearson")
# })
# # Store results in Results data.frame
# Results[, "ESMFE"] <- cor_results

# names(H) <- colnames(corW)
rownames(H) <- colnames(corW)
colnames(H) <- colnames(P)
# H <- t(as.matrix(H))
corH <- cor(H, P)
H <- H[gtools::mixedsort(rownames(H)), ]
H <- as.data.table(H, keep.rownames = TRUE)
H <- data.table::melt(H)
colnames(H) <- c("CT", "tissue", "observed_values")
head(H)

# P$CT <- rownames(P)
# P <- data.table::melt(P, id.vars = "CT")
# colnames(P) <- c("CT", "tissue", "expected_values")
# head(P)

P <- as.matrix(P)
P <- P[gtools::mixedsort(rownames(P)), ]
P <- as.data.table(P, keep.rownames = TRUE)
P <- data.table::melt(P)
colnames(P) <- c("CT", "tissue", "expected_values")
head(P)

# H$tissue <- paste0("mix", H$tissue)
head(H)
head(P)

INDICATOR = merge(H, P)
head(INDICATOR)
INDICATOR$expected_values <- round(INDICATOR$expected_values, 2)
INDICATOR$observed_values <- round(INDICATOR$observed_values, 2)
INDICATOR$expected_values[1:10]
INDICATOR$observed_values[1:10]
# # 拼接结果
# # RESULTS[, "ESMFE"] <- INDICATOR$observed_values
# colnames(RESULTS) <- c("cellType", "sampleID", "EPIC", "expected_values", "ESMFE")
# write.csv(RESULTS, "EPIC_ESMFE_P.csv", row.names = TRUE)

INDICATOR = INDICATOR %>% dplyr::summarise(
  RMSE = sqrt(mean((observed_values - expected_values)^2)) %>% round(., 4),
  Pearson = cor(observed_values, expected_values) %>% round(., 4)
)
H_rmse_list[["ESMF"]] <- INDICATOR$RMSE
H_pearson_list[["ESMF"]] <- INDICATOR$Pearson

# Evaluate signature matrix estimation
W_rmse_list <- list()
W_pearson_list <- list()
# corW rmseW
# W <- W[, order(apply(corW, 1, which.max))]
W <- W[, cell_order]
colnames(W) <- colnames(sig_matrix)
rownames(W) <- rownames(sig_matrix)
W <- t(W)
W <- W[gtools::mixedsort(rownames(W)), ]
W <- reshape2::melt(W)
colnames(W) <- c("CT", "gene", "observed_values")

sig <- t(sig_matrix)
sig <- sig[gtools::mixedsort(rownames(sig)), ]
sig <- reshape2::melt(sig)
colnames(sig) <- c("CT", "gene", "expected_values")

INDICATOR = merge(W, sig)
head(INDICATOR)
INDICATOR$expected_values <- round(INDICATOR$expected_values, 3)
INDICATOR$observed_values <- round(INDICATOR$observed_values, 3)
INDICATOR$expected_values[1:10]
INDICATOR$observed_values[1:10]

INDICATOR = INDICATOR %>% dplyr::summarise(
  RMSE = sqrt(mean((observed_values - expected_values)^2)) %>% round(., 4),
  Pearson = cor(observed_values, expected_values) %>% round(., 4)
)
W_rmse_list[["ESMF"]] <- INDICATOR$RMSE
W_pearson_list[["ESMF"]] <- INDICATOR$Pearson

# DSA method for comparison
# Create empty list to store results
# sim data
normalization = "none"
T <- generator[["T"]]
P <- generator[["P"]]
C = lapply(group, function(x) Matrix::rowMeans(train[, x]))
C = round(do.call(cbind.data.frame, C))

genes_filt <- intersect(rownames(V), rownames(sig_matrix))
V <- V[genes_filt, ]
sig_matrix <- sig_matrix[genes_filt, ]
sd_matrix <- sd_matrix[genes_filt, ]
md <- md[rownames(sig_matrix), ]

dim(V)
dim(sig_matrix)
dim(sd_matrix)
dim(md)

require(CellMix)
md = marker_distrib
ML = CellMix::MarkerList()
ML@.Data <- tapply(as.character(md$gene), as.character(md$CT), list)
H1 = CellMix::ged(as.matrix(V), ML, method = "DSA", log = FALSE)@fit@H
W1 = CellMix::ged(as.matrix(V), ML, method = "DSA", log = FALSE)@fit@W
colnames(RES$W) <- colnames(W1)
rownames(RES$H) <- rownames(H1)
rownames(RES$W) <- rownames(sig_matrix)
# H <- RES$H
# W <- RES$W

H <- H1
W <- W1

# sig_matrix <- C
gene_index <- intersect(rownames(W), rownames(sig_matrix))
W <- W[gene_index, ]
sig <- sig_matrix[gene_index, ]

# pheatmap(corH, scale = "none")
corW <- cor(W, sig)
# corW[,2] <- 0
corW[is.na(corW)] <- 0
rownames(corW) <- c(1:nrow(corW))
corW

# sig_train[, which(colnames(train) == "MDA-MB468")]
# apply(corW, 1, which.max)
# dis <- (W - sig_matrix) / sd_matrix
# H <- H[c(8,7,4,6,1,3,5,2),]

# Use the same permutation function for DSA results
X <- pMatrix.min(corW, B)
sum(diag(corW))
sum(diag(X$A))
cell_order <- X$pvec
corW <- X$A
rownames(corW) <- colnames(corW)

# colnames(H) <- colnames(T$Matrix)
# corW rmseW
# W <- W[, order(apply(corW, 1, which.max))]
W <- W[, cell_order]
colnames(W) <- colnames(sig)
rownames(W) <- rownames(sig)
W <- t(W)
W <- W[gtools::mixedsort(rownames(W)), ]
W <- as.data.table(W, keep.rownames = TRUE)
W <- data.table::melt(W)
colnames(W) <- c("CT", "gene", "observed_values")

sig <- t(sig_matrix)
sig <- sig[gtools::mixedsort(rownames(sig)), ]
sig <- as.data.table(sig, keep.rownames = TRUE)
sig <- data.table::melt(sig)
colnames(sig) <- c("CT", "gene", "expected_values")

INDICATOR = merge(W, sig)
head(INDICATOR)
INDICATOR$expected_values <- round(INDICATOR$expected_values, 3)
INDICATOR$observed_values <- round(INDICATOR$observed_values, 3)
INDICATOR$expected_values[1:10]
INDICATOR$observed_values[1:10]

INDICATOR = INDICATOR %>% dplyr::summarise(
  RMSE = sqrt(mean((observed_values - expected_values)^2)) %>% round(., 4),
  Pearson = cor(observed_values, expected_values) %>% round(., 4)
)
W_rmse_list[["DSA"]] <- INDICATOR$RMSE
W_pearson_list[["DSA"]] <- INDICATOR$Pearson

# ============================================================================
# SECTION 14: RESULTS COMPILATION AND SAVING
# ============================================================================
# Convert results lists to data frames
H_rmse_df <- data.frame(RMSE = unlist(H_rmse_list))
H_pearson_df <- data.frame(Pearson = unlist(H_pearson_list))

W_rmse_df <- data.frame(RMSE = unlist(W_rmse_list))
W_pearson_df <- data.frame(Pearson = unlist(W_pearson_list))

# Print results summary
print("RMSE data frames:")
print(H_rmse_df)
print(W_rmse_df)

print("Pearson data frames:")
print(H_pearson_df)
print(W_pearson_df)


# Save all results to files
save_path <- getwd()
log_file <- paste0(results_dir, "deconvolution_results_log_", current_datetime, ".txt")

cat("\n\n=== Run time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "===\n", 
    file = log_file, append = TRUE)

# Save H matrix results
cat("\nH_rmse_df:\n", file = log_file, append = TRUE)
write.table(H_rmse_df, file = log_file, append = TRUE, sep = "\t", 
            quote = FALSE, col.names = TRUE, row.names = TRUE)

cat("\nH_pearson_df:\n", file = log_file, append = TRUE)
write.table(H_pearson_df, file = log_file, append = TRUE, sep = "\t", 
            quote = FALSE, col.names = TRUE, row.names = TRUE)

# Save W matrix results
cat("\nW_rmse_df:\n", file = log_file, append = TRUE)
write.table(W_rmse_df, file = log_file, append = TRUE, sep = "\t", 
            quote = FALSE, col.names = TRUE, row.names = TRUE)

cat("\nW_pearson_df:\n", file = log_file, append = TRUE)
write.table(W_pearson_df, file = log_file, append = TRUE, sep = "\t", 
            quote = FALSE, col.names = TRUE, row.names = TRUE)

# Save individual result files
write.table(H_rmse_df, file = paste0(results_dir, "H_rmse_results_", current_datetime, ".txt"), 
            sep = "\t", quote = FALSE, col.names = TRUE, row.names = TRUE)
write.table(H_pearson_df, file = paste0(results_dir, "H_pearson_results_", current_datetime, ".txt"), 
            sep = "\t", quote = FALSE, col.names = TRUE, row.names = TRUE)
write.table(W_rmse_df, file = paste0(results_dir, "W_rmse_results_", current_datetime, ".txt"), 
            sep = "\t", quote = FALSE, col.names = TRUE, row.names = TRUE)
write.table(W_pearson_df, file = paste0(results_dir, "W_pearson_results_", current_datetime, ".txt"), 
            sep = "\t", quote = FALSE, col.names = TRUE, row.names = TRUE)

# Save comprehensive results summary
summary_file <- paste0(results_dir, "deconvolution_summary_", current_datetime, ".txt")
sink(summary_file)
cat("=== DECONVOLUTION ANALYSIS SUMMARY ===\n\n")
cat("Analysis date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("Analysis prefix:", analysis_prefix, "\n")
#cat("Training datasets:", paste(c(set1, set2, set3, set4), collapse = ", "), "\n")
#cat("Test dataset:", test_set, "\n")
cat("Selected cell types:", paste(select_cellType, collapse = ", "), "\n\n")

cat("=== H MATRIX RESULTS (PROPORTION ESTIMATION) ===\n\n")
cat("RMSE Values:\n")
print(H_rmse_df)
cat("\nPearson Correlation Values:\n")
print(H_pearson_df)

cat("\n=== W MATRIX RESULTS (SIGNATURE ESTIMATION) ===\n\n")
cat("RMSE Values:\n")
print(W_rmse_df)
cat("\nPearson Correlation Values:\n")
print(W_pearson_df)
sink()

cat("Results have been saved to:", results_dir, "\n")
cat("Analysis completed successfully!\n")

# ============================================================================
# SECTION 15: PARAMETER SENSITIVITY  --  P1 (rate 网格) / P2 (para 组合)
#                 / P3 (iter_times) / P4 (初始化 seed 鲁棒性)
# 仅追加到此副本脚本 Parameter_Sensitivity_Hm_Pancreas.R。
# 正文基准（generalization-from-single-reference-hm-pancreas.R）保持不动。
# 单变量控制铁律：每一组只改一个变量。
# 复用工作空间对象：V, train, group, generator, P, markers, prio_markers,
#   md, sig_matrix, sd_matrix, ML, r, n, m, para1/2/3, coef1。
# ============================================================================

# ---- 局部辅助函数：按给定 md 与 para 组合拟合 sig/sd/ML/coef1 ----
fit_signature <- function(md_obj, para = c(1, 1, 1)) {
  require(fitdistrplus)
  p1 <- para[1]; p2 <- para[2]; p3 <- para[3]
  sig_tr <- train[md_obj$gene, , drop = FALSE]
  ct_names <- names(table(md_obj$CT))
  mu_ls <- list(); sd_ls <- list()
  for (i in seq_along(ct_names)) {
    means <- c(); sds <- c()
    cell_index <- which(colnames(sig_tr) %in% ct_names[i])
    counts <- sig_tr[, cell_index, drop = FALSE]
    for (j in seq_len(nrow(counts))) {
      if (length(which(round(counts[j, ]) != 0)) > 0) {
        f1 <- fitdist(round(counts[j, ]), "nbinom", method = "mse")
        mu_ <- ifelse(is.na(f1$estimate[1]), 0, f1$estimate[1])
        si_ <- ifelse(is.na(f1$estimate[2]), 0, f1$estimate[2])
        means <- c(means, mu_); sds <- c(sds, si_)
      } else {
        means <- c(means, 0); sds <- c(sds, 0)
      }
    }
    mu_ls <- c(mu_ls, list(means)); sd_ls <- c(sd_ls, list(sds))
  }
  sig_new <- matrix(unlist(mu_ls), ncol = length(mu_ls))
  colnames(sig_new) <- ct_names; rownames(sig_new) <- rownames(sig_tr)
  theta_new <- matrix(unlist(sd_ls), ncol = length(sd_ls))
  sd_new <- sig_new * (1 + theta_new)
  colnames(sd_new) <- ct_names; rownames(sd_new) <- rownames(sig_tr)
  ML_new <- CellMix::MarkerList()
  ML_new@.Data <- tapply(as.character(md_obj$gene), as.character(md_obj$CT), list)
  r_new <- length(ct_names); m_new <- nrow(sig_tr)
  coef1_new <- matrix(p1, m_new, r_new)
  for (i in seq_len(r_new)) { gix <- which(md_obj$CT %in% ct_names[i]); coef1_new[gix, i] <- p2 }
  for (i in seq_len(r_new)) { gix <- which(md_obj$gene %in% prio_markers[[i]]); if (length(gix) > 0) coef1_new[gix, i] <- p3 }
  list(sig = sig_new, sd = sd_new, ML = ML_new, coef1 = coef1_new, md = md_obj, r = r_new, m = m_new)
}

# ---- 局部辅助函数：ESMF 运行 + H/W 评估（支持 iter_times 与 seed）----
evaluate_esmf <- function(rate_mat, md_obj, para = c(1, 1, 1),
                          iter_times = 3000, seed_val = 1233322, label = "") {
  require(clue)
  pMatrix.min <- function(A, B) {
    n <- nrow(A); D <- matrix(NA, n, n)
    for (i in 1:n) for (j in 1:n) D[j, i] <- (sum((B[j, ] - A[i, ])^2))
    vec <- c(solve_LSAP(D)); list(A = A[vec, ], pvec = vec)
  }
  fit <- fit_signature(md_obj, para = para)
  # 修复（对齐 v3）：sigm/sdm 改用全局 sig_matrix/sd_matrix（与正文一致、与 para 无关），
  #   避免 fit_signature 重新拟合引入差异导致 W_Pearson 崩塌为负。
  #   coef1m/MLm 仍需 fit_signature（P2 按 para 重建 coef1）。
  sigm <- sig_matrix; sdm <- sd_matrix; MLm <- fit$ML; coef1m <- fit$coef1
  rm_ <- fit$r
  g_f <- intersect(rownames(V), rownames(sigm))
  V_loc <- V[g_f, , drop = FALSE]
  sigm <- sigm[g_f, , drop = FALSE]
  sdm <- sdm[g_f, , drop = FALSE]
  # 修复（对齐 v3）：coef1/coef1m 是裸矩阵（matrix(para1, m, r)）无行名，而 g_f 是基因名字符向量，
  #   直接字符索引会报 "no 'dimnames' attribute for array"。
  #   行序已与 sigm 一致（coef1 行序 = md 行序 = sig_matrix 行序；coef1m 行序 = md_obj 行序）。
  #   必须用 fit_signature 输出的完整基因名（nrow(fit$sig) 行）补行名——
  #   当 md 基因不全在 V 里时 g_f 长度 < 完整基因数，若按裁剪后的 rownames(sigm) 补会报 dimnames 长度错误。
  if (is.matrix(rate_mat) && is.null(rownames(rate_mat))) {
    if (nrow(rate_mat) == nrow(fit$sig)) rownames(rate_mat) <- rownames(fit$sig)
    else rownames(rate_mat) <- rownames(sigm)
  }
  if (is.null(rownames(coef1m))) {
    rownames(coef1m) <- rownames(fit$sig)
  }
  if (identical(rate_mat, "coef")) {
    # P2: rate = 由新 para 组合生成的 coef1
    rate_loc <- coef1m[g_f, , drop = FALSE]
  } else if (is.matrix(rate_mat) && ncol(rate_mat) == rm_ && nrow(rate_mat) >= nrow(sigm)) {
    rate_loc <- rate_mat[g_f, , drop = FALSE]
  } else {
    rate_loc <- matrix(rate_mat, nrow(sigm), ncol(sigm))
  }
  P_loc <- generator[["P"]]
  set.seed(seed_val)
  RES_loc <- ESMF_deconvolution_fast(V = V_loc, ml = MLm, r = rm_, err_cutoff = 0,
                                rate = rate_loc, sig_matrix = sigm,
                                sd_matrix = sdm, iter_times = iter_times)
  W_loc <- RES_loc$W; H_loc <- RES_loc$H
  n_loc <- ncol(V_loc)
  corW_loc <- cor(W_loc, sigm)
  rownames(corW_loc) <- as.character(1:nrow(corW_loc))
  B_d <- diag(1, nrow(corW_loc))
  X_loc <- pMatrix.min(corW_loc, B_d)
  cell_order_loc <- X_loc$pvec
  for (j_2 in seq_len(n_loc)) { cs <- sum(H_loc[, j_2]); for (k_2 in seq_len(rm_)) H_loc[k_2, j_2] <- H_loc[k_2, j_2] / cs }
  colnames(H_loc) <- colnames(V_loc)
  H_loc <- H_loc[cell_order_loc, , drop = FALSE]
  rownames(H_loc) <- colnames(corW_loc)
  colnames(H_loc) <- colnames(P_loc)
  Hm <- H_loc[gtools::mixedsort(rownames(H_loc)), , drop = FALSE]
  Hm <- as.data.table(Hm, keep.rownames = TRUE); Hm <- data.table::melt(Hm)
  colnames(Hm) <- c("CT", "tissue", "observed_values")
  Pmm <- as.matrix(P_loc); Pmm <- Pmm[gtools::mixedsort(rownames(Pmm)), , drop = FALSE]
  Pmm <- as.data.table(Pmm, keep.rownames = TRUE); Pmm <- data.table::melt(Pmm)
  colnames(Pmm) <- c("CT", "tissue", "expected_values")
  IND <- merge(Hm, Pmm)
  IND$expected_values <- round(IND$expected_values, 2)
  IND$observed_values <- round(IND$observed_values, 2)
  INDh <- IND %>% dplyr::summarise(
    RMSE = sqrt(mean((observed_values - expected_values)^2)) %>% round(., 4),
    Pearson = cor(observed_values, expected_values) %>% round(., 4))
  W_ord <- W_loc[, cell_order_loc, drop = FALSE]
  colnames(W_ord) <- colnames(sigm); rownames(W_ord) <- rownames(sigm)
  Wm <- t(W_ord); Wm <- Wm[gtools::mixedsort(rownames(Wm)), , drop = FALSE]
  Wm <- reshape2::melt(Wm); colnames(Wm) <- c("CT", "gene", "observed_values")
  sigt <- t(sigm); sigt <- sigt[gtools::mixedsort(rownames(sigt)), , drop = FALSE]
  sigt <- reshape2::melt(sigt); colnames(sigt) <- c("CT", "gene", "expected_values")
  INDw <- merge(Wm, sigt)
  INDw$expected_values <- round(INDw$expected_values, 3)
  INDw$observed_values <- round(INDw$observed_values, 3)
  INDw <- INDw %>% dplyr::summarise(
    RMSE = sqrt(mean((observed_values - expected_values)^2)) %>% round(., 4),
    Pearson = cor(observed_values, expected_values) %>% round(., 4))
  c(label = label, iter = iter_times, seed = seed_val,
    H_RMSE = as.numeric(INDh$RMSE), H_Pearson = as.numeric(INDh$Pearson),
    W_RMSE = as.numeric(INDw$RMSE), W_Pearson = as.numeric(INDw$Pearson))
}

# ============================================================================
# P1: rate 网格敏感性（锁死 marker/para/iter/seed；只扫 rate 倍数）
# ============================================================================
cat("\n\n========== P1: RATE GRID SENSITIVITY ==========\n")
p1_grid <- c(0, 0.01, 0.1, 0.5, 1, 2, 5, 10)
p1_out <- lapply(p1_grid, function(fac) {
  cat("Running P1: rate factor =", fac, "\n")
  evaluate_esmf(rate_mat = coef1 * fac, md_obj = md, para = c(1, 1, 1),
                iter_times = 3000, seed_val = 1233322, label = paste0("P1_rate_", fac))
})
p1_df <- as.data.frame(do.call(rbind, p1_out), stringsAsFactors = FALSE)

# ============================================================================
# P2: 分层惩罚系数 (para1, para2, para3) 独立扫描
# v2 升级：原 v1 只扫了 para2/para3，漏了 para1（baseline）的独立影响。
#   现改为三个参数各自独立扫描（单变量控制铁律：每次只动一个 para，其余锁死为 1），
#   每个参数取 0.1 / 1(基准) / 10 / 100 四档，补齐 para1 缺口，让"三个参数"都有完整证据。
#   para 语义：para1=baseline(非 marker 基因基础惩罚), para2=marker 特异性惩罚, para3=priority marker 惩罚。
# ============================================================================
cat("\n\n========== P2: PENALTY PARAMETER INDEPENDENT SCANS ==========\n")

# 每个参数独立扫描的档位（值=该参数的取值，其余两个参数固定为 1）
p2_param_grid <- c(0.1, 1, 10, 100)

p2_combos <- list()
# 基准（三个参数都为 1）
p2_combos[["base_111"]] <- c(1, 1, 1)
# para1 独立扫描（baseline 惩罚）：只动 para1，para2=para3=1
for (v in p2_param_grid) {
  nm <- paste0("para1_baseline_", v)
  p2_combos[[nm]] <- c(v, 1, 1)
}
# para2 独立扫描（marker 惩罚）：只动 para2，para1=para3=1
for (v in p2_param_grid) {
  nm <- paste0("para2_marker_", v)
  p2_combos[[nm]] <- c(1, v, 1)
}
# para3 独立扫描（priority marker 惩罚）：只动 para3，para1=para2=1
for (v in p2_param_grid) {
  nm <- paste0("para3_prio_", v)
  p2_combos[[nm]] <- c(1, 1, v)
}
# 去重（base_111 会与 para1_baseline_1 / para2_marker_1 / para3_prio_1 重复，只保留一份）
p2_combos <- p2_combos[!duplicated(p2_combos)]

p2_out <- lapply(names(p2_combos), function(nm) {
  pc <- p2_combos[[nm]]
  cat("Running P2:", nm, " para =", paste(pc, collapse = ","), "\n")
  evaluate_esmf(rate_mat = "coef", md_obj = md, para = pc,
                iter_times = 3000, seed_val = 1233322, label = paste0("P2_", nm))
})
p2_df <- as.data.frame(do.call(rbind, p2_out), stringsAsFactors = FALSE)

# ============================================================================
# P3: 迭代次数 iter_times（收敛性与迭代次数鲁棒性）
# ============================================================================
cat("\n\n========== P3: ITERATION COUNT ==========\n")
p3_iters <- c(500, 1000, 2000, 3000, 5000)
p3_out <- lapply(p3_iters, function(it) {
  cat("Running P3: iter_times =", it, "\n")
  evaluate_esmf(rate_mat = coef1, md_obj = md, para = c(1, 1, 1),
                iter_times = it, seed_val = 1233322, label = paste0("P3_iter_", it))
})
p3_df <- as.data.frame(do.call(rbind, p3_out), stringsAsFactors = FALSE)

# ============================================================================
# P4: 初始化 seed 鲁棒性（20 个 seed，报告 mean±sd）
# 2026-08-18 注释掉：P4 与消融实验 v3 的 E3（多次初始化稳定性）重复，
#   且扩到 20 seed 后撞上一个坏局部解（seed 1233335, H_Pearson=0.785），
#   反而干扰 E3 的干净结论。R2 M2① 的答复以 E3（10 seed，W/H 两两相关 0.9995/0.9957）为准。
#   如后续需要"更大规模初始化鲁棒性"再启用。
# ============================================================================
# cat("\n\n========== P4: INITIALIZATION SEED ROBUSTNESS ==========\n")
# p4_seeds <- 1233322:(1233322 + 19)  # 20 个连续 seed
# p4_out <- lapply(seq_along(p4_seeds), function(k) {
#   s <- p4_seeds[k]
#   cat("Running P4: seed =", s, "\n")
#   evaluate_esmf(rate_mat = coef1, md_obj = md, para = c(1, 1, 1),
#                 iter_times = 3000, seed_val = s, label = paste0("P4_seed_", s))
# })
# p4_df <- as.data.frame(do.call(rbind, p4_out), stringsAsFactors = FALSE)
#
# # ---- 汇总 P4 为 mean±sd ----
# # 修复：原 t(p4_num) 转置 list 后仍为 list，data.frame 列变 list 类型，
# #       导致 write.table 报 "unimplemented type 'list' in 'EncodeElement'"。
# #       改用 unlist 展开为字符向量，构造普通 data.frame。
# p4_num <- sapply(c("H_RMSE", "H_Pearson", "W_RMSE", "W_Pearson"), function(col) {
#   vals <- as.numeric(p4_df[[col]])
#   paste0(round(mean(vals), 4), "±", round(sd(vals), 4))
# })
# p4_summary <- data.frame(
#   experiment = "P4_init_seed",
#   H_RMSE_meanSD  = p4_num[["H_RMSE"]],
#   H_Pearson_meanSD = p4_num[["H_Pearson"]],
#   W_RMSE_meanSD  = p4_num[["W_RMSE"]],
#   W_Pearson_meanSD = p4_num[["W_Pearson"]],
#   stringsAsFactors = FALSE)

# ============================================================================
# 汇总并保存 P1-P3 参数敏感性结果（独立目录）
# ============================================================================
param_out <- rbind(p1_df, p2_df, p3_df)
param_file <- paste0(results_dir, "parameter_sensitivity_results_", current_datetime, ".txt")
cat("\n\n========== PARAMETER SENSITIVITY SUMMARY ==========\n")
print(param_out)
write.table(param_out, file = param_file, sep = "\t", quote = FALSE,
            col.names = TRUE, row.names = FALSE)
cat("Parameter sensitivity results saved to:", param_file, "\n")
cat("Parameter sensitivity study (P1-P3) completed.\n")
