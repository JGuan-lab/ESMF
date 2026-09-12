# ============================================================================
# COMPREHENSIVE DECONVOLUTION ANALYSIS FOR MOUSE BRAIN SINGLE-CELL DATA
# Multi-Source Reference Integration and Method Comparison
# ============================================================================
# DESCRIPTION:
# This script performs a comprehensive benchmarking of multiple deconvolution
# methods for mouse brain single-cell RNA-seq data. It integrates data from
# multiple sources, performs quality control, marker gene selection, generates
# pseudo-bulk mixtures, and evaluates 12 deconvolution methods including both
# bulk and single-cell approaches. The analysis focuses on three brain cell types:
# Astrocytes, Endothelial cells, and Microglia.

# ANALYTICAL WORKFLOW:
# 1. Data loading and preprocessing from multiple PanglaoDB datasets
# 2. Quality control and data filtering
# 3. Training-test splitting and reference matrix construction
# 4. Marker gene selection using limma-voom
# 5. Generation of 1000 pseudo-bulk mixtures for validation
# 6. Evaluation of 10 bulk deconvolution methods
# 7. Evaluation of 2 single-cell deconvolution methods (MuSiC, CIBERSORTx)
# 8. Signature matrix estimation using negative binomial distribution
# 9. Non-negative matrix factorization deconvolution
# 10. Results compilation and saving

# OUTPUT:
# Results are saved to: /media/desk16/tjn050/LSC/LSC_deconvolution/LSC/results/MultiSourceReferenceIntegration_Mm_Brain/

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
analysis_prefix <- "crossPlatformConsistency_Hm_Adipose_threshold200"
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
# SECTION 3: LOAD AND PREPROCESS TRAINING DATA FROM HUSCH DATASET
# ============================================================================
# This section loads and preprocesses human adipose tissue single-cell RNA-seq data
# from the HUSCH database. The training dataset is HU_0303_Adipose_GSE153643.
# Selected cell types: CD4T, Fibroblast, Mono_Macro
# ============================================================================

# ----------------------------------------------------------------------------
# Training Dataset: HU_0303_Adipose_GSE153643
# ----------------------------------------------------------------------------
# Set training data directory and load metadata
train_dir = "/media/desk16/tjn050/LSC/HUSCH/Adipose/HU_0303_Adipose_GSE153643/"
train_set = strsplit(train_dir, "/")[[1]][8]
cat("Loading training dataset:", train_set, "\n")

# Load cell phenotype data (metadata)
phenoData_train <- read.csv(paste0(train_dir, train_set, "_meta.txt"), 
                            sep = "\t", header = TRUE) 

# Load sparse expression matrix from HDF5 file
train_mat <- h5read(paste0(train_dir, train_set, "_gene_count.h5"), "matrix")

# Clean cell type names: replace special characters for consistency
phenoData_train$Celltype <- sub("/", "_", phenoData_train$Celltype)
phenoData_train$Celltype <- sub(" ", "_", phenoData_train$Celltype)
phenoData_train$Celltype <- sub(" ", "_", phenoData_train$Celltype)

# Standardize column names
colnames(phenoData_train) = c("cellID", "cellType", "Platform")

# Display initial cell type distribution
cat("Initial cell type distribution in training data:\n")
print(table(phenoData_train$cellType))

# ----------------------------------------------------------------------------
# Process Sparse Expression Matrix
# ----------------------------------------------------------------------------
# Extract components from HDF5 sparse matrix structure
library(rhdf5)
barcodes <- train_mat$barcodes
counts <- train_mat$data
gene_id <- train_mat$features$id
indices <- train_mat$indices
indptr <- train_mat$indptr
shape <- train_mat$shape

# Reconstruct sparse matrix using Matrix package
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

# Convert to dense matrix for downstream analysis
train_data <- as.matrix(mat)
rm(mat)  # Remove temporary sparse matrix to free memory

# ----------------------------------------------------------------------------
# Standardize Sample Identifiers and Filter Cell Types
# ----------------------------------------------------------------------------
# Add "train" prefix to sample identifiers for consistency
colnames(train_data) <- paste0("train", colnames(train_data))
phenoData_train$cellID <- paste0("train", phenoData_train$cellID)

# Ensure phenotype and expression data samples match
phenoData_train <- phenoData_train[which(phenoData_train$cellID %in% colnames(train_data)), ]
train_data <- train_data[, intersect(phenoData_train$cellID, colnames(train_data))]

cat("Cell type distribution after initial matching:\n")
print(table(phenoData_train$cellType))

# ----------------------------------------------------------------------------
# Select Specific Cell Types for Analysis
# ----------------------------------------------------------------------------
# For this adipose tissue analysis, we focus on three major cell types:
# CD4T (CD4+ T cells), Fibroblast, and Mono_Macro (Monocytes/Macrophages)
select_cellType <- c("CD4T", "Fibroblast", "Mono_Macro")

# Filter training data to selected cell types
phenoData_train <- phenoData_train[which(phenoData_train$cellType %in% select_cellType), ]
phenoData_train <- phenoData_train[which(phenoData_train$cellID %in% colnames(train_data)), ]
train_data <- train_data[, intersect(phenoData_train$cellID, colnames(train_data))]

cat("Cell type distribution after filtering to selected types:\n")
print(table(phenoData_train$cellType))

# ----------------------------------------------------------------------------
# Finalize Training Data Structure
# ----------------------------------------------------------------------------
# Ensure matrix format and prepare phenotype data structure
train_data <- as(train_data, "matrix")
phenoData_train$name <- phenoData_train$cellID

# Reorder columns to match expected format: cellType, cellID, name
phenoData_train <- phenoData_train[, c(2, 1, 4)]

cat("Final training data dimensions:\n")
cat("Expression matrix:", dim(train_data)[1], "genes x", dim(train_data)[2], "cells\n")
cat("Phenotype data:", dim(phenoData_train)[1], "cells x", dim(phenoData_train)[2], "attributes\n")

# ============================================================================
# SECTION 4: LOAD AND PREPROCESS TEST DATA FROM HUSCH DATASET
# ============================================================================
# This section loads and preprocesses independent test data for validation
# Test dataset: HU_0273_Adipose_GSE134355
# Selected cell types: CD4T, Fibroblast, Mono_Macro (same as training)
# ============================================================================

# ----------------------------------------------------------------------------
# Test Dataset: HU_0273_Adipose_GSE134355
# ----------------------------------------------------------------------------
# Set test data directory and load metadata
test_dir = "/media/desk16/tjn050/LSC/HUSCH/Adipose/HU_0273_Adipose_GSE134355/"
test_set = strsplit(test_dir, "/")[[1]][8]
cat("\nLoading test dataset:", test_set, "\n")

# Load cell phenotype data for test set
phenoData_test <- read.csv(paste0(test_dir, test_set, "_meta.txt"), 
                           sep = "\t", header = TRUE)

# Load sparse expression matrix from HDF5 file
test_mat <- h5read(paste0(test_dir, test_set, "_gene_count.h5"), "matrix")

# Clean cell type names using same procedure as training data
phenoData_test$Celltype <- sub("/", "_", phenoData_test$Celltype)
phenoData_test$Celltype <- sub(" ", "_", phenoData_test$Celltype)
phenoData_test$Celltype <- sub(" ", "_", phenoData_test$Celltype)

# Standardize column names
colnames(phenoData_test) = c("cellID", "cellType", "Platform")

cat("Initial cell type distribution in test data:\n")
print(table(phenoData_test$cellType))

# ----------------------------------------------------------------------------
# Process Sparse Expression Matrix for Test Data
# ----------------------------------------------------------------------------
# Extract components from HDF5 sparse matrix structure
barcodes <- test_mat$barcodes
counts <- test_mat$data
gene_id <- test_mat$features$id
indices <- test_mat$indices
indptr <- test_mat$indptr
shape <- test_mat$shape

# Reconstruct sparse matrix
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
test_data <- as.matrix(mat)
rm(mat)  # Remove temporary sparse matrix

# ----------------------------------------------------------------------------
# Standardize Sample Identifiers for Test Data
# ----------------------------------------------------------------------------
# Add "test" prefix to sample identifiers for consistency
colnames(test_data) <- paste0("test", colnames(test_data))
phenoData_test$cellID <- paste0("test", phenoData_test$cellID)

# Ensure phenotype and expression data samples match
phenoData_test <- phenoData_test[which(phenoData_test$cellID %in% colnames(test_data)), ]
test_data <- test_data[, intersect(phenoData_test$cellID, colnames(test_data))]

cat("Cell type distribution after initial matching in test data:\n")
print(table(phenoData_test$cellType))

# ----------------------------------------------------------------------------
# Filter Test Data to Selected Cell Types
# ----------------------------------------------------------------------------
# Use the same cell type selection as training data
phenoData_test <- phenoData_test[which(phenoData_test$cellType %in% select_cellType), ]
phenoData_test <- phenoData_test[which(phenoData_test$cellID %in% colnames(test_data)), ]
test_data <- test_data[, intersect(phenoData_test$cellID, colnames(test_data))]

cat("Cell type distribution after filtering to selected types in test data:\n")
print(table(phenoData_test$cellType))

# ----------------------------------------------------------------------------
# Finalize Test Data Structure
# ----------------------------------------------------------------------------
# Ensure matrix format and prepare phenotype data structure
test_data <- as(test_data, "matrix")
phenoData_test$name <- phenoData_test$cellID

# Reorder columns to match training data format: cellType, cellID, name
phenoData_test <- phenoData_test[, c(2, 1, 4)]

cat("Final test data dimensions:\n")
cat("Expression matrix:", dim(test_data)[1], "genes x", dim(test_data)[2], "cells\n")
cat("Phenotype data:", dim(phenoData_test)[1], "cells x", dim(phenoData_test)[2], "attributes\n")

# ============================================================================
# DATA INTEGRATION AND GENE NAME STANDARDIZATION
# ============================================================================
# This section integrates training and test data, ensuring gene name consistency
# for downstream deconvolution analysis
# ============================================================================

# ----------------------------------------------------------------------------
# Standardize Gene Names
# ----------------------------------------------------------------------------
# Extract primary gene symbol (part before underscore) and convert to uppercase
rownames(train_data) <- unlist(lapply(strsplit(rownames(train_data), "_"), function(x) x[1]))
rownames(test_data) <- unlist(lapply(strsplit(rownames(test_data), "_"), function(x) x[1]))

# Convert all gene symbols to uppercase for case-insensitive matching
rownames(train_data) <- toupper(rownames(train_data))
rownames(test_data) <- toupper(rownames(test_data))

# Keep only genes present in both training and test datasets
common_genes <- intersect(rownames(train_data), rownames(test_data))
train_data <- train_data[common_genes, ]
test_data <- test_data[common_genes, ]

cat("\nGene intersection between training and test data:\n")
cat(length(common_genes), "common genes retained\n")
cat("Training data dimensions after intersection:", dim(train_data)[1], "genes x", dim(train_data)[2], "cells\n")
cat("Test data dimensions after intersection:", dim(test_data)[1], "genes x", dim(test_data)[2], "cells\n")

# ----------------------------------------------------------------------------
# Combine Training and Test Data
# ----------------------------------------------------------------------------
# Create combined expression matrix and phenotype data for quality control
data <- as.matrix(cbind(train_data, test_data))
full_phenoData <- rbind(phenoData_train, phenoData_test)

cat("\nFinal combined dataset for deconvolution analysis:\n")
cat("Expression matrix (data):", dim(data)[1], "genes x", dim(data)[2], "cells\n")
cat("Phenotype data (full_phenoData):", dim(full_phenoData)[1], "cells x", dim(full_phenoData)[2], "attributes\n")

cat("\nCell type distribution in combined dataset:\n")
print(table(full_phenoData$cellType))

cat("\nTraining set cell type distribution:\n")
print(table(phenoData_train$cellType))

cat("\nTest set cell type distribution:\n")
print(table(phenoData_test$cellType))

# ============================================================================
# DATA INTEGRATION COMPLETE - READY FOR QUALITY CONTROL (SECTION 5)
# ============================================================================
# The variables 'data' and 'full_phenoData' are now prepared for quality control
# and subsequent deconvolution analysis. The data structure is consistent with
# the original pipeline expectations.
# ============================================================================

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
to_keep = names(cell_counts)[cell_counts >= 200]
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


#"HU_0303_Adipose_GSE153643"/"HU_0273_Adipose_GSE134355" sim22

prio_markers=list(
  CD4T=c("CCR7","CD28","CD3","CD4","IL7R","TCF7"),
  Fibroblast =c("ACTA2","COL1A1","COL1A2","COL3A1","COL6A1","COL6A2","COL6A3","DCN","FBLN1","FBLN2","FBN1","FN1","GSN","LY6A","MGP","MYL9","PDGFRA","PI16","POSTN"),
  Mono_Macro =c("CD14","CD68","CSF1R","CST3","FABP4","MMP9")
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
methods_to_use <- c("nnls", "FARDEEP", "RLR", "DCQ", "elastic_net", "lasso", "ridge", "OLS", "EPIC", "DSA")

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

methods_to_use <- c("nnls", "FARDEEP", "RLR", "DCQ", "elastic_net", "lasso", "ridge", "OLS", "DSA", "EPIC")

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
results <- read.table("/media/desk16/tjn050/LSC/LSC_deconvolution/LSC/results/CIBERSORTx_HU_0303_Adipose_GSE153643_HU_0273_Adipose_GSE134355_Results.csv", 
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
# SECTION 12: ESMF DECONVOLUTION
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


set.seed(1233322)


RES <- ESMF_deconvolution(V = V, ml = ML, r = r, err_cutoff = 0, rate = coef1, sig_matrix = sig_matrix, sd_matrix = sd_matrix, iter_times = 3000)

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
