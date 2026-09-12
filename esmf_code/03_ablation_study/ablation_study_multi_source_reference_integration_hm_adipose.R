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
analysis_prefix <- "Ablation_Study_MultiSourceReferenceIntegration_Hm_Adipose"
ablation_version <- "v3"
results_dir <- file.path("/media/desk16/tjn050/LSC/LSC_deconvolution/LSC/results",
                         analysis_prefix, ablation_version, "")
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
# SECTION 3: LOAD AND PREPROCESS TRAINING DATA FROM HUSCH ADIPOSE DATASETS
# ============================================================================
# This section loads and preprocesses four human adipose tissue single-cell
# RNA-seq datasets from the HUSCH database (GSE134355).
# Training datasets: HU_0271_Adipose, HU_0272_Adipose,
#                    HU_0273_Adipose, HU_0274_Adipose
# Test dataset: HU_0303_Adipose (GSE153643)
# Selected cell types: CD4T, Fibroblast, Mono_Macro
# Data format: H5 (gene count matrix) and TXT (metadata)
# ============================================================================

cat("Loading and integrating HUSCH human adipose tissue datasets...\n")

# ----------------------------------------------------------------------------
# Helper Function: Load a Single HUSCH Dataset
# ----------------------------------------------------------------------------
# Function to load expression matrix and metadata from HUSCH H5/TXT files
# Args:
#   data_dir: Directory containing the dataset
#   prefix: Unique prefix to assign to cell IDs (e.g., "train1")
# Returns: A list containing the expression matrix and phenotype data frame
load_husch_dataset <- function(data_dir, prefix) {
  # Extract dataset name from directory path
  dataset_name <- tail(strsplit(data_dir, "/")[[1]], 1)
  cat("  Processing dataset:", dataset_name, "\n")
  
  # Define file paths
  meta_file <- paste0(data_dir, dataset_name, "_meta.txt")
  h5_file <- paste0(data_dir, dataset_name, "_gene_count.h5")
  
  # Check if files exist
  if (!file.exists(meta_file)) stop("Metadata file not found: ", meta_file)
  if (!file.exists(h5_file)) stop("H5 file not found: ", h5_file)
  
  # --------------------------------------------------------------------
  # Step 1: Load and clean metadata
  # --------------------------------------------------------------------
  pheno_data <- read.csv(meta_file, sep = "\t", header = TRUE)
  
  # Standardize column names
  colnames(pheno_data) <- c("cellID", "cellType", "Platform")
  
  # Clean cell type names: replace special characters
  pheno_data$cellType <- gsub("/", "_", pheno_data$cellType)
  pheno_data$cellType <- gsub(" ", "_", pheno_data$cellType)
  
  # --------------------------------------------------------------------
  # Step 2: Load sparse expression matrix from H5 file
  # --------------------------------------------------------------------
  h5_data <- rhdf5::h5read(h5_file, "matrix")
  
  # Construct sparse matrix using Matrix package
  sparse_mat <- Matrix::sparseMatrix(
    i = h5_data$indices[] + 1,      # Convert from 0-based to 1-based indexing
    p = h5_data$indptr[],
    x = as.numeric(h5_data$data[]),
    dims = h5_data$shape[],
    repr = "T"
  )
  
  # Set row and column names
  rownames(sparse_mat) <- h5_data$features$id
  colnames(sparse_mat) <- h5_data$barcodes
  
  # Convert to dense matrix for compatibility with downstream analysis
  expr_mat <- as.matrix(sparse_mat)
  rm(sparse_mat, h5_data)  # Free memory
  gc()
  
  # --------------------------------------------------------------------
  # Step 3: Add prefix to cell identifiers and ensure consistency
  # --------------------------------------------------------------------
  colnames(expr_mat) <- paste0(prefix, colnames(expr_mat))
  pheno_data$cellID <- paste0(prefix, pheno_data$cellID)
  
  # Subset to cells present in both metadata and expression matrix
  common_cells <- intersect(pheno_data$cellID, colnames(expr_mat))
  pheno_data <- pheno_data[pheno_data$cellID %in% common_cells, ]
  expr_mat <- expr_mat[, common_cells, drop = FALSE]
  
  cat("    Loaded", ncol(expr_mat), "cells and", nrow(expr_mat), "genes.\n")
  cat("    Cell type distribution:\n")
  print(table(pheno_data$cellType))
  
  return(list(expr = expr_mat, pheno = pheno_data))
}

# ----------------------------------------------------------------------------
# Step 1: Load Individual Training Datasets
# ----------------------------------------------------------------------------
cat("\n--- Loading Training Datasets ---\n")

# Define paths to four training datasets
train_dirs <- c(
  "/media/desk16/tjn050/LSC/HUSCH/Adipose/HU_0271_Adipose_GSE134355/",
  "/media/desk16/tjn050/LSC/HUSCH/Adipose/HU_0272_Adipose_GSE134355/",
  "/media/desk16/tjn050/LSC/HUSCH/Adipose/HU_0273_Adipose_GSE134355/",
  "/media/desk16/tjn050/LSC/HUSCH/Adipose/HU_0274_Adipose_GSE134355/"
)

train_prefixes <- paste0("train", 1:4)

# Initialize lists to store data
train_data_list <- list()
pheno_data_list <- list()

# Load each dataset
for (i in seq_along(train_dirs)) {
  dataset <- load_husch_dataset(train_dirs[i], train_prefixes[i])
  train_data_list[[i]] <- dataset$expr
  pheno_data_list[[i]] <- dataset$pheno
}
names(train_data_list) <- train_prefixes
names(pheno_data_list) <- train_prefixes

# ----------------------------------------------------------------------------
# Step 2: Integrate Training Datasets Sequentially
# ----------------------------------------------------------------------------
cat("\n--- Integrating Training Datasets ---\n")

# Start with first two datasets
cat("1. Merging train1 and train2...\n")
train_data_temp <- train_data_list[[1]]
pheno_data_temp <- pheno_data_list[[1]]

# Convert gene names to uppercase for consistency
rownames(train_data_temp) <- toupper(rownames(train_data_temp))
rownames(train_data_list[[2]]) <- toupper(rownames(train_data_list[[2]]))

# Keep only common genes
common_genes <- intersect(rownames(train_data_temp), rownames(train_data_list[[2]]))
train_data_temp <- train_data_temp[common_genes, , drop = FALSE]
train_data_list[[2]] <- train_data_list[[2]][common_genes, , drop = FALSE]

# Merge expression matrices and phenotype data
train_data_temp <- cbind(train_data_temp, train_data_list[[2]])
pheno_data_temp <- rbind(pheno_data_temp, pheno_data_list[[2]])

cat("   After merge: ", ncol(train_data_temp), "cells, ", nrow(train_data_temp), "genes.\n")

# Sequentially merge remaining datasets
for (i in 3:length(train_data_list)) {
  cat(paste0(i, ". Merging with train", i, "...\n"))
  
  # Convert gene names to uppercase
  rownames(train_data_list[[i]]) <- toupper(rownames(train_data_list[[i]]))
  rownames(train_data_temp) <- toupper(rownames(train_data_temp))
  
  # Keep only common genes
  common_genes <- intersect(rownames(train_data_temp), rownames(train_data_list[[i]]))
  train_data_temp <- train_data_temp[common_genes, , drop = FALSE]
  train_data_list[[i]] <- train_data_list[[i]][common_genes, , drop = FALSE]
  
  # Merge
  train_data_temp <- cbind(train_data_temp, train_data_list[[i]])
  pheno_data_temp <- rbind(pheno_data_temp, pheno_data_list[[i]])
  
  cat("   After merge: ", ncol(train_data_temp), "cells, ", nrow(train_data_temp), "genes.\n")
}

# Assign final training data
train_data <- train_data_temp
phenoData_train <- pheno_data_temp
rm(train_data_temp, pheno_data_temp, train_data_list, pheno_data_list)
gc()

cat("\nFinal training data dimensions:", dim(train_data), "\n")
cat("Cell type distribution in training data:\n")
print(table(phenoData_train$cellType))

# ============================================================================
# SECTION 4: LOAD AND PREPROCESS TEST DATA
# ============================================================================
# Load independent test dataset for validation
# ============================================================================

cat("\n--- Loading Test Dataset ---\n")

# Define test dataset path
test_dir <- "/media/desk16/tjn050/LSC/HUSCH/Adipose/HU_0303_Adipose_GSE153643/"

# Load test data using the same helper function
test_dataset <- load_husch_dataset(test_dir, "test")
test_data <- test_dataset$expr
phenoData_test <- test_dataset$pheno
rm(test_dataset)

# ----------------------------------------------------------------------------
# Step 1: Select Cell Types of Interest
# ----------------------------------------------------------------------------
# For adipose tissue analysis, we focus on three major cell types:
# CD4T (CD4+ T cells), Fibroblast, and Mono_Macro (Monocytes/Macrophages)
select_cellType <- c("CD4T", "Fibroblast", "Mono_Macro")

cat("\n--- Selecting Cell Types of Interest ---\n")
cat("Selected cell types:", paste(select_cellType, collapse = ", "), "\n")

# Filter training data
cat("\nFiltering training data...\n")
original_train_cells <- ncol(train_data)
phenoData_train <- phenoData_train[phenoData_train$cellType %in% select_cellType, ]
train_data <- train_data[, phenoData_train$cellID, drop = FALSE]
cat("  Retained", ncol(train_data), "of", original_train_cells, "cells (", 
    round(ncol(train_data)/original_train_cells*100, 1), "%).\n")
cat("  Training cell type distribution:\n")
print(table(phenoData_train$cellType))

# Filter test data
cat("\nFiltering test data...\n")
original_test_cells <- ncol(test_data)
phenoData_test <- phenoData_test[phenoData_test$cellType %in% select_cellType, ]
test_data <- test_data[, phenoData_test$cellID, drop = FALSE]
cat("  Retained", ncol(test_data), "of", original_test_cells, "cells (", 
    round(ncol(test_data)/original_test_cells*100, 1), "%).\n")
cat("  Test cell type distribution:\n")
print(table(phenoData_test$cellType))

# ----------------------------------------------------------------------------
# Step 2: Integrate Training and Test Data
# ----------------------------------------------------------------------------
cat("\n--- Final Data Integration ---\n")

# Ensure gene name consistency
rownames(train_data) <- toupper(rownames(train_data))
rownames(test_data) <- toupper(rownames(test_data))

# Keep only genes present in both training and test data
common_genes <- intersect(rownames(train_data), rownames(test_data))
cat("Number of common genes between training and test:", length(common_genes), "\n")

train_data <- train_data[common_genes, , drop = FALSE]
test_data <- test_data[common_genes, , drop = FALSE]

# Combine all data for potential downstream use
data <- cbind(train_data, test_data)
full_phenoData <- rbind(phenoData_train, phenoData_test)

# ----------------------------------------------------------------------------
# Step 3: Final Data Summary
# ----------------------------------------------------------------------------
cat("\n=== DATA LOADING SUMMARY ===\n")
cat("Training data (reference):", ncol(train_data), "cells,", nrow(train_data), "genes\n")
cat("Test data (validation):", ncol(test_data), "cells,", nrow(test_data), "genes\n")
cat("Combined data:", ncol(data), "cells,", nrow(data), "genes\n\n")

cat("Final cell type distribution:\n")
print(table(full_phenoData$cellType))

# Remove temporary variables to clean workspace
rm(common_genes, original_train_cells, original_test_cells, 
   train_dirs, train_prefixes, select_cellType, test_dir, 
   load_husch_dataset, i)
gc()

cat("\nSections 3 and 4 completed successfully.\n")
# ============================================================================


# ============================================================================
# SECTION 5: QUALITY CONTROL
# ============================================================================
# Continue with the original workflow...
# Note: The following sections remain unchanged from the original script



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


#sim4 Adipose (人类同组织不同测序方法)
prio_markers=list(
  CD4T=c("CCR7","CD28","CD4","IL7R","TCF7"),
  #CD8T=c("CD8A","CD8B"),
  #NK=c("GNLY","IFNG","KLRB1","KLRC1","KLRC2","KLRD1","NKG7")
  Fibroblast =c("ACTA2","COL1A1","COL1A2","COL3A1","COL6A1","COL6A2","COL6A3","DCN","FBLN1","FBLN2","FBN1","FN1","GSN","LY6A","MGP","MYL9","PDGFRA","PI16","POSTN"),
  Mono_Macro=c("CD14","CD68","CSF1R","CST3","FABP4","MMP9")
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
methods_to_use <- c("nnls", "OLS", "ridge")  # 消融/调参脚本：仅跑 3 个确定性方法用于复现验证

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
if (FALSE) {  # 调参/消融脚本：跳过 single-cell 方法（MuSiC/CIBERSORTx），不参与调参/消融
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
results <- read.table("/media/desk16/tjn050/LSC/LSC_deconvolution/LSC/results/CIBERSORTx_Job119_sim42_Results.csv", 
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
}  # end if(FALSE) — 跳过 Section 10 single-cell 方法
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

# ============================================================================
# SECTION 12.5: INPUT REPRODUCIBILITY CHECK（输入复现校验）
# ESMF 重跑（原 3000 迭代，~100 分钟）已跳过以节省时间。
# 本段双重验证输入与正文基准一致，全部通过才继续消融：
#   1) 3 个确定性方法（nnls/OLS/ridge）结果 vs 基准 deconvolution_summary_20260512_000914.txt
#   2) ESMF 输入（V/sig_matrix/sd_matrix/coef1/marker_distrib）vs 基准参数文件 *_20260816_002201.txt
# 任一不一致立即终止，避免基于漂移输入做消融。
# ============================================================================
bench_ok <- TRUE
baseline_methods <- list(
  nnls  = c(RMSE = 0.0953, Pearson = 0.9444),
  OLS   = c(RMSE = 0.1749, Pearson = 0.8503),
  ridge = c(RMSE = 0.1258, Pearson = 0.8795))
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
baseline_ts <- "20260512_002853"   # 正文基准 ESMF 参数 timestamp（ESMF 基准 0.0731/0.9629 那次；注意 20260513_111602 是 252 行基因、非正文最终版）
baseline_param_dir <- file.path(
  "/media/desk16/tjn050/LSC/LSC_deconvolution/LSC/results/MultiSourceReferenceIntegration_Hm_Adipose",
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
  cat(sprintf("[DEBUG] %-12s ref dim=%s cur dim=%s\n", obj_name,
              paste(dim(ref), collapse="x"), paste(dim(cur), collapse="x")))
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
  cat("\n>>> INPUT REPRODUCED: benchmark + ESMF 输入均与基准一致，开始消融实验。\n")
} else {
  cat("\n>>> FATAL: 输入与基准不一致！数据/代码可能已变动，消融实验已终止。\n")
  stop("Input reproducibility check failed. Aborting ablation study.")
}





# ============================================================================
# SECTION 15: ABLATION STUDY  --  A1 (正则化消融) + A2 (marker 面板消融) + E3 (多次初始化稳定性)
# 仅追加到此副本脚本 ablation-study-multi-source-reference-integration-hm-adipose.R。
# 正文基准（multi-source-reference-integration-hm-adipose.R）保持不动。
# 单变量控制铁律：每一组只改一个变量，其余全部照抄基准。
# 本段复用上方已计算的工作空间对象：V, train, group, generator, P,
#   markers, prio_markers, md, sig_matrix, sd_matrix, ML, r, n, m, para1/2/3, coef1。
# ============================================================================

# ---- 局部辅助函数：按给定 md 拟合 sig_matrix/sd_matrix/ML/coef1 ----
fit_signature <- function(md_obj) {
  require(fitdistrplus)
  # 修复: 按 gene 去重。prio_markers 中 SCGB2A1 同时属于 Alpha 和 Beta，导致 md_obj$gene
  #       含同一基因两次，train[md_obj$gene, ] 会返回重复行（重复行名），
  #       使 sig_tr/sig_new 产生重复 rownames，后续 rownames(W_ord) <- rownames(sigm) 维度不匹配。
  #       同一基因保留第一次出现（首个 CT）。
  md_obj <- md_obj[!duplicated(md_obj$gene), , drop = FALSE]
  sig_tr <- train[md_obj$gene, , drop = FALSE]
  # 二次保险: 若 train 行本身含重复（极端情况），按行名去重
  if (anyDuplicated(rownames(sig_tr))) {
    sig_tr <- sig_tr[!duplicated(rownames(sig_tr)), , drop = FALSE]
    md_obj <- md_obj[!duplicated(md_obj$gene), , drop = FALSE]
    md_obj <- md_obj[match(rownames(sig_tr), md_obj$gene), , drop = FALSE]
  }
  ct_names <- names(table(md_obj$CT))
  mu_ls <- list(); sd_ls <- list()
  for (i in seq_along(ct_names)) {
    means <- c(); sds <- c()
    cell_index <- which(colnames(sig_tr) %in% ct_names[i])
    counts <- sig_tr[, cell_index, drop = FALSE]
    for (j in seq_len(nrow(counts))) {
      if (length(which(round(counts[j, ]) != 0)) > 0) {
        f1 <- fitdist(round(counts[j, ]), "nbinom", method = "mse")
        mu_ <- ifelse(is.na(f1$estimate[["mu"]]), 0, f1$estimate[["mu"]])
        si_ <- ifelse(is.na(f1$estimate[["size"]]), 0, f1$estimate[["size"]])
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
  # coef1 重建
  r_new <- length(ct_names); m_new <- nrow(sig_tr)
  coef1_new <- matrix(para1, m_new, r_new)
  for (i in seq_len(r_new)) { gix <- which(md_obj$CT %in% ct_names[i]); coef1_new[gix, i] <- para2 }
  for (i in seq_len(r_new)) { gix <- which(md_obj$gene %in% prio_markers[[i]]); if (length(gix) > 0) coef1_new[gix, i] <- para3 }
  list(sig = sig_new, sd = sd_new, ML = ML_new, coef1 = coef1_new, md = md_obj, r = r_new, m = m_new)
}

# ---- 局部辅助函数：向量化 ESMF 内核（与 esmf_deconvolution.R 数学等价，仅提速）----
# 说明：W/H 更新公式与原版逐元素一致（相同 seed 下 50 次迭代已验证逐位 identical）。
#       err 只用于提前终止判断；err_cutoff=0 且 err 单调下降时不触发终止，
#       因此 3000 次迭代的 W/H 与原版一致。去掉每迭代 3 行 message 日志，实测提速约 20 倍。
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

# ---- 局部辅助函数：ESMF 运行 + H/W 评估（照搬 SECTION 13 逻辑）----
# v3 修改说明（恢复 W 指标，口径对齐正文）：
#   正文 generalization-from-single-reference-hm-pancreas.R 第 1264-1290 行报告 W_Pearson=0.8641，
#   是 ESMF 优于对照方法（DSA 0.6928）的卖点，消融不应缺失。
#   v1 的 W 指标为负（-0.0677）是因为未做列对齐（W 列与 sigm 列错位），非"尺度不可比"。
#   v3 用已列对齐的 W_ord（按 cell_order_loc 对齐、命名 sigm 行/列名）与 sigm 逐元素
#   做 RMSE/Pearson（melt + merge），严格复现正文口径。W_ord 未做列归一化（与正文一致）。
#   evaluate_esmf 返回：H 指标 + W_RMSE/W_Pearson + 列对齐/列归一化后的 W 与 H（供 E3 复用）。
evaluate_esmf <- function(rate_mat, sig_mat, sd_mat, md_obj, label,
                          seed = 1233322, return_matrices = FALSE) {
  require(clue)
  pMatrix.min <- function(A, B) {
    n <- nrow(A); D <- matrix(NA, n, n)
    for (i in 1:n) for (j in 1:n) D[j, i] <- (sum((B[j, ] - A[i, ])^2))
    vec <- c(solve_LSAP(D)); list(A = A[vec, ], pvec = vec)
  }
  fit <- fit_signature(md_obj)
  # v3 修复：A1（md_obj 基因集合与全局 sig_mat 一致，即 benchmark/rate 档）应直接用传入的
  #   sig_mat/sd_mat（= 全局 sig_matrix/sd_matrix，与正文完全一致），避免 fit_signature 重新
  #   拟合引入数值差异，导致 W 与 sig 逐元素相关崩塌为负（W_Pearson=-0.0677）。
  #   仅当 md_obj 换 marker（A2，基因集合与全局 sig 不同）时才用 fit_signature 重新拟合签名。
  if (setequal(md_obj$gene, rownames(sig_mat))) {
    sigm <- sig_mat; sdm <- sd_mat
  } else {
    sigm <- fit$sig; sdm <- fit$sd
  }
  MLm <- fit$ML; coef1m <- fit$coef1
  rm_ <- fit$r; mm_ <- fit$m
  # 对齐 V 与 sig
  g_f <- intersect(rownames(V), rownames(sigm))
  V_loc <- V[g_f, , drop = FALSE]
  sigm <- sigm[g_f, , drop = FALSE]
  sdm <- sdm[g_f, , drop = FALSE]
  # rate 需与 sigm 同维度
  # 修复: coef1/coef1m 是裸矩阵（matrix(para1, m, r)）无行名，而 g_f 是基因名字符向量，
  #       直接字符索引会报 "no 'dimnames' attribute for array"。
  #       行序已与 sigm 一致（coef1 行序 = md 行序 = sig_matrix 行序；coef1m 行序 = md_obj 行序）。
  #       注意：必须用 fit_signature 输出的完整基因名（nrow(fit$sig) 行）补行名——
  #       当 md 基因不全在 V 里时 g_f 长度 < 完整基因数，若按裁剪后的 rownames(sigm) 补，
  #       会报 "length of 'dimnames' [1] not equal to array extent"。
  if (is.matrix(rate_mat) && is.null(rownames(rate_mat))) {
    if (nrow(rate_mat) == nrow(fit$sig)) rownames(rate_mat) <- rownames(fit$sig)
    else rownames(rate_mat) <- rownames(sigm)
  }
  if (is.null(rownames(coef1m))) {
    rownames(coef1m) <- rownames(fit$sig)
  }
  if (is.character(rate_mat) && length(rate_mat) == 1 && rate_mat == "coef") {
    # A2: 使用按新 md 重建的 coef1（fit_signature 内部生成，保证单变量控制）
    rate_loc <- coef1m[g_f, , drop = FALSE]
  } else if (identical(dim(rate_mat), dim(sigm))) {
    rate_loc <- rate_mat[g_f, , drop = FALSE]
  } else if (is.matrix(rate_mat) && nrow(rate_mat) == mm_ && ncol(rate_mat) == rm_) {
    rate_loc <- rate_mat[g_f, , drop = FALSE]
  } else {
    rate_loc <- matrix(rate_mat, nrow(sigm), ncol(sigm))
  }
  P_loc <- generator[["P"]]
  set.seed(seed)
  RES_loc <- ESMF_deconvolution_fast(V = V_loc, ml = MLm, r = rm_, err_cutoff = 0,
                                rate = rate_loc, sig_matrix = sigm,
                                sd_matrix = sdm, iter_times = 3000)
  W_loc <- RES_loc$W; H_loc <- RES_loc$H
  n_loc <- ncol(V_loc)
  # ---- H 比例评估 ----
  corW_loc <- cor(W_loc, sigm)
  rownames(corW_loc) <- as.character(1:nrow(corW_loc))
  B_d <- diag(1, nrow(corW_loc))
  X_loc <- pMatrix.min(corW_loc, B_d)
  cell_order_loc <- X_loc$pvec
  # 归一化 H
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
  # ---- W 列对齐（供 W vs sigm 评估 + E3 跨 seed 稳定性复用）----
  W_ord <- W_loc[, cell_order_loc, drop = FALSE]
  if (nrow(W_ord) != length(rownames(sigm))) {
    stop(sprintf("[evaluate_esmf] dimension mismatch in '%s': nrow(W_ord)=%d, length(rownames(sigm))=%d, anyDuplicated(rownames(sigm))=%d, g_f_len=%d",
                 label, nrow(W_ord), length(rownames(sigm)),
                 anyDuplicated(rownames(sigm)), length(g_f)))
  }
  colnames(W_ord) <- colnames(sigm); rownames(W_ord) <- rownames(sigm)
  # ---- v3：恢复 W vs sigm 逐元素评估（口径对齐正文第 1264-1290 行）----
  # 用已列对齐、已命名的 W_ord（未归一化）与 sigm 对比，melt + merge 后算 RMSE/Pearson
  Wm <- t(W_ord)
  Wm <- as.data.table(Wm, keep.rownames = TRUE)
  Wm <- data.table::melt(Wm, id.vars = "rn")
  colnames(Wm) <- c("CT", "gene", "observed_values")
  Sm <- t(sigm)
  Sm <- as.data.table(Sm, keep.rownames = TRUE)
  Sm <- data.table::melt(Sm, id.vars = "rn")
  colnames(Sm) <- c("CT", "gene", "expected_values")
  INDw <- merge(Wm, Sm, by = c("CT", "gene"))
  INDw$observed_values <- round(INDw$observed_values, 3)
  INDw$expected_values <- round(INDw$expected_values, 3)
  INDw_sum <- INDw %>% dplyr::summarise(
    RMSE = sqrt(mean((observed_values - expected_values)^2)) %>% round(., 4),
    Pearson = cor(observed_values, expected_values) %>% round(., 4))
  # 列归一化（每列和=1）：使 W 落到与 H 同构的比例尺度，消除尺度不确定性
  csW <- colSums(W_ord)
  csW[csW == 0] <- 1
  W_norm <- sweep(W_ord, 2, csW, FUN = "/")
  out <- c(label = label,
           H_RMSE = as.numeric(INDh$RMSE), H_Pearson = as.numeric(INDh$Pearson),
           W_RMSE = as.numeric(INDw_sum$RMSE), W_Pearson = as.numeric(INDw_sum$Pearson))
  if (return_matrices) {
    attr(out, "W_norm") <- W_norm
    attr(out, "H_norm") <- H_loc[gtools::mixedsort(rownames(H_loc)), , drop = FALSE]
    attr(out, "cell_order") <- cell_order_loc
    attr(out, "sigm") <- sigm
  }
  out
}

# ============================================================================
# A1: 正则化消融（锁死 marker = 正文基准 md；只改 rate）
# ============================================================================
cat("\n\n========== A1: REGULARIZATION ABLATION ==========\n")
a1_factor <- c("benchmark_rate1" = 1, "pure_NMF_rate0" = 0, "weak_rate0.1" = 0.1, "strong_rate10" = 10)
a1_out <- lapply(names(a1_factor), function(nm) {
  fac <- a1_factor[[nm]]
  cat("Running A1:", nm, " rate factor =", fac, "\n")
  rate_mat <- coef1 * fac
  evaluate_esmf(rate_mat = rate_mat, sig_mat = sig_matrix, sd_mat = sd_matrix,
                md_obj = md, label = paste0("A1_", nm))
})
a1_df <- as.data.frame(do.call(rbind, a1_out), stringsAsFactors = FALSE)

# ============================================================================
# A2: marker 面板消融（锁死模型 rate=coef1；只换 md 并据此重建 sig/sd/ML）
# 档位: 基准 / 文献增强(基准+prio) / DE top genes / 随机基因(固定 seed)
# ============================================================================
cat("\n\n========== A2: MARKER PANEL ABLATION ==========\n")
a2_out <- list()

# A2-1 基准 marker 面板（= 正文，rate 保持 coef1 不缩放）
cat("Running A2: benchmark marker panel (正文)\n")
a2_out[["benchmark_marker"]] <- evaluate_esmf(
  rate_mat = "coef", sig_mat = sig_matrix, sd_mat = sd_matrix,
  md_obj = md, label = "A2_benchmark_marker")

# A2-2 文献 marker 增强：在基准 DE 面板 md 之上，额外并入 prio_markers（文献先验）
#       关键修正：不再把文献 marker 单独当签名矩阵用（基因数过少必然崩），
#       而是"基准 + 文献先验增强"，单变量只多加了文献基因，直接测文献先验是否有增量。
cat("Running A2: literature marker augmented (基准 md + prio_markers)\n")
# 注意：md 有 4 列（gene/log2FC/CT/AveExpr），rbind 新增行必须补齐同列，否则报"列数不匹配"。
#       fit_signature 只用 gene/CT 两列，log2FC/AveExpr 填 NA 不影响拟合。
lit_aug_md <- md
for (ct in names(prio_markers)) {
  extra <- setdiff(intersect(prio_markers[[ct]], rownames(train)),
                   lit_aug_md$gene[lit_aug_md$CT == ct])
  if (length(extra) > 0) {
    add <- data.frame(gene = extra, log2FC = NA_real_, CT = ct, AveExpr = NA_real_,
                      stringsAsFactors = FALSE)
    lit_aug_md <- rbind(lit_aug_md, add)
  }
}
# 显式去重：一个基因只保留首个 CT（SCGB2A1 在 Alpha/Beta 重复，fit_signature 内部也会去重，这里先做保证干净）
lit_aug_md <- lit_aug_md[!duplicated(lit_aug_md$gene), , drop = FALSE]
cat(sprintf("[A2] 文献增强面板大小: %d 基因 (基准 md=%d)\n", nrow(lit_aug_md), nrow(md)))
a2_out[["literature_augmented"]] <- evaluate_esmf(
  rate_mat = "coef", sig_mat = sig_matrix, sd_mat = sd_matrix,
  md_obj = lit_aug_md, label = "A2_literature_augmented")

# A2-3 DE top genes：取现有 markers 每型 top 50（统计差异前 N，削弱生物先验）
cat("Running A2: DE top genes (每型 top 50)\n")
de_md <- do.call(rbind, lapply(unique(markers$CT), function(ct) {
  idx <- which(markers$CT == ct)
  if (length(idx) > 50) idx <- idx[1:50]
  markers[idx, c("gene", "CT")]
}))
a2_out[["DE_top50"]] <- evaluate_esmf(
  rate_mat = "coef", sig_mat = sig_matrix, sd_mat = sd_matrix,
  md_obj = de_md, label = "A2_DE_top50")

# A2-4 随机基因（每型随机抽取与基准 md 等量基因，固定 seed 复现）——阴性对照
cat("Running A2: random genes (negative control)\n")
set.seed(999)
rand_md <- do.call(rbind, lapply(unique(md$CT), function(ct) {
  n_ct <- sum(md$CT == ct)
  pool <- setdiff(rownames(train), md$gene[md$CT == ct])
  if (length(pool) < n_ct) pool <- rownames(train)
  data.frame(gene = sample(pool, min(n_ct, length(pool))), CT = ct, stringsAsFactors = FALSE)
}))
a2_out[["random_genes"]] <- evaluate_esmf(
  rate_mat = "coef", sig_mat = sig_matrix, sd_mat = sd_matrix,
  md_obj = rand_md, label = "A2_random_genes")

a2_df <- as.data.frame(do.call(rbind, a2_out), stringsAsFactors = FALSE)

# ============================================================================
# E3: 多次初始化（multi-seed）稳定性 —— 回应 R2 M2①"报告 W/H 变异性"
# 路线 A：W 不作为"W vs sigm 绝对误差"对外报告（尺度不可比，无信息量），
#        改报 W 自身的跨 seed 稳定性：不同随机初始化下，列归一化后的 W 两两相关。
#        同时报告 H 的跨 seed 稳定性。锁死 marker=基准 md、rate=coef1，只变 seed。
# ============================================================================
cat("\n\n========== E3: MULTI-SEED STABILITY (W/H 变异性) ==========\n")
e3_seeds <- 1233322 + 0:9   # 10 个不同随机初始化 seed
e3_n <- length(e3_seeds)
# 逐个 seed 跑基准配置，收集列归一化后的 W 与 H
e3_W <- vector("list", e3_n)
e3_H <- vector("list", e3_n)
for (s in seq_along(e3_seeds)) {
  cat(sprintf("Running E3 seed %d/%d (seed=%d)\n", s, e3_n, e3_seeds[s]))
  res <- evaluate_esmf(rate_mat = "coef", sig_mat = sig_matrix, sd_mat = sd_matrix,
                       md_obj = md, label = paste0("E3_seed_", e3_seeds[s]),
                       seed = e3_seeds[s], return_matrices = TRUE)
  e3_W[[s]] <- attr(res, "W_norm")
  e3_H[[s]] <- attr(res, "H_norm")
}

# 计算跨 seed 两两相关（Pearson）。W 按 (基因 × 细胞型) 展平成向量后比较；
# 由于 W 列已按 cell_order 对齐 + 列归一化，跨 seed 的列对应关系一致，可直接展平比较。
e3_pairwise <- function(mat_list) {
  n <- length(mat_list)
  M <- matrix(NA, n, n, dimnames = list(seq_len(n), seq_len(n)))
  for (i in seq_len(n)) for (j in seq_len(n)) {
    M[i, j] <- cor(as.numeric(mat_list[[i]]), as.numeric(mat_list[[j]]),
                   use = "complete.obs")
  }
  M
}
W_cor <- e3_pairwise(e3_W)
H_cor <- e3_pairwise(e3_H)

# 汇总统计：下三角（i<j）的两两相关均值/最小值（稳健地反映一致性）
lower_tri <- function(M) M[lower.tri(M)]
W_mean <- round(mean(lower_tri(W_cor)), 4)
W_min  <- round(min(lower_tri(W_cor)), 4)
H_mean <- round(mean(lower_tri(H_cor)), 4)
H_min  <- round(min(lower_tri(H_cor)), 4)

e3_summary <- data.frame(
  label = c("E3_W_stability", "E3_H_stability"),
  H_RMSE = NA, H_Pearson = NA,
  metric = c("W_Pearson_mean", "H_Pearson_mean"),
  value = c(W_mean, H_mean),
  min_pairwise = c(W_min, H_min),
  n_seeds = e3_n,
  stringsAsFactors = FALSE)

cat("\n--- E3 W 跨 seed 稳定性（列归一化后两两 Pearson）---\n")
print(round(W_cor, 4))
cat(sprintf("W 两两相关: mean=%.4f, min=%.4f\n", W_mean, W_min))
cat("\n--- E3 H 跨 seed 稳定性（两两 Pearson）---\n")
print(round(H_cor, 4))
cat(sprintf("H 两两相关: mean=%.4f, min=%.4f\n", H_mean, H_min))

# ============================================================================
# 汇总并保存 A1/A2/E3 消融结果（写入独立目录，不覆盖基准）
# ============================================================================
ablation_out <- rbind(a1_df, a2_df)
ablation_file <- paste0(results_dir, "ablation_results_", current_datetime, ".txt")
e3_file <- paste0(results_dir, "e3_multiseed_stability_", current_datetime, ".txt")
cat("\n\n========== ABLATION STUDY SUMMARY (A1 + A2) ==========\n")
print(ablation_out)
write.table(ablation_out, file = ablation_file, sep = "\t", quote = FALSE,
            col.names = TRUE, row.names = FALSE)
cat("Ablation results saved to:", ablation_file, "\n")

cat("\n========== E3 MULTI-SEED STABILITY SUMMARY ==========\n")
print(e3_summary)
write.table(e3_summary, file = e3_file, sep = "\t", quote = FALSE,
            col.names = TRUE, row.names = FALSE)
cat("E3 stability results saved to:", e3_file, "\n")
cat("Ablation study (A1 + A2 + E3) completed.\n")

