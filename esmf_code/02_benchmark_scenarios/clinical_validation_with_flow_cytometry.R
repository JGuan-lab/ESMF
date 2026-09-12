# ============================================================================
# Script: ClinicalValidationWithFlowCytometry_RealBulk_Deconvolution.R
# Description: Deconvolution analysis of real whole blood RNA-seq data
#              using multiple methods and validation with flow cytometry ground truth.
# Author: Your Name
# Date: 2026-04-26
# ============================================================================

# ============================================================================
# SECTION 1: ENVIRONMENT SETUP AND PACKAGE LOADING
# ============================================================================

# Set working directory
setwd("/media/desk16/tjn050/LSC/LSC_deconvolution/LSC/")

# Load required CRAN packages
library(NMF)
library(CellMix)
library(MASS)
library(fitdistrplus)
library(scater)
library(Matrix)
library(Seurat)
library(descend)
library(data.table)
library(cowplot)
library(rhdf5)
library(clusterProfiler)
library(limma)
library(pheatmap)
library(gtools)
library(xbioc)
library(matrixStats)
library(dplyr)
library(SingleCellExperiment)
library(MuSiC)

# Load custom functions for deconvolution analysis
source("esmf_deconvolution.R")
source("esmf_helper_functions.R")
source("bulk_deconvolution.R")
source('normalize_col_in_matrix.R')
source("basic_functions.R")
source("cibersort.R")

# ============================================================================
# SECTION 2: ANALYSIS PARAMETER CONFIGURATION
# ============================================================================

# Define core analysis parameters
marker_strategy <- "all"           # Marker gene selection strategy
number_cells <- 100                # Number of cells per pseudo-bulk mixture
to_remove <- "none"                # Cell types to exclude
dataset <- "example"               # Dataset identifier
transformation <- "none"           # Data transformation method
deconv_type <- "bulk"              # Deconvolution framework
method <- "OLS"                    # Deconvolution method for initial run
num_cores <- 1                     # Number of cores for parallel processing

# Normalization parameters
if (deconv_type == "bulk") {
  normalization <- "none"          # Changed to "none" instead of "TPM"
  marker_strategy <- "all"
} else if (deconv_type == "sc") {
  normalization_scC <- "TMM"
  normalization_scT <- "TMM"
} else {
  stop("Please enter a valid deconvolution framework")
}

# Ensure normalization_scC and normalization_scT are always defined
# These variables are needed regardless of deconv_type
normalization_scC <- "TMM"
normalization_scT <- "TMM"

# If normalization is not defined, define it
if (!exists("normalization")) {
  normalization <- "none"
}

# Define cell types of interest
select_cellType <- c("CD4T", "CD8T", "NK")

# Priority marker genes for specific cell types
prio_markers <- list(
  CD4T = c("CCR7", "CD28", "CD4", "IL7R", "TCF7"),
  CD8T = c("CD8A", "CD8B"),
  NK = c("IFNG", "KLRB1", "KLRC1", "KLRC2", "KLRD1")
)

# Define analysis name prefix
analysis_prefix <- "ClinicalValidationWithFlowCytometry_threshold200"

# ============================================================================
# SECTION 3: CREATE RESULTS DIRECTORY
# ============================================================================

# Create results directory based on analysis_prefix
results_dir <- file.path("/media/desk16/tjn050/LSC/LSC_deconvolution/LSC/results", analysis_prefix, "")

# Create results directory (if it doesn't exist)
if (!dir.exists(results_dir)) {
  dir.create(results_dir, recursive = TRUE)
  cat("Results directory created:", results_dir, "\n")
}

# ============================================================================
# SECTION 4: SINGLE-CELL REFERENCE DATA LOADING AND PROCESSING
# ============================================================================

cat("Loading single-cell reference datasets...\n")

# List of training dataset directories
train_dirs <- c(
  "/media/desk16/tjn050/LSC/HUSCH/Blood/HU_0036_Blood_10x/",
  "/media/desk16/tjn050/LSC/HUSCH/Blood/HU_0037_Blood_10x/",
  "/media/desk16/tjn050/LSC/HUSCH/Blood/HU_0038_Blood_10x/",
  "/media/desk16/tjn050/LSC/HUSCH/Blood/HU_0039_Blood_10x/",
  "/media/desk16/tjn050/LSC/HUSCH/Blood/HU_0040_Blood_10x/",
  "/media/desk16/tjn050/LSC/HUSCH/Blood/HU_0041_Blood_10x/"
)

# Initialize objects for combined data
sc_data_combined <- NULL
phenoData_sc_combined <- NULL

# Function to load a single HUSCH dataset
load_husch_dataset <- function(data_dir, dataset_id) {
  dataset_name <- strsplit(data_dir, "/")[[1]][8]
  
  # Load metadata
  phenoData <- read.csv(
    paste0(data_dir, dataset_name, "_meta.txt"),
    sep = "\t",
    header = TRUE
  )
  
  # Clean cell type names
  phenoData$Celltype <- gsub("/", "_", phenoData$Celltype)
  phenoData$Celltype <- gsub(" ", "_", phenoData$Celltype)
  colnames(phenoData) <- c("cellID", "cellType", "Platform")
  
  # Load expression data from HDF5 format
  data_mat <- h5read(paste0(data_dir, dataset_name, "_gene_count.h5"), "matrix")
  
  # Construct sparse matrix
  barcodes <- data_mat$barcodes
  counts <- data_mat$data
  gene_id <- data_mat$features$id
  indices <- data_mat$indices
  indptr <- data_mat$indptr
  shape <- data_mat$shape
  
  mat <- sparseMatrix(
    i = indices[] + 1,
    p = indptr[],
    x = as.numeric(x = counts[]),
    dims = shape[],
    repr = "T"
  )
  colnames(mat) <- barcodes
  rownames(mat) <- gene_id
  expr_data <- as.matrix(mat)
  
  # Add dataset-specific prefix
  colnames(expr_data) <- paste0("dataset", dataset_id, colnames(expr_data))
  phenoData$cellID <- paste0("dataset", dataset_id, phenoData$cellID)
  
  # Ensure consistent cell ordering
  phenoData <- phenoData[phenoData$cellID %in% colnames(expr_data), ]
  expr_data <- expr_data[, intersect(phenoData$cellID, colnames(expr_data))]
  
  return(list(expr_data = expr_data, phenoData = phenoData))
}

# Load and combine all datasets
for (i in seq_along(train_dirs)) {
  cat("Loading dataset", i, ":", train_dirs[i], "\n")
  dataset_result <- load_husch_dataset(train_dirs[i], i)
  
  if (is.null(sc_data_combined)) {
    sc_data_combined <- dataset_result$expr_data
    phenoData_sc_combined <- dataset_result$phenoData
  } else {
    # Standardize gene names to uppercase
    rownames(sc_data_combined) <- toupper(rownames(sc_data_combined))
    rownames(dataset_result$expr_data) <- toupper(rownames(dataset_result$expr_data))
    
    # Keep only genes present in both existing and new datasets
    common_genes <- intersect(rownames(sc_data_combined), rownames(dataset_result$expr_data))
    sc_data_combined <- sc_data_combined[common_genes, ]
    dataset_result$expr_data <- dataset_result$expr_data[common_genes, ]
    
    # Merge with existing data
    sc_data_combined <- cbind(sc_data_combined, dataset_result$expr_data)
    phenoData_sc_combined <- rbind(phenoData_sc_combined, dataset_result$phenoData)
  }
  
  cat("  Cumulative dimensions:", dim(sc_data_combined), "\n")
  cat("  Cell type distribution:\n")
  print(table(phenoData_sc_combined$cellType))
}

cat("Final single-cell reference data dimensions:", dim(sc_data_combined), "\n")

# Select cell types of interest
cat("\nSelecting cell types of interest:", paste(select_cellType, collapse = ", "), "\n")
phenoData_sc_combined <- phenoData_sc_combined[phenoData_sc_combined$cellType %in% select_cellType, ]
phenoData_sc_combined <- phenoData_sc_combined[phenoData_sc_combined$cellID %in% colnames(sc_data_combined), ]
sc_data_combined <- sc_data_combined[, intersect(phenoData_sc_combined$cellID, colnames(sc_data_combined))]
cat("Single-cell data after cell type selection:", dim(sc_data_combined), "\n")

# ============================================================================
# SECTION 5: BULK EXPRESSION DATA LOADING
# ============================================================================

cat("\nLoading bulk whole blood RNA-seq data...\n")

# Load bulk expression data
bulk_expr_file <- "/media/desk16/tjn050/LSC/realbulk/Fig2b-WholeBlood_RNAseq.txt"
bulk_expr_data <- read.table(bulk_expr_file, sep = "\t", header = TRUE)
rownames(bulk_expr_data) <- bulk_expr_data[, 1]
bulk_expr_data <- bulk_expr_data[, -1]

# Load ground truth proportions from flow cytometry
ground_truth_file <- "/media/desk16/tjn050/LSC/realbulk/Fig2b_ground_truth_whole_blood.txt"
ground_truth <- read.table(ground_truth_file, sep = "\t", header = TRUE)
ground_truth <- t(ground_truth)
colnames(ground_truth) <- ground_truth[1, ]
ground_truth <- ground_truth[-1, ]
ground_truth <- ground_truth[-4, ]  # Remove unnecessary row
ground_truth <- apply(ground_truth, 2, as.numeric)
ground_truth <- as.data.frame(ground_truth)
rownames(ground_truth) <- c("Neutrophil", "Lymphatic_Endothelial", "Mono_Macro", 
                            "CD8T", "CD4T", "B", "NK")

# Select relevant cell types for comparison
ground_truth_selected <- ground_truth[c("CD4T", "CD8T", "NK"), , drop = FALSE]
rownames(ground_truth_selected) <- c("CD4T", "CD8T", "NK")

cat("Bulk expression data dimensions:", dim(bulk_expr_data), "\n")
cat("Ground truth proportions dimensions:", dim(ground_truth_selected), "\n")

# ============================================================================
# SECTION 6: DATA PREPROCESSING AND QUALITY CONTROL
# ============================================================================

cat("\nAligning single-cell and bulk datasets...\n")

# Standardize gene names
rownames(sc_data_combined) <- toupper(rownames(sc_data_combined))
rownames(bulk_expr_data) <- toupper(rownames(bulk_expr_data))

# Keep only genes present in both datasets
common_genes <- intersect(rownames(sc_data_combined), rownames(bulk_expr_data))
sc_data_aligned <- sc_data_combined[common_genes, ]
bulk_expr_aligned <- bulk_expr_data[common_genes, ]

cat("Aligned single-cell data dimensions:", dim(sc_data_aligned), "\n")
cat("Aligned bulk data dimensions:", dim(bulk_expr_aligned), "\n")

# Quality control for single-cell data
cat("\nPerforming quality control on single-cell data...\n")

# Function to filter outlier cells
filter_outlier_cells <- function(filter_param) {
  cells_to_remove <- which(
    filter_param > median(filter_param) + 3 * mad(filter_param) |
      filter_param < median(filter_param) - 3 * mad(filter_param)
  )
  return(cells_to_remove)
}

# Calculate library sizes and gene content metrics
lib_sizes <- colSums(sc_data_aligned)
gene_names <- rownames(sc_data_aligned)

# Identify mitochondrial and ribosomal genes
mt_genes <- grepl("^MT-|_MT-", gene_names, ignore.case = TRUE)
rb_genes <- grepl("^RPL|^RPS|_RPL|_RPS", gene_names, ignore.case = TRUE)

mt_percent <- colSums(sc_data_aligned[mt_genes, ]) / lib_sizes
rb_percent <- colSums(sc_data_aligned[rb_genes, ]) / lib_sizes

# Identify outlier cells
cells_to_remove <- lapply(
  list(lib_size = lib_sizes, mt_percent = mt_percent, rb_percent = rb_percent),
  filter_outlier_cells
) %>% unlist() %>% unique()

# Remove outlier cells
if (length(cells_to_remove) != 0) {
  sc_data_aligned <- sc_data_aligned[, -cells_to_remove]
  phenoData_sc_combined <- phenoData_sc_combined[-cells_to_remove, ]
}

# Filter lowly expressed genes
keep_genes <- which(Matrix::rowSums(sc_data_aligned > 0) >= round(0.05 * ncol(sc_data_aligned)))
sc_data_filtered <- sc_data_aligned[keep_genes, ]

cat("Dimensions after quality control:", dim(sc_data_filtered), "\n")

# Prepare training data
cat("\nPreparing training data for marker selection...\n")

# Map cell type annotations
original_cell_names <- colnames(sc_data_filtered)
colnames(sc_data_filtered) <- as.character(
  phenoData_sc_combined$cellType[match(colnames(sc_data_filtered), phenoData_sc_combined$cellID)]
)

# Keep cell types with sufficient cell counts
cell_counts <- table(colnames(sc_data_filtered))
celltypes_to_keep <- names(cell_counts)[cell_counts >= 200]
pData_train <- phenoData_sc_combined[phenoData_sc_combined$cellType %in% celltypes_to_keep, ]
sc_data_filtered <- sc_data_filtered[, colnames(sc_data_filtered) %in% celltypes_to_keep]
original_cell_names <- original_cell_names[colnames(sc_data_filtered) %in% celltypes_to_keep]

# All data is used for training in this real bulk analysis
training_idx <- which(original_cell_names %in% original_cell_names)
train_data <- sc_data_filtered[, training_idx]
original_cell_names_train <- original_cell_names[training_idx]

cat("Training data dimensions:", dim(train_data), "\n")

# Create training data with cell IDs
train_data_cellID <- train_data
colnames(train_data_cellID) <- original_cell_names_train

# ============================================================================
# SECTION 7: REFERENCE MATRIX CONSTRUCTION AND MARKER GENE SELECTION
# ============================================================================

cat("\nConstructing reference expression matrix...\n")

cell_types <- colnames(train_data)
cell_type_groups <- lapply(unique(cell_types), function(x) which(cell_types %in% x))
names(cell_type_groups) <- unique(cell_types)

# Calculate mean expression per cell type
ref_matrix <- lapply(cell_type_groups, function(x) Matrix::rowMeans(train_data[, x]))
ref_matrix <- round(do.call(cbind.data.frame, ref_matrix))

# Calculate expression variance per cell type
ref_var_matrix <- lapply(cell_type_groups, function(x) train_data[, x])
ref_var_matrix <- lapply(ref_var_matrix, function(x) matrixStats::rowSds(Matrix::as.matrix(x)))
ref_var_matrix <- round(do.call(cbind.data.frame, ref_var_matrix))
rownames(ref_var_matrix) <- rownames(train_data)

cat("Reference matrix dimensions:", dim(ref_matrix), "\n")

cat("\nSelecting marker genes...\n")

# Filter genes: expressed in at least 30% of cells per cell type
keep_genes_markers <- sapply(unique(cell_types), function(x) {
  ct_indices <- which(cell_types %in% x)
  min_cells <- ceiling(0.3 * length(ct_indices))
  Matrix::rowSums(train_data[, ct_indices, drop = FALSE] != 0) >= min_cells
})
train_data_filtered <- train_data[Matrix::rowSums(keep_genes_markers) > 0, ]

# Normalize for differential expression analysis
train_data_normalized <- Normalization(train_data_filtered)

# Perform differential expression analysis using limma
annotation <- factor(colnames(train_data_normalized))
design <- model.matrix(~0 + annotation)
colnames(design) <- unlist(lapply(strsplit(colnames(design), "annotation"), function(x) x[2]))

# Create contrast matrix for pairwise comparisons
cont.matrix <- matrix((-1 / ncol(design)), nrow = ncol(design), ncol = ncol(design))
colnames(cont.matrix) <- colnames(design)
diag(cont.matrix) <- (ncol(design) - 1) / ncol(design)

# Differential expression analysis
v <- limma::voom(train_data_normalized, design = design, plot = FALSE)
fit <- limma::lmFit(v, design)
fit_contrasts <- limma::contrasts.fit(fit, cont.matrix)
fit_ebayes <- limma::eBayes(fit_contrasts, trend = TRUE)

# Ensure cont.matrix object exists before calling marker.fc
if (!exists("cont.matrix")) {
  stop("Error: cont.matrix object not found. Please check the contrast matrix creation.")
}

# Identify marker genes
marker_genes <- marker.fc(fit_ebayes, log2.threshold = log2(1))
cat("Initial marker gene counts:\n")
print(table(marker_genes$CT))

# Filter and select top markers
expressed_markers <- train_data[marker_genes$gene, ]
keep_markers <- which(Matrix::rowSums(expressed_markers > 0) >= round(0.5 * ncol(expressed_markers)))
marker_genes <- marker_genes[keep_markers, ]

# Select top 200 markers per cell type
top_markers <- data.frame()
for (ct in unique(marker_genes$CT)) {
  ct_indices <- which(marker_genes$CT == ct)
  if (length(ct_indices) > 200) {
    ct_indices <- ct_indices[1:200]
  }
  top_markers <- rbind(top_markers, marker_genes[ct_indices, ])
}
marker_genes <- top_markers

# Apply marker selection strategy
marker_distribution <- marker_strategies(marker_genes, marker_strategy, ref_matrix)
final_markers <- marker_distribution

cat("Final marker distribution:\n")
print(table(final_markers$CT))

# ============================================================================
# SECTION 8: SIGNATURE MATRIX ESTIMATION
# ============================================================================

cat("\nEstimating signature matrix using negative binomial distribution...\n")

# Extract marker gene expression
marker_expr <- train_data[final_markers$gene, ]

# Initialize parameter lists
mean_list <- list()
dispersion_list <- list()

# Fit negative binomial distribution for each cell type
for (i in 1:length(unique(final_markers$CT))) {
  cell_type_name <- names(table(final_markers$CT))[i]
  means <- vector()
  dispersions <- vector()
  
  cell_indices <- which(colnames(marker_expr) %in% cell_type_name)
  cell_counts <- marker_expr[, cell_indices]
  
  for (j in 1:nrow(cell_counts)) {
    non_zero_cells <- which(round(cell_counts[j, ]) != 0)
    
    if (length(non_zero_cells) > 0) {
      # Fit negative binomial distribution
      nb_fit <- tryCatch({
        fitdist(round(cell_counts[j, ]), "nbinom", method = "mse")
      }, error = function(e) {
        NULL
      })
      
      if (!is.null(nb_fit) && !is.na(nb_fit$estimate[1]) && !is.na(nb_fit$estimate[2])) {
        means <- c(means, nb_fit$estimate["mu"])
        dispersions <- c(dispersions, nb_fit$estimate["size"])
      } else {
        means <- c(means, 0)
        dispersions <- c(dispersions, 0)
      }
    } else {
      means <- c(means, 0)
      dispersions <- c(dispersions, 0)
    }
  }
  
  mean_list <- c(mean_list, list(means))
  dispersion_list <- c(dispersion_list, list(dispersions))
}

# Construct signature matrices
signature_matrix <- matrix(unlist(mean_list), ncol = length(mean_list))
colnames(signature_matrix) <- names(table(final_markers$CT))
rownames(signature_matrix) <- rownames(marker_expr)

dispersion_matrix <- matrix(unlist(dispersion_list), ncol = length(dispersion_list))
sd_matrix <- signature_matrix * (1 + dispersion_matrix)
colnames(sd_matrix) <- names(table(final_markers$CT))
rownames(sd_matrix) <- rownames(marker_expr)

cat("Signature matrix dimensions:", dim(signature_matrix), "\n")

# Create marker list for constrained methods
marker_list <- CellMix::MarkerList()
marker_list@.Data <- tapply(as.character(final_markers$gene), as.character(final_markers$CT), list)


# ============================================================================
# SECTION 9: BULK DECONVOLUTION METHODS EVALUATION
# ============================================================================

cat("\nEvaluating bulk deconvolution methods...\n")

# 9.1 Prepare data for bulk deconvolution
# ----------------------------------------------------------------------------

# Get ground truth proportions
P <- ground_truth_selected
cat("Ground truth proportions dimensions:", dim(P), "\n")

# Align bulk expression data with signature matrix
common_genes <- intersect(rownames(bulk_expr_aligned), rownames(signature_matrix))
T_aligned <- bulk_expr_aligned[common_genes, ]
C_aligned <- signature_matrix[common_genes, ]

cat("Aligned bulk expression data dimensions:", dim(T_aligned), "\n")
cat("Aligned signature matrix dimensions:", dim(C_aligned), "\n")

# Transform and normalize data
T_transformed <- Transformation(T_aligned, transformation)
C_transformed <- Transformation(C_aligned, transformation)

T_scaled <- Scaling(T_transformed, normalization)
C_scaled <- Scaling(C_transformed, normalization)

# Apply marker selection strategy
marker_distrib <- marker_strategies(marker_genes, marker_strategy, C_scaled)

# 9.2 Evaluate multiple deconvolution methods
# ----------------------------------------------------------------------------
methods_to_use <- c("nnls", "FARDEEP", "RLR", "DCQ", "elastic_net", 
                    "lasso", "ridge", "OLS", "EPIC", "DSA")

# Initialize results data frame
samples <- colnames(P)
results_df <- as.data.frame(matrix(NA, nrow = length(samples), ncol = length(methods_to_use)))
rownames(results_df) <- samples
colnames(results_df) <- methods_to_use

cat("\nStarting deconvolution methods evaluation...\n")
cat(strrep("-", 60), "\n", sep = "")

# Get pDataC (phenotype data for reference)
pDataC <- pData_train

# Set additional parameters
elem <- NULL
STRING <- "ClinicalValidationWithFlowCytometry"

# Evaluate each deconvolution method
for (method in methods_to_use) {
  cat("Evaluating method:", method, "\n")
  
  tryCatch({
    set.seed(123456)
    
    # Execute deconvolution
    deconv_result <- deconv_lsc(
      T = T_scaled,
      C = C_scaled,
      method = method,
      phenoDataC = pDataC,
      P = P,
      elem = elem,
      STRING = STRING,
      marker_distrib = marker_distrib,
      refProfiles.var = ref_var_matrix
    )
    
    if (is.null(deconv_result) || nrow(deconv_result) == 0) {
      cat("  Warning: Deconvolution returned empty result\n")
      next
    }
    
    # Calculate sample-wise Pearson correlations
    estimated_matrix <- deconv_result
    true_matrix <- P
    
    # Ensure matrices are aligned
    common_samples <- intersect(colnames(estimated_matrix), colnames(true_matrix))
    
    if (length(common_samples) == 0) {
      cat("  Warning: No common samples between estimated and true matrices\n")
      next
    }
    
    # Align matrices
    estimated_matrix <- estimated_matrix[, common_samples, drop = FALSE]
    true_matrix <- true_matrix[, common_samples, drop = FALSE]
    
    # Calculate correlations for each sample
    sample_correlations <- numeric(length(common_samples))
    names(sample_correlations) <- common_samples
    
    for (i in seq_along(common_samples)) {
      sample_id <- common_samples[i]
      estimated_vector <- estimated_matrix[, sample_id]
      true_vector <- true_matrix[, sample_id]
      
      # Calculate Pearson correlation
      cor_result <- tryCatch({
        cor(estimated_vector, true_vector, method = "pearson", use = "complete.obs")
      }, error = function(e) {
        cat("  Warning: Correlation calculation failed for sample", sample_id, ":", e$message, "\n")
        NA_real_
      })
      
      sample_correlations[i] <- cor_result
    }
    
    # Store correlations in results data frame
    for (sample_id in names(sample_correlations)) {
      if (sample_id %in% rownames(results_df)) {
        results_df[sample_id, method] <- sample_correlations[sample_id]
      }
    }
    
    # Print summary
    non_na_cor <- sample_correlations[!is.na(sample_correlations)]
    if (length(non_na_cor) > 0) {
      cat("  Successfully calculated", length(non_na_cor), "/", length(sample_correlations), "correlations\n")
      cat("  Mean correlation:", round(mean(non_na_cor, na.rm = TRUE), 4), "\n")
    } else {
      cat("  All correlations are NA\n")
    }
    
  }, error = function(e) {
    cat("  Error in deconvolution for", method, ":", e$message, "\n")
  })
  
  cat(strrep("-", 60), "\n", sep = "")
}

# 9.3 Create final sample_results data frame
# ----------------------------------------------------------------------------
cat("\nCreating final results data frame...\n")

# Initialize sample_results with sample names
sample_results <- data.frame(
  sample = rownames(results_df),
  stringsAsFactors = FALSE
)

# Add method columns
for (method in methods_to_use) {
  sample_results[[method]] <- results_df[, method]
}

# Check for NA values
na_counts <- colSums(is.na(sample_results[, -1, drop = FALSE]))
cat("\nNA counts per method:\n")
for (method in methods_to_use) {
  cat(sprintf("  %-15s: %d / %d samples are NA\n", 
              method, na_counts[method], nrow(sample_results)))
}

# Replace NA with 0 (optional, as in original code)
# sample_results[is.na(sample_results)] <- 0

cat("\nFirst few rows of sample correlations:\n")
print(head(sample_results))


# ============================================================================
# SECTION 10: SINGLE-CELL DECONVOLUTION METHOD EVALUATION (MuSiC)
# ============================================================================

cat("\nEvaluating single-cell deconvolution method (MuSiC)...\n")

# Prepare single-cell data for MuSiC
sc_data_for_music <- train_data_cellID
phenoData_music <- pData_train

# Ensure appropriate sample identifier column
if (length(grep("[N-n]ame", colnames(phenoData_music))) > 0) {
  sample_col <- grep("[N-n]ame", colnames(phenoData_music))
} else {
  sample_col <- grep("[S-s]ample|[S-s]ubject", colnames(phenoData_music))
}
colnames(phenoData_music)[sample_col] <- "SubjectName"
rownames(phenoData_music) <- phenoData_music$cellID

# Ensure normalization_scC and normalization_scT exist
if (!exists("normalization_scC")) {
  normalization_scC <- "TMM"
}
if (!exists("normalization_scT")) {
  normalization_scT <- "TMM"
}

# Transform and normalize
sc_data_trans <- Transformation(sc_data_for_music, transformation)
bulk_expr_trans_music <- Transformation(bulk_expr_aligned, transformation)

sc_data_scaled <- Scaling(sc_data_trans, normalization_scC)
bulk_expr_scaled_music <- Scaling(bulk_expr_trans_music, normalization_scT)

# Align genes
common_genes_music <- intersect(rownames(sc_data_scaled), rownames(bulk_expr_scaled_music))
sc_data_aligned_music <- sc_data_scaled[common_genes_music, ]
bulk_expr_aligned_music <- bulk_expr_scaled_music[common_genes_music, ]

# Create SingleCellExperiment object
sce_object <- SingleCellExperiment(
  assays = list(counts = sc_data_aligned_music),
  colData = phenoData_music
)

# Run MuSiC deconvolution
set.seed(123456)
music_results <- t(MuSiC::music_prop(
  bulk.mtx = bulk_expr_aligned_music,
  sc.sce = sce_object,
  clusters = 'cellType',
  samples = 'cellID',
  markers = NULL,
  normalize = FALSE,
  verbose = FALSE
)$Est.prop.weighted)

# Calculate sample-wise correlations
estimated_music <- music_results[gtools::mixedsort(rownames(music_results)), ]
true_prop_music <- ground_truth_selected[gtools::mixedsort(rownames(ground_truth_selected)), ]

common_samples_music <- intersect(colnames(estimated_music), colnames(true_prop_music))
music_correlations <- sapply(common_samples_music, function(sample_id) {
  cor(estimated_music[, sample_id], true_prop_music[, sample_id], method = "pearson")
})

# Store MuSiC results
sample_results[["MuSiC"]] <- NA
sample_results[match(common_samples_music, sample_results$sample), "MuSiC"] <- music_correlations
# ============================================================================
# SECTION 11: CIBERSORTx DECONVOLUTION EVALUATION
# ============================================================================

cat("\nEvaluating CIBERSORTx deconvolution...\n")

# Load CIBERSORTx results file
cibersort_file <- "/media/desk16/tjn050/LSC/LSC_deconvolution/LSC/results/CIBERSORTx_wbl_Results.txt"
if (file.exists(cibersort_file)) {
  cibersort_results <- read.table(cibersort_file, header = TRUE, sep = "\t")
  
  # Extract proportion matrix
  n_cols <- ncol(cibersort_results)
  prop_matrix <- as.matrix(cibersort_results[, 2:(n_cols - 3)])
  rownames(prop_matrix) <- cibersort_results$Mixture
  prop_matrix <- t(prop_matrix)
  
  # Calculate sample-wise correlations
  common_samples_cibersort <- intersect(colnames(prop_matrix), colnames(ground_truth_selected))
  prop_matrix_aligned <- prop_matrix[, common_samples_cibersort, drop = FALSE]
  true_prop_aligned <- as.matrix(ground_truth_selected[, common_samples_cibersort, drop = FALSE])
  
  cibersort_correlations <- sapply(common_samples_cibersort, function(sample_id) {
    cor(prop_matrix_aligned[, sample_id], true_prop_aligned[, sample_id], method = "pearson")
  })
  
  # Store CIBERSORTx results
  sample_results[["CIBERSORTx"]] <- NA
  sample_results[match(common_samples_cibersort, sample_results$sample), "CIBERSORTx"] <- cibersort_correlations
} else {
  cat("  CIBERSORTx results file not found. Skipping CIBERSORTx evaluation.\n")
  sample_results[["CIBERSORTx"]] <- NA
}
# ============================================================================
# SECTION 12: ESMF DECONVOLUTION ANALYSIS
# ============================================================================

cat("\nPerforming ESMF deconvolution analysis...\n")

# 12.1 Prepare data for ESMF
# ----------------------------------------------------------------------------
# Ensure correct V matrix is used
# In original code, V was defined as:
# V <- as.matrix(T)[genes_filt, ]

# Use consistent V matrix construction as original code
P_selected <- ground_truth_selected
V_matrix <- as.matrix(bulk_expr_aligned)

# Ensure V is numeric matrix
V_matrix <- as.matrix(V_matrix)
if (!is.numeric(V_matrix)) {
  V_matrix <- apply(V_matrix, 2, as.numeric)
  rownames(V_matrix) <- rownames(bulk_expr_aligned)
}

# Align genes
genes_filt_esmf <- intersect(rownames(V_matrix), rownames(signature_matrix))
V_filt <- V_matrix[genes_filt_esmf, ]
signature_matrix_filt <- signature_matrix[genes_filt_esmf, ]
sd_matrix_filt <- sd_matrix[genes_filt_esmf, ]
final_markers_filt <- final_markers[rownames(signature_matrix_filt), ]

cat("  ESMF input dimensions:\n")
cat("    V matrix:", dim(V_filt)[1], "genes x", dim(V_filt)[2], "samples\n")
cat("    Signature matrix:", dim(signature_matrix_filt)[1], "genes x", dim(signature_matrix_filt)[2], "cell types\n")
cat("    Standard deviation matrix:", dim(sd_matrix_filt)[1], "genes x", dim(sd_matrix_filt)[2], "cell types\n")
cat("    Marker distribution:", nrow(final_markers_filt), "genes\n\n")

# 12.2 Check matrix types and dimensions
# ----------------------------------------------------------------------------
cat("  Checking matrix types and dimensions...\n")

# Check V matrix
cat("    V matrix type:", class(V_filt), "\n")
cat("    V matrix mode:", mode(V_filt), "\n")
cat("    V matrix contains NA:", any(is.na(V_filt)), "\n")
cat("    V matrix contains NaN:", any(is.nan(V_filt)), "\n")
cat("    V matrix contains Inf:", any(is.infinite(V_filt)), "\n")

# Check H matrix (to be initialized later)
cat("    Signature matrix type:", class(signature_matrix_filt), "\n")
cat("    Standard deviation matrix type:", class(sd_matrix_filt), "\n")

# 12.3 Generate coefficient matrices
# ----------------------------------------------------------------------------
m <- dim(V_filt)[1]  # Number of genes
n <- dim(V_filt)[2]  # Number of samples
r <- length(unique(final_markers_filt$CT))  # Number of cell types

cat("  Matrix dimensions:\n")
cat("    m (genes):", m, "\n")
cat("    n (samples):", n, "\n")
cat("    r (cell types):", r, "\n")

# Initialize random matrices
set.seed(1111)
W_init <- NMF::rmatrix(m, r)
set.seed(1112)
H_init <- NMF::rmatrix(r, n)

cat("  Initialized W matrix dimensions:", dim(W_init), "\n")
cat("  Initialized H matrix dimensions:", dim(H_init), "\n")

# Define penalty parameters
para1 <- 0.01  # General penalty
para2 <- 1     # Penalty for marker genes
para3 <- 1     # Penalty for priority marker genes

# Create coefficient matrix
coef_matrix <- matrix(data = para1, nrow = m, ncol = r)

# Apply enhanced penalty to marker genes
for (i in 1:length(unique(final_markers_filt$CT))) {
  gene_indices <- which(final_markers_filt$CT %in% names(table(final_markers_filt$CT))[i])
  coef_matrix[gene_indices, i] <- para2
}

# Apply highest penalty to priority marker genes
for (i in 1:length(unique(final_markers_filt$CT))) {
  cell_type_name <- names(table(final_markers_filt$CT))[i]
  if (cell_type_name %in% names(prio_markers)) {
    gene_indices <- which(final_markers_filt$gene %in% prio_markers[[cell_type_name]])
    if (length(gene_indices) > 0) {
      coef_matrix[gene_indices, i] <- para3
    }
  }
}

# 12.4 Check matrix multiplication compatibility
# ----------------------------------------------------------------------------
# Test V and H matrix multiplication compatibility
cat("  Testing matrix multiplication compatibility...\n")
cat("    V dimensions:", dim(V_filt), "\n")
cat("    t(H_init) dimensions:", dim(t(H_init)), "\n")

# Try computing V %*% t(H_init) to check for errors
test_result <- tryCatch({
  V_filt %*% t(H_init)
}, error = function(e) {
  cat("    Matrix multiplication error:", e$message, "\n")
  NULL
})

if (!is.null(test_result)) {
  cat("    Matrix multiplication successful. Result dimensions:", dim(test_result), "\n")
}

# 12.5 Execute ESMF deconvolution
# ----------------------------------------------------------------------------
cat("  Executing ESMF deconvolution...\n")

# Ensure correct function and parameters are used
# Original code used NMF_lsc_punish1 function
# Ensure function parameter names match original code

# Set sdfold parameter
sdfold <- 3
err_cutoff <- 0
iter_times <- 3000

# Create MarkerList object
ML_esmf <- CellMix::MarkerList()
ML_esmf@.Data <- tapply(as.character(final_markers_filt$gene), as.character(final_markers_filt$CT), list)

# Execute ESMF
set.seed(11)
esmf_results <- ESMF_deconvolution(
  V = V_filt,
  ml = ML_esmf,
  r = r,
  err_cutoff = err_cutoff,
  rate = coef_matrix,
  sig_matrix = signature_matrix_filt,
  sd_matrix = sd_matrix_filt,
  iter_times = iter_times
)

cat("  ESMF deconvolution completed.\n")

W_esmf <- esmf_results$W
H_esmf <- esmf_results$H

# 12.6 Process ESMF results
# ----------------------------------------------------------------------------
cat("  Processing ESMF results...\n")

# Calculate correlations
corW_esmf <- cor(W_esmf, signature_matrix_filt)
rownames(corW_esmf) <- 1:nrow(corW_esmf)

# Find optimal cell type alignment
find_optimal_permutation <- function(A, B) {
  n <- nrow(A)
  D <- matrix(NA, n, n)
  for (i in 1:n) {
    for (j in 1:n) {
      D[j, i] <- sum((B[j, ] - A[i, ])^2)
    }
  }
  permutation <- c(clue::solve_LSAP(D))
  return(list(A = A[permutation, ], permutation = permutation))
}

B_diag <- diag(1, nrow(corW_esmf))
alignment <- find_optimal_permutation(corW_esmf, B_diag)
cell_type_order <- alignment$permutation
corW_aligned <- alignment$A
rownames(corW_aligned) <- colnames(corW_aligned)

# Normalize H matrix column sums to 1
for (j in 1:n) {
  col_sum <- sum(H_esmf[, j])
  H_esmf[, j] <- H_esmf[, j] / col_sum
}

# Reorder H matrix according to cell type matching
H_esmf <- H_esmf[cell_type_order, ]
rownames(H_esmf) <- colnames(corW_aligned)
colnames(H_esmf) <- colnames(V_filt)

# 12.7 Calculate sample correlations
# ----------------------------------------------------------------------------
# Align true proportions
common_samples_esmf <- intersect(colnames(H_esmf), colnames(P_selected))
H_esmf_aligned <- H_esmf[, common_samples_esmf, drop = FALSE]
true_prop_esmf <- as.matrix(P_selected[, common_samples_esmf, drop = FALSE])

# Ensure matrix dimensions match
if (nrow(H_esmf_aligned) != nrow(true_prop_esmf)) {
  cat("  Warning: Row dimensions do not match. Transposing matrices...\n")
  H_esmf_aligned <- t(H_esmf_aligned)
}

# Calculate sample correlations
esmf_correlations <- sapply(common_samples_esmf, function(sample_id) {
  cor(H_esmf_aligned[, sample_id], true_prop_esmf[, sample_id], method = "pearson")
})

# Store ESMF results
sample_results[["ESMF"]] <- NA
sample_results[match(common_samples_esmf, sample_results$sample), "ESMF"] <- esmf_correlations

cat("  ESMF analysis completed successfully.\n")

# ============================================================================
# SECTION 12.8: SAVE DETAILED EPIC AND ESMF ESTIMATES
# ============================================================================

cat("\nSaving detailed EPIC and ESMF estimates...\n")

# 12.8.1 Re-calculate EPIC detailed estimates
# ----------------------------------------------------------------------------
cat("  Re-calculating EPIC detailed estimates...\n")

# Prepare EPIC required data
require(EPIC)
C_EPIC <- list()

# Get common genes
common_genes_epic <- intersect(rownames(signature_matrix_filt), rownames(V_filt))

# Prepare marker gene list
marker_distrib_epic <- marker_distrib[marker_distrib$gene %in% common_genes_epic, ]
markers_epic <- as.character(marker_distrib_epic$gene)

# Get common cell types
common_CTs <- intersect(colnames(signature_matrix_filt), colnames(ref_var_matrix))

# Build EPIC reference data
C_EPIC[["sigGenes"]] <- rownames(signature_matrix_filt[markers_epic, common_CTs])
C_EPIC[["refProfiles"]] <- as.matrix(signature_matrix_filt[markers_epic, common_CTs])
C_EPIC[["refProfiles.var"]] <- ref_var_matrix[markers_epic, common_CTs]

# Run EPIC
epic_detailed <- t(EPIC::EPIC(
  bulk = as.matrix(V_filt),
  reference = C_EPIC,
  withOtherCells = TRUE,
  scaleExprs = FALSE
)$cellFractions)

# Remove "otherCells" row
epic_detailed <- subset(epic_detailed, rownames(epic_detailed) != "otherCells")

cat("  EPIC detailed estimates calculated. Dimensions:", dim(epic_detailed), "\n")

# 12.8.2 Prepare ESMF detailed estimates
# ----------------------------------------------------------------------------
cat("  Preparing ESMF detailed estimates...\n")

# Use H_esmf_aligned calculated in SECTION 12
esmf_detailed <- H_esmf_aligned

# Ensure ESMF matrix row names match EPIC
if (nrow(esmf_detailed) != nrow(epic_detailed)) {
  cat("  Warning: EPIC and ESMF matrix dimensions do not match.\n")
  cat("  EPIC:", nrow(epic_detailed), "rows, ESMF:", nrow(esmf_detailed), "rows\n")
  
  # Try transpose
  if (nrow(esmf_detailed) == ncol(epic_detailed)) {
    esmf_detailed <- t(esmf_detailed)
    cat("  Transposed ESMF matrix to match dimensions.\n")
  }
}

# 12.8.3 Create combined EPIC_ESMF_RESULTS data frame
# ----------------------------------------------------------------------------
cat("  Creating combined EPIC_ESMF_RESULTS data frame...\n")

# Convert EPIC matrix to long format
epic_long <- reshape2::melt(
  as.matrix(epic_detailed),
  varnames = c("cellType", "sampleID"),
  value.name = "EPIC"
)
epic_long <- epic_long[order(epic_long$cellType, epic_long$sampleID), ]

# Convert ESMF matrix to long format
esmf_long <- reshape2::melt(
  as.matrix(esmf_detailed),
  varnames = c("cellType", "sampleID"),
  value.name = "ESMF"
)
esmf_long <- esmf_long[order(esmf_long$cellType, esmf_long$sampleID), ]

# Merge EPIC and ESMF
epic_esmf_combined <- merge(epic_long, esmf_long, by = c("cellType", "sampleID"))

# 12.8.4 Add ground truth values
# ----------------------------------------------------------------------------
cat("  Adding ground truth values...\n")

# Prepare ground truth matrix
true_matrix <- as.matrix(P_selected)

# Convert ground truth to long format
true_long <- reshape2::melt(
  true_matrix,
  varnames = c("cellType", "sampleID"),
  value.name = "expected_values"
)
true_long <- true_long[order(true_long$cellType, true_long$sampleID), ]

# Merge all data
EPIC_ESMF_RESULTS <- merge(epic_esmf_combined, true_long, by = c("cellType", "sampleID"))

# Reorder columns
EPIC_ESMF_RESULTS <- EPIC_ESMF_RESULTS[, c("cellType", "sampleID", "EPIC", "expected_values", "ESMF")]

cat("  EPIC_ESMF_RESULTS created. Dimensions:", dim(EPIC_ESMF_RESULTS), "\n")
cat("  First few rows:\n")
print(head(EPIC_ESMF_RESULTS))

# ============================================================================
# SECTION 12.8.5: SAVE EPIC_ESMF_RESULTS
# ============================================================================

cat("  Creating combined EPIC_ESMF_RESULTS data frame...\n")

# Convert EPIC matrix to long format
epic_long <- reshape2::melt(
  as.matrix(epic_detailed),
  varnames = c("cellType", "sampleID"),
  value.name = "EPIC"
)
epic_long <- epic_long[order(epic_long$cellType, epic_long$sampleID), ]

# Convert ESMF matrix to long format
esmf_long <- reshape2::melt(
  as.matrix(esmf_detailed),
  varnames = c("cellType", "sampleID"),
  value.name = "ESMF"
)
esmf_long <- esmf_long[order(esmf_long$cellType, esmf_long$sampleID), ]

# Merge EPIC and ESMF
epic_esmf_combined <- merge(epic_long, esmf_long, by = c("cellType", "sampleID"))

# 12.8.4 Add ground truth values
# ----------------------------------------------------------------------------
cat("  Adding ground truth values...\n")

# Prepare ground truth matrix
true_matrix <- as.matrix(P_selected)

# Convert ground truth to long format
true_long <- reshape2::melt(
  true_matrix,
  varnames = c("cellType", "sampleID"),
  value.name = "expected_values"
)
true_long <- true_long[order(true_long$cellType, true_long$sampleID), ]

# Merge all data
EPIC_ESMF_RESULTS <- merge(epic_esmf_combined, true_long, by = c("cellType", "sampleID"))

# Reorder columns
EPIC_ESMF_RESULTS <- EPIC_ESMF_RESULTS[, c("cellType", "sampleID", "EPIC", "expected_values", "ESMF")]

cat("  EPIC_ESMF_RESULTS created. Dimensions:", dim(EPIC_ESMF_RESULTS), "\n")
cat("  First few rows:\n")
print(head(EPIC_ESMF_RESULTS))

# 12.8.5 Create timestamp
# ----------------------------------------------------------------------------
current_timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")

# 12.8.6 Save EPIC_ESMF_RESULTS
# ----------------------------------------------------------------------------
epic_esmf_filename <- file.path(results_dir, paste0(analysis_prefix, "_EPIC_ESMF_detailed_", current_timestamp, ".csv"))
write.csv(EPIC_ESMF_RESULTS, file = epic_esmf_filename, row.names = FALSE)
cat("  Saved EPIC_ESMF_RESULTS:", epic_esmf_filename, "\n")

# ============================================================================
# SECTION 13: RESULTS SAVING
# ============================================================================

cat("\nSaving results...\n")

# Replace all NA values with 0 in sample_results
sample_results[is.na(sample_results)] <- 0
cat("  Replaced all NA values with 0 in sample_results\n")

# Save sample-wise correlation results
results_filename <- file.path(results_dir, paste0(analysis_prefix, "_sample_correlations_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv"))
write.csv(sample_results, file = results_filename, row.names = FALSE)
cat("Results saved to:", results_filename, "\n")
cat("Analysis completed.\n")
