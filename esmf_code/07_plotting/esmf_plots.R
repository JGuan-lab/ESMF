# ============================================================================
# Heatmap Visualization for Deconvolution Performance Evaluation
# 
# This script generates publication-quality heatmaps for evaluating 
# deconvolution methods using Pearson correlation and RMSE metrics.
# Features: Golden ratio layout, Dunhuang-inspired color palettes, 
#           Arial fonts, and high-resolution output for SCI papers.
# ============================================================================

# Load required libraries
library(ggplot2)
library(reshape2)

# ============================================================================
# CONFIGURATION SECTION
# ============================================================================

# Define golden ratio dimensions (mm) - optimized for publication
plot_width <- 180
plot_height <- plot_width / ((1 + sqrt(5)) / 2)  # Golden ratio: ~61.8mm

# Color palettes for different metrics
warm_palette <- list(
  low = "#6B4226",    
  mid = "#E5D0A3",    
  high = "#BF472C",   
  grid = "#8E735B",
  text_light = "#FFFFFF",
  text_dark = "#4A2C1A"
)

cool_palette <- list(
  low = "#3A5F73",    
  mid = "#8CA9B7",    
  high = "#E1EDF2",   
  grid = "#8E735B",
  text_light = "#FFFFFF", 
  text_dark = "#2B414D"
)

# Output directory
output_dir <- "/media/desk16/tjn050/LSC/LSC_deconvolution/LSC/plot"

# ============================================================================
# FUNCTION: create_heatmap
# ============================================================================
create_heatmap <- function(data_path, metric_type = c("Pearson", "RMSE")) {
  # Create output directory if needed
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Set parameters based on metric type
  metric_type <- match.arg(metric_type)
  if (metric_type == "Pearson") {
    palette <- warm_palette
    text_color <- palette$text_light
    color_scale_low <- palette$low
    color_scale_mid <- palette$mid
    color_scale_high <- palette$high
  } else {
    palette <- cool_palette
    text_color <- "#1A3447"
    color_scale_low <- palette$high
    color_scale_mid <- palette$mid
    color_scale_high <- palette$low
  }
  
  # Read and prepare data
  data <- read.csv(data_path, sep = ",", header = TRUE)
  methods <- data[, 1]
  data_matrix <- as.matrix(apply(data[, -1], 2, as.numeric))
  rownames(data_matrix) <- methods
  data_melt <- melt(data_matrix)
  
  # Generate output filename
  output_png <- file.path(
    output_dir,
    paste0("GOLDEN_HEATMAP_", metric_type, "_", 
           tools::file_path_sans_ext(basename(data_path)),
           "_", plot_width, "x", plot_height, "mm.png")
  )
  
  # Create heatmap
  heatmap_plot <- ggplot(data_melt, aes(x = Var2, y = Var1, fill = value)) +
    geom_tile(color = palette$grid, linewidth = 0.3) +
    geom_text(
      aes(label = sprintf("%.2f", value)),
      color = text_color,
      size = 3.5,
      family = "Arial",
      fontface = "bold"
    ) +
    scale_fill_gradient2(
      low = color_scale_low,
      mid = color_scale_mid,
      high = color_scale_high,
      midpoint = 0,
      space = "Lab",
      trans = "pseudo_log",
      guide = guide_colorbar(
        title = metric_type,
        title.position = "top",
        barwidth = unit(2, "mm"),
        barheight = unit(40, "mm"),
        direction = "vertical",
        title.theme = element_text(
          size = 8,
          family = "Arial",
          face = "bold",
          margin = margin(b = 2, unit = "mm")
        ),
        label.theme = element_text(
          size = 7,
          family = "Arial",
          margin = margin(t = 1, unit = "mm")
        )
      )
    ) +
    labs(x = NULL, y = NULL) +  # Remove axis titles
    theme(
      legend.position = "right",
      legend.key = element_rect(color = NA),
      axis.title.x = element_blank(),  # Remove x-axis title
      axis.text.x = element_text(
        angle = 45,
        hjust = 1,
        vjust = 1.05,  # Adjusted to 1.05 for slight upward movement
        family = "Arial",
        size = 9,
        face = "bold"  # MODIFIED: Bold x-axis labels
      ),
      axis.text.y = element_text(
        family = "Arial",
        size = 9,
        color = palette$text_dark,
        face = "bold"  # MODIFIED: Bold y-axis labels
      ),
      axis.ticks.y = element_line(color = palette$grid),
      axis.ticks.x = element_line(color = palette$grid),
      plot.margin = margin(5, 5, 5, 5, "mm")
    )
  
  # Save high-resolution plot
  ggsave(
    filename = output_png,
    plot = heatmap_plot,
    device = png,
    type = "cairo",
    width = plot_width / 25.4,
    height = plot_height / 25.4,
    units = "in",
    dpi = 600,
    bg = "white"
  )
  
  cat("Heatmap saved to:", output_png, "\n")
}

# ============================================================================
# BATCH PROCESSING FOR PEARSON FILES
# ============================================================================

pearson_files <- c(
  "/media/desk16/tjn050/LSC/deconvolution result/PLOTS/crossPlatformConsistency-PEARSON.csv",
  "/media/desk16/tjn050/LSC/deconvolution result/PLOTS/MultiSourceReferenceIntegration-PEARSON.csv",
  "/media/desk16/tjn050/LSC/deconvolution result/PLOTS/GeneralizationFromSingleReference-PEARSON.csv"
)

cat("Processing Pearson correlation heatmaps...\n")
for (file in pearson_files) {
  if (file.exists(file)) {
    cat("Processing:", basename(file), "...\n")
    create_heatmap(file, "Pearson")
  } else {
    warning("File not found: ", file)
  }
}
cat("Pearson heatmaps processing completed.\n")

# ============================================================================
# BATCH PROCESSING FOR RMSE FILES
# ============================================================================

rmse_files <- c(
  "/media/desk16/tjn050/LSC/deconvolution result/PLOTS/crossPlatformConsistency-RMSE.csv",
  "/media/desk16/tjn050/LSC/deconvolution result/PLOTS/MultiSourceReferenceIntegration-RMSE.csv",
  "/media/desk16/tjn050/LSC/deconvolution result/PLOTS/GeneralizationFromSingleReference-RMSE.csv"
)

cat("\nProcessing RMSE heatmaps...\n")
for (file in rmse_files) {
  if (file.exists(file)) {
    cat("Processing:", basename(file), "...\n")
    create_heatmap(file, "RMSE")
  } else {
    warning("File not found: ", file)
  }
}
cat("RMSE heatmaps processing completed.\n")
cat("\nAll heatmaps have been generated successfully!\n")

# ============================================================================
# DIRECT EXECUTION EXAMPLE (Uncomment to run)
# ============================================================================

# # Process a single file
input_path <- "/media/desk16/tjn050/LSC/deconvolution result/PLOTS/crossPlatformConsistency-PEARSON.csv"
if (file.exists(input_path)) {
  create_heatmap(input_path, "Pearson")
} else {
  cat("Input file does not exist. Please check the path.\n")
  cat("Current working directory:", getwd(), "\n")
  cat("Please select one of the following options:\n")
  cat("1. Pearson correlation: crossPlatformConsistency-PEARSON.csv, MultiSourceReferenceIntegration-PEARSON.csv, GeneralizationFromSingleReference-PEARSON.csv\n")
  cat("2. RMSE: crossPlatformConsistency-RMSE.csv, MultiSourceReferenceIntegration-RMSE.csv, GeneralizationFromSingleReference-RMSE.csv\n")
}




# ============================================================================
# VIOLIN PLOT FOR DECONVOLUTION METHODS COMPARISON
# ============================================================================
# Description: Creates publication-ready violin plot comparing deconvolution
#              method performance. ESMF is highlighted as the primary method.
#              The plot meets SCI journal requirements with harmonious
#              Dunhuang-inspired color palette and proper typography.
# Dependencies: ggplot2, dplyr, tidyr, ggsignif
# Output: High-resolution PNG (18×6 cm, 600 dpi) in publication format
# ============================================================================

# ============================================================================
# SCI-READY VIOLIN PLOT FOR DECONVOLUTION METHODS COMPARISON
# ============================================================================
# Description: Creates publication-ready violin plot comparing deconvolution
#              method performance. ESMF is highlighted as the primary method.
#              The plot meets SCI journal requirements with harmonious
#              Dunhuang-inspired color palette and proper typography.
# Dependencies: ggplot2, dplyr, tidyr, ggsignif
# Output: High-resolution PNG (18×6 cm, 600 dpi) in publication format
# ============================================================================

# LOAD REQUIRED LIBRARIES --------------------------------------------------
library(ggplot2)    # For plotting
library(dplyr)      # For data manipulation
library(tidyr)      # For data reshaping
library(ggsignif)   # For statistical significance bars

# CONFIGURATION ------------------------------------------------------------
input_path <- "/media/desk16/tjn050/LSC/deconvolution result/PLOTS/ClinicalValidationWithFlowCytometry_sample_correlations.csv"
output_dir <- "/media/desk16/tjn050/LSC/LSC_deconvolution/LSC/plot/"

# DUNHUANG-INSPIRED COLOR PALETTE FOR SCI PUBLICATIONS ---------------------
# Optimized for visual harmony and print compatibility
# ESMF: Primary method in distinctive red, others in balanced complementary colors
method_colors <- c(
  "ESMF"         = "#9E2A2B",  # Primary: Deep crimson (slightly desaturated for print)
  "CIBERSORTx"   = "#2D4C7A",  # Rich navy blue
  "MuSiC"        = "#C17B3C",  # Warm terracotta
  "DCQ"          = "#006C5B",  # Deep teal
  "DSA"          = "#4A5D34",  # Forest green
  "ElasticNet"   = "#3E606F",  # Slate blue
  "EPIC"         = "#B3570B",  # Burnt orange
  "FARDEEP"      = "#6D3D14",  # Chocolate brown
  "Lasso"        = "#7A1B1D",  # Maroon red
  "NNLS"         = "#005C73",  # Deep cyan
  "OLS"          = "#2E5C4A",  # Jade green
  "Ridge"        = "#556B2F",  # Olive green
  "RLR"          = "#4B3D33"   # Dark taupe
)

# DATA PREPROCESSING -------------------------------------------------------
# Read and transform correlation data
data_long <- read.csv(input_path) %>%
  # Convert wide format to long format
  pivot_longer(
    cols = -1,
    names_to = "Method",
    values_to = "Correlation"
  ) %>%
  # Filter valid correlation values
  filter(
    !is.na(Correlation),
    between(Correlation, -1, 1)
  ) %>%
  # Standardize method names
  mutate(
    Method = recode(
      Method,
      "elastic_net" = "ElasticNet",
      "nnls" = "NNLS",
      "ESMF" = "ESMF"
    ),
    # Order: ESMF first, then others in consistent order
    Method = factor(Method, levels = names(method_colors))
  )

# Calculate summary statistics
stats <- data_long %>% 
  group_by(Method) %>% 
  summarise(
    Mean = mean(Correlation, na.rm = TRUE),
    SD = sd(Correlation, na.rm = TRUE),
    .groups = "drop"
  )

# SCI PUBLICATION THEME ----------------------------------------------------
# Based on Nature/Science journal requirements
sci_theme <- theme_classic(base_size = 10) +
  theme(
    # Typography
    text = element_text(family = "Arial", color = "black"),
    plot.title = element_text(
      hjust = 0.5,
      face = "bold",
      size = 12,
      margin = margin(b = 8)
    ),
    plot.subtitle = element_text(
      hjust = 0.5,
      size = 9,
      color = "gray40",
      margin = margin(b = 12)
    ),
    
    # Axes
    axis.title.x = element_blank(),  # Remove x-axis title for cleaner look
    axis.title.y = element_text(
      size = 10,
      face = "bold",
      margin = margin(r = 8)
    ),
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      vjust = 1,
      size = 9,
      face = "bold"  # MODIFIED: Bold x-axis labels
    ),
    axis.text.y = element_text(
      size = 9,
      face = "bold"  # MODIFIED: Bold y-axis labels
    ),
    axis.line = element_line(linewidth = 0.4),
    axis.ticks = element_line(linewidth = 0.3),
    
    # Legend
    legend.position = "right",
    legend.title = element_text(size = 9, face = "bold"),
    legend.text = element_text(size = 8),
    legend.key.size = unit(0.7, "lines"),
    
    # Layout
    plot.margin = margin(12, 12, 12, 12)
  )

# CREATE VIOLIN PLOT -------------------------------------------------------
violin_plot <- ggplot(data_long, aes(x = Method)) + 
  # Violin density plot
  geom_violin(
    aes(y = Correlation, fill = Method),
    alpha = 0.85, 
    trim = FALSE,
    linewidth = 0.4,
    color = "gray30"
  ) +
  # Boxplot overlay
  geom_boxplot(
    aes(y = Correlation),
    width = 0.15,
    fill = "white",
    outlier.shape = NA,
    linewidth = 0.3
  ) +
  # Individual data points
  geom_jitter(
    aes(y = Correlation, color = Method),
    width = 0.15,
    size = 1.2,
    alpha = 0.5
  ) +
  # Mean ± SD error bars
  geom_errorbar(
    data = stats,
    aes(
      x = Method, 
      y = Mean,
      ymin = Mean - SD,
      ymax = Mean + SD
    ),
    width = 0.1,
    color = "black",
    linewidth = 0.6
  ) +
  # Mean point indicator
  geom_point(
    data = stats,
    aes(x = Method, y = Mean),
    size = 3,
    color = "white",
    fill = "black",
    shape = 21,
    stroke = 0.8
  ) +
  # Statistical significance annotation (ESMF vs EPIC)
  geom_signif(
    aes(y = Correlation),
    comparisons = list(c("ESMF", "EPIC")),
    map_signif_level = TRUE,
    y_position = 1.25,
    tip_length = 0.01,
    textsize = 3.2,
    vjust = -0.3
  ) +
  # Color scales
  scale_x_discrete(limits = names(method_colors)) +
  scale_fill_manual(values = method_colors, name = "Method") +
  scale_color_manual(values = method_colors, name = "Method", guide = "none") +
  # Axis limits
  coord_cartesian(ylim = c(-0.8, 1.3)) +
  # Labels
  labs(
    title = "ESMF Demonstrates Superior Correlation Accuracy",
    subtitle = "Comparison with EPIC Deconvolution Method",
    y = "Pearson Correlation Coefficient (r)"
  ) +
  # Apply SCI theme
  sci_theme

# DISPLAY PLOT -------------------------------------------------------------
print(violin_plot)

# SAVE PLOT FUNCTION -------------------------------------------------------
# Saves plot in SCI-compliant format
save_sci_plot <- function(
    plot_obj, 
    base_name = "violin_plot",
    output_path = output_dir,  # Use output_path parameter instead of default
    width_cm = 18,    # Standard 2-column width
    height_cm = 6,    # 3:1 aspect ratio
    dpi = 600
) {
  # Create filename with timestamp
  timestamp <- format(Sys.time(), "%Y%m%d")
  filename <- sprintf(
    "%s_%s_%ddpi.png",
    base_name,
    timestamp,
    dpi
  )
  
  full_path <- file.path(output_path, filename)
  
  # Save with publication settings
  ggsave(
    full_path,
    plot = plot_obj,
    device = "png",
    width = width_cm,
    height = height_cm,
    units = "cm",    # Use centimeters for SCI publications
    dpi = dpi,
    bg = "white"
  )
  
  # Confirmation message
  message(sprintf(
    "✓ Plot saved:\n  %s\n  Dimensions: %d × %d cm\n  Resolution: %d dpi",
    full_path,
    width_cm,
    height_cm,
    dpi
  ))
}

# SAVE THE PLOT ------------------------------------------------------------
save_sci_plot(violin_plot)

# ============================================================================
# Heatmap Visualization for Deconvolution Performance Evaluation
# 
# This script generates publication-quality heatmaps for evaluating 
# deconvolution methods using Pearson correlation and RMSE metrics.
# Features: Golden ratio layout, Dunhuang-inspired color palettes, 
#           Arial fonts, and high-resolution output for SCI papers.
# ============================================================================

# Load required libraries
library(ggplot2)
library(reshape2)

# ============================================================================
# CONFIGURATION SECTION
# ============================================================================

# Define golden ratio dimensions (mm) - optimized for publication
plot_width <- 180
plot_height <- plot_width / ((1 + sqrt(5)) / 2)  # Golden ratio: ~61.8mm

# Color palettes for different metrics
warm_palette <- list(
  low = "#6B4226",    
  mid = "#E5D0A3",    
  high = "#BF472C",   
  grid = "#8E735B",
  text_light = "#FFFFFF",
  text_dark = "#4A2C1A"
)

cool_palette <- list(
  low = "#3A5F73",    
  mid = "#8CA9B7",    
  high = "#E1EDF2",   
  grid = "#8E735B",
  text_light = "#FFFFFF", 
  text_dark = "#2B414D"
)

# Output directory
output_dir <- "/media/desk16/tjn050/LSC/LSC_deconvolution/LSC/plot"

# ============================================================================
# FUNCTION: create_heatmap
# ============================================================================
create_heatmap <- function(data_path, metric_type = c("Pearson", "RMSE")) {
  # Create output directory if needed
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Set parameters based on metric type
  metric_type <- match.arg(metric_type)
  if (metric_type == "Pearson") {
    palette <- warm_palette
    text_color <- palette$text_light
    color_scale_low <- palette$low
    color_scale_mid <- palette$mid
    color_scale_high <- palette$high
  } else {
    palette <- cool_palette
    text_color <- "#1A3447"
    color_scale_low <- palette$high
    color_scale_mid <- palette$mid
    color_scale_high <- palette$low
  }
  
  # Read and prepare data
  data <- read.csv(data_path, sep = ",", header = TRUE)
  methods <- data[, 1]
  data_matrix <- as.matrix(apply(data[, -1], 2, as.numeric))
  rownames(data_matrix) <- methods
  data_melt <- melt(data_matrix)
  
  # Generate output filename
  output_png <- file.path(
    output_dir,
    paste0("GOLDEN_HEATMAP_", metric_type, "_", 
           tools::file_path_sans_ext(basename(data_path)),
           "_", plot_width, "x", plot_height, "mm.png")
  )
  
  # Create heatmap
  heatmap_plot <- ggplot(data_melt, aes(x = Var2, y = Var1, fill = value)) +
    geom_tile(color = palette$grid, linewidth = 0.3) +
    geom_text(
      aes(label = sprintf("%.2f", value)),
      color = text_color,
      size = 3.5,
      family = "Arial",
      fontface = "bold"
    ) +
    scale_fill_gradient2(
      low = color_scale_low,
      mid = color_scale_mid,
      high = color_scale_high,
      midpoint = 0,
      space = "Lab",
      trans = "pseudo_log",
      guide = guide_colorbar(
        title = metric_type,
        title.position = "top",
        barwidth = unit(2, "mm"),
        barheight = unit(40, "mm"),
        direction = "vertical",
        title.theme = element_text(
          size = 8,
          family = "Arial",
          face = "bold",
          margin = margin(b = 2, unit = "mm")
        ),
        label.theme = element_text(
          size = 7,
          family = "Arial",
          margin = margin(t = 1, unit = "mm")
        )
      )
    ) +
    labs(x = NULL, y = NULL) +  # Remove axis titles
    theme(
      legend.position = "right",
      legend.key = element_rect(color = NA),
      axis.title.x = element_blank(),  # Remove x-axis title
      axis.text.x = element_text(
        angle = 45,
        hjust = 1,
        vjust = 1.05,  # Adjusted to 1.05 for slight upward movement
        family = "Arial",
        size = 9,
        face = "bold"  # MODIFIED: Bold x-axis labels
      ),
      axis.text.y = element_text(
        family = "Arial",
        size = 9,
        color = palette$text_dark,
        face = "bold"  # MODIFIED: Bold y-axis labels
      ),
      axis.ticks.y = element_line(color = palette$grid),
      axis.ticks.x = element_line(color = palette$grid),
      plot.margin = margin(5, 5, 5, 5, "mm")
    )
  
  # Save high-resolution plot
  ggsave(
    filename = output_png,
    plot = heatmap_plot,
    device = png,
    type = "cairo",
    width = plot_width / 25.4,
    height = plot_height / 25.4,
    units = "in",
    dpi = 600,
    bg = "white"
  )
  
  cat("Heatmap saved to:", output_png, "\n")
}

# ============================================================================
# BATCH PROCESSING FOR PEARSON FILES
# ============================================================================

pearson_files <- c(
  "/media/desk16/tjn050/LSC/deconvolution result/PLOTS/crossPlatformConsistency-PEARSON.csv",
  "/media/desk16/tjn050/LSC/deconvolution result/PLOTS/MultiSourceReferenceIntegration-PEARSON.csv",
  "/media/desk16/tjn050/LSC/deconvolution result/PLOTS/GeneralizationFromSingleReference-PEARSON.csv"
)

cat("Processing Pearson correlation heatmaps...\n")
for (file in pearson_files) {
  if (file.exists(file)) {
    cat("Processing:", basename(file), "...\n")
    create_heatmap(file, "Pearson")
  } else {
    warning("File not found: ", file)
  }
}
cat("Pearson heatmaps processing completed.\n")

# ============================================================================
# BATCH PROCESSING FOR RMSE FILES
# ============================================================================

rmse_files <- c(
  "/media/desk16/tjn050/LSC/deconvolution result/PLOTS/crossPlatformConsistency-RMSE.csv",
  "/media/desk16/tjn050/LSC/deconvolution result/PLOTS/MultiSourceReferenceIntegration-RMSE.csv",
  "/media/desk16/tjn050/LSC/deconvolution result/PLOTS/GeneralizationFromSingleReference-RMSE.csv"
)

cat("\nProcessing RMSE heatmaps...\n")
for (file in rmse_files) {
  if (file.exists(file)) {
    cat("Processing:", basename(file), "...\n")
    create_heatmap(file, "RMSE")
  } else {
    warning("File not found: ", file)
  }
}
cat("RMSE heatmaps processing completed.\n")
cat("\nAll heatmaps have been generated successfully!\n")

# ============================================================================
# DIRECT EXECUTION EXAMPLE (Uncomment to run)
# ============================================================================

# # Process a single file
input_path <- "/media/desk16/tjn050/LSC/deconvolution result/PLOTS/crossPlatformConsistency-PEARSON.csv"
if (file.exists(input_path)) {
  create_heatmap(input_path, "Pearson")
} else {
  cat("Input file does not exist. Please check the path.\n")
  cat("Current working directory:", getwd(), "\n")
  cat("Please select one of the following options:\n")
  cat("1. Pearson correlation: crossPlatformConsistency-PEARSON.csv, MultiSourceReferenceIntegration-PEARSON.csv, GeneralizationFromSingleReference-PEARSON.csv\n")
  cat("2. RMSE: crossPlatformConsistency-RMSE.csv, MultiSourceReferenceIntegration-RMSE.csv, GeneralizationFromSingleReference-RMSE.csv\n")
}

# ============================================================================
# VIOLIN PLOT FOR DECONVOLUTION METHODS COMPARISON
# ============================================================================
# Description: Creates publication-ready violin plot comparing deconvolution
#              method performance. ESMF is highlighted as the primary method.
#              The plot meets SCI journal requirements with harmonious
#              Dunhuang-inspired color palette and proper typography.
# Dependencies: ggplot2, dplyr, tidyr, ggsignif
# Output: High-resolution PNG (18×6 cm, 600 dpi) in publication format
# ============================================================================

# ============================================================================
# SCI-READY VIOLIN PLOT FOR DECONVOLUTION METHODS COMPARISON
# ============================================================================
# Description: Creates publication-ready violin plot comparing deconvolution
#              method performance. ESMF is highlighted as the primary method.
#              The plot meets SCI journal requirements with harmonious
#              Dunhuang-inspired color palette and proper typography.
# Dependencies: ggplot2, dplyr, tidyr, ggsignif
# Output: High-resolution PNG (18×6 cm, 600 dpi) in publication format
# ============================================================================

# LOAD REQUIRED LIBRARIES --------------------------------------------------
library(ggplot2)    # For plotting
library(dplyr)      # For data manipulation
library(tidyr)      # For data reshaping
library(ggsignif)   # For statistical significance bars

# CONFIGURATION ------------------------------------------------------------
input_path <- "/media/desk16/tjn050/LSC/LSC_deconvolution/LSC/wbl_cor_sample_groundtruth_20251028.csv"
output_dir <- "/media/desk16/tjn050/LSC/LSC_deconvolution/LSC/plot/"

# DUNHUANG-INSPIRED COLOR PALETTE FOR SCI PUBLICATIONS ---------------------
# Optimized for visual harmony and print compatibility
# ESMF: Primary method in distinctive red, others in balanced complementary colors
method_colors <- c(
  "ESMF"         = "#9E2A2B",  # Primary: Deep crimson (slightly desaturated for print)
  "CIBERSORTx"   = "#2D4C7A",  # Rich navy blue
  "MuSiC"        = "#C17B3C",  # Warm terracotta
  "DCQ"          = "#006C5B",  # Deep teal
  "DSA"          = "#4A5D34",  # Forest green
  "ElasticNet"   = "#3E606F",  # Slate blue
  "EPIC"         = "#B3570B",  # Burnt orange
  "FARDEEP"      = "#6D3D14",  # Chocolate brown
  "Lasso"        = "#7A1B1D",  # Maroon red
  "NNLS"         = "#005C73",  # Deep cyan
  "OLS"          = "#2E5C4A",  # Jade green
  "Ridge"        = "#556B2F",  # Olive green
  "RLR"          = "#4B3D33"   # Dark taupe
)

# DATA PREPROCESSING -------------------------------------------------------
# Read and transform correlation data
data_long <- read.csv(input_path) %>%
  # Convert wide format to long format
  pivot_longer(
    cols = -1,
    names_to = "Method",
    values_to = "Correlation"
  ) %>%
  # Filter valid correlation values
  filter(
    !is.na(Correlation),
    between(Correlation, -1, 1)
  ) %>%
  # Standardize method names
  mutate(
    Method = recode(
      Method,
      "elastic_net" = "ElasticNet",
      "nnls" = "NNLS",
      "ESMF" = "ESMF"
    ),
    # Order: ESMF first, then others in consistent order
    Method = factor(Method, levels = names(method_colors))
  )

# Calculate summary statistics
stats <- data_long %>% 
  group_by(Method) %>% 
  summarise(
    Mean = mean(Correlation, na.rm = TRUE),
    SD = sd(Correlation, na.rm = TRUE),
    .groups = "drop"
  )

# SCI PUBLICATION THEME ----------------------------------------------------
# Based on Nature/Science journal requirements
sci_theme <- theme_classic(base_size = 10) +
  theme(
    # Typography
    text = element_text(family = "Arial", color = "black"),
    plot.title = element_text(
      hjust = 0.5,
      face = "bold",
      size = 12,
      margin = margin(b = 8)
    ),
    plot.subtitle = element_text(
      hjust = 0.5,
      size = 9,
      color = "gray40",
      margin = margin(b = 12)
    ),
    
    # Axes
    axis.title.x = element_blank(),  # Remove x-axis title for cleaner look
    axis.title.y = element_text(
      size = 10,
      face = "bold",
      margin = margin(r = 8)
    ),
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      vjust = 1,
      size = 9,
      face = "bold"  # MODIFIED: Bold x-axis labels
    ),
    axis.text.y = element_text(
      size = 9,
      face = "bold"  # MODIFIED: Bold y-axis labels
    ),
    axis.line = element_line(linewidth = 0.4),
    axis.ticks = element_line(linewidth = 0.3),
    
    # Legend
    legend.position = "right",
    legend.title = element_text(size = 9, face = "bold"),
    legend.text = element_text(size = 8),
    legend.key.size = unit(0.7, "lines"),
    
    # Layout
    plot.margin = margin(12, 12, 12, 12)
  )

# CREATE VIOLIN PLOT -------------------------------------------------------
violin_plot <- ggplot(data_long, aes(x = Method)) + 
  # Violin density plot
  geom_violin(
    aes(y = Correlation, fill = Method),
    alpha = 0.85, 
    trim = FALSE,
    linewidth = 0.4,
    color = "gray30"
  ) +
  # Boxplot overlay
  geom_boxplot(
    aes(y = Correlation),
    width = 0.15,
    fill = "white",
    outlier.shape = NA,
    linewidth = 0.3
  ) +
  # Individual data points
  geom_jitter(
    aes(y = Correlation, color = Method),
    width = 0.15,
    size = 1.2,
    alpha = 0.5
  ) +
  # Mean ± SD error bars
  geom_errorbar(
    data = stats,
    aes(
      x = Method, 
      y = Mean,
      ymin = Mean - SD,
      ymax = Mean + SD
    ),
    width = 0.1,
    color = "black",
    linewidth = 0.6
  ) +
  # Mean point indicator
  geom_point(
    data = stats,
    aes(x = Method, y = Mean),
    size = 3,
    color = "white",
    fill = "black",
    shape = 21,
    stroke = 0.8
  ) +
  # Statistical significance annotation (ESMF vs EPIC)
  geom_signif(
    aes(y = Correlation),
    comparisons = list(c("ESMF", "EPIC")),
    map_signif_level = TRUE,
    y_position = 1.25,
    tip_length = 0.01,
    textsize = 3.2,
    vjust = -0.3
  ) +
  # Color scales
  scale_x_discrete(limits = names(method_colors)) +
  scale_fill_manual(values = method_colors, name = "Method") +
  scale_color_manual(values = method_colors, name = "Method", guide = "none") +
  # Axis limits
  coord_cartesian(ylim = c(-0.8, 1.3)) +
  # Labels
  labs(
    title = "ESMF Demonstrates Superior Correlation Accuracy",
    subtitle = "Comparison with EPIC Deconvolution Method",
    y = "Pearson Correlation Coefficient (r)"
  ) +
  # Apply SCI theme
  sci_theme

# DISPLAY PLOT -------------------------------------------------------------
print(violin_plot)

# SAVE PLOT FUNCTION -------------------------------------------------------
# Saves plot in SCI-compliant format
save_sci_plot <- function(
    plot_obj, 
    base_name = "violin_plot",
    output_path = output_dir,  # Use output_path parameter instead of default
    width_cm = 18,    # Standard 2-column width
    height_cm = 6,    # 3:1 aspect ratio
    dpi = 600
) {
  # Create filename with timestamp
  timestamp <- format(Sys.time(), "%Y%m%d")
  filename <- sprintf(
    "%s_%s_%ddpi.png",
    base_name,
    timestamp,
    dpi
  )
  
  full_path <- file.path(output_path, filename)
  
  # Save with publication settings
  ggsave(
    full_path,
    plot = plot_obj,
    device = "png",
    width = width_cm,
    height = height_cm,
    units = "cm",    # Use centimeters for SCI publications
    dpi = dpi,
    bg = "white"
  )
  
  # Confirmation message
  message(sprintf(
    "✓ Plot saved:\n  %s\n  Dimensions: %d × %d cm\n  Resolution: %d dpi",
    full_path,
    width_cm,
    height_cm,
    dpi
  ))
}

# SAVE THE PLOT ------------------------------------------------------------
save_sci_plot(violin_plot)

# ============================================================================
# Radar Charts for Comparison of EPIC and ESMF Methods
# ============================================================================
# Radar charts comparing deconvolution methods (EPIC vs ESMF) with ground truth
# Optimized for scientific publication with unified font styling
# Dependencies: fmsb, dplyr, tidyr
# ============================================================================

# Load data
data <- read.csv("/media/desk16/tjn050/LSC/deconvolution result/PLOTS/ClinicalValidationWithFlowCytometry_EPIC_ESMF.csv")

# Rename the last column to ESMF
colnames(data)[ncol(data)] <- "ESMF"

library(fmsb)
library(tidyverse)

# Dunhuang-inspired color palette optimized for SCI publications
dunhuang_colors <- c(
  "EPIC" = "#9A3B3B",        # Cinnabar red
  "ESMF" = "#2A5D3D",        # Emerald green
  "True Value" = "#4A708B"   # Deep blue-grey
)

# Calculate mean values for each cell type
avg_data <- data %>%
  group_by(cellType) %>%
  summarise(
    EPIC = mean(EPIC, na.rm = TRUE),
    ESMF = mean(ESMF, na.rm = TRUE),
    TrueValue = mean(expected_values, na.rm = TRUE)
  ) %>%
  ungroup()

# ====================== 1. EPIC vs True Value ======================
# Prepare EPIC data
epic_data <- avg_data %>% 
  select(cellType, EPIC) %>%
  pivot_wider(names_from = cellType, values_from = EPIC)

# Prepare true value data
true_data <- avg_data %>% 
  select(cellType, TrueValue) %>%
  pivot_wider(names_from = cellType, values_from = TrueValue)

# Create radar chart data frame
epic_radar <- rbind(
  rep(max(avg_data[, -1]) * 1.1, ncol(epic_data)),  # Upper limit
  rep(0, ncol(epic_data)),                          # Lower limit
  epic_data,
  true_data
)
rownames(epic_radar) <- c("max", "min", "EPIC", "True Value")

# Save EPIC radar chart
png("/media/desk16/tjn050/LSC/LSC_deconvolution/LSC/plot/EPIC_Radar.png", 
    width = 800, height = 800, res = 150)
par(mar = c(1, 1, 2, 1), family = "Arial", bg = "white")

# Draw radar chart
radarchart(
  epic_radar,
  axistype = 1,
  pcol = c(dunhuang_colors["EPIC"], dunhuang_colors["True Value"]),
  pfcol = c(adjustcolor(dunhuang_colors["EPIC"], alpha.f = 0.4), 
            adjustcolor(dunhuang_colors["True Value"], alpha.f = 0.4)),
  plwd = 2.5,
  cglcol = "grey50",
  cglty = 1,
  cglwd = 0.8,
  vlcex = 1.2,  # Increase axis label font size
  title = "EPIC vs True Value",
  axislabcol = "black",
  vfont = c("sans serif", "bold")  # MODIFIED: Use sans serif font with bold
)

# Add legend with bold font
legend("topright", 
       legend = c("EPIC", "True Value"),
       bty = "n", 
       pch = 15, 
       col = c(dunhuang_colors["EPIC"], dunhuang_colors["True Value"]),
       pt.cex = 1.8,
       cex = 1.2,  # Increase legend font size
       text.font = 2)  # Bold font (2 indicates bold)

dev.off()

# ====================== 2. ESMF vs True Value ======================
# Prepare ESMF data
esmf_data <- avg_data %>% 
  select(cellType, ESMF) %>%
  pivot_wider(names_from = cellType, values_from = ESMF)

# Create radar chart data frame
esmf_radar <- rbind(
  rep(max(avg_data[, -1]) * 1.1, ncol(esmf_data)),  # Upper limit
  rep(0, ncol(esmf_data)),                          # Lower limit
  esmf_data,
  true_data
)
rownames(esmf_radar) <- c("max", "min", "ESMF", "True Value")

# Save ESMF radar chart
png("/media/desk16/tjn050/LSC/LSC_deconvolution/LSC/plot/ESMF_Radar.png", 
    width = 800, height = 800, res = 150)
par(mar = c(1, 1, 2, 1), family = "Arial", bg = "white")

# Draw radar chart
radarchart(
  esmf_radar,
  axistype = 1,
  pcol = c(dunhuang_colors["ESMF"], dunhuang_colors["True Value"]),
  pfcol = c(adjustcolor(dunhuang_colors["ESMF"], alpha.f = 0.4), 
            adjustcolor(dunhuang_colors["True Value"], alpha.f = 0.4)),
  plwd = 2.5,
  cglcol = "grey50",
  cglty = 1,
  cglwd = 0.8,
  vlcex = 1.2,  # Increase axis label font size
  title = "ESMF vs True Value",
  axislabcol = "black",
  vfont = c("sans serif", "bold")  # MODIFIED: Use sans serif font with bold
)

# Add legend with bold font
legend("topright", 
       legend = c("ESMF", "True Value"),
       bty = "n", 
       pch = 15, 
       col = c(dunhuang_colors["ESMF"], dunhuang_colors["True Value"]),
       pt.cex = 1.8,
       cex = 1.2,  # Increase legend font size
       text.font = 2)  # Bold font (2 indicates bold)

dev.off()
# ============================================================================
# Script: Comparative Bar Plots for Deconvolution Metrics
# Purpose: Generate publication-quality grouped bar plots comparing two methods
#          (ESMF vs. MuSic) on six samples. The script processes two metrics:
#          1. Pearson Correlation (crossPlatformConsistency/MultiSourceReferenceIntegration/GeneralizationFromSingleReference-PEARSON.csv)
#          2. Root Mean Square Error (crossPlatformConsistency/MultiSourceReferenceIntegration/GeneralizationFromSingleReference-RMSE.csv)
#          Each metric plot includes paired t-test significance annotations.
# ============================================================================

# Load required libraries
library(ggplot2)
library(reshape2)
library(ggpubr)

# ============================================================================
# 1. FUNCTION TO PROCESS DATA AND GENERATE BAR PLOT
# ============================================================================
# Args:
#   file_path: Path to the input CSV file containing the data matrix
#   metric_name: Metric name to be used in the plot y-axis label ("Pearson Correlation" or "RMSE Value")
#   colors: Named vector of two colors for the two methods (ESMF and MuSic)
#   output_suffix: Suffix to append to the output file name
#
# Returns: A ggplot2 bar plot object

generate_bar_plot <- function(file_path, metric_name, colors, output_suffix) {
  
  # --- Data Loading and Preprocessing ---
  # Read the CSV file. Rownames are method names, columns are samples.
  data_mat <- read.csv(file_path, sep = ",", header = TRUE)
  methods <- data_mat[, 1]
  data_mat <- data_mat[, -1]
  data_mat <- apply(data_mat, 2, as.numeric)
  data_mat <- as.matrix(data_mat)
  rownames(data_mat) <- methods
  
  # Extract values for the two methods of interest
  # Assumes methods "MuSic" and "ESMF" are present in the data
  method1_vals <- data_mat["MuSic", ]
  method2_vals <- data_mat["ESMF", ]
  
  # Combine into a long-format dataframe for plotting
  plot_df <- rbind(
    data.frame(Sample = colnames(data_mat), Value = method1_vals, Method = "MuSic"),
    data.frame(Sample = colnames(data_mat), Value = method2_vals, Method = "ESMF")
  )
  # Ensure sample order is preserved and method order is fixed
  plot_df$Sample <- factor(plot_df$Sample, levels = colnames(data_mat))
  plot_df$Method <- factor(plot_df$Method, levels = c("MuSic", "ESMF"))
  
  n_samples <- ncol(data_mat)  # Number of samples (should be 6)
  
  # --- Statistical Comparison (Paired t-test) ---
  # Perform paired t-test between the two methods across all samples
  stat_test <- compare_means(
    Value ~ Method, 
    data = plot_df,
    method = "t.test",
    paired = TRUE
  )
  
  # Determine significance annotation symbol
  p_val <- stat_test$p
  sig_symbol <- ifelse(p_val < 0.001, "***",
                       ifelse(p_val < 0.01, "**",
                              ifelse(p_val < 0.05, "*", "ns")))
  
  # --- Plotting ---
  bar_plot <- ggplot(plot_df, aes(x = Sample, y = Value, fill = Method)) +
    geom_col(
      position = position_dodge(0.8),
      width = 0.7,
      color = "#5D4C46",  # Bar border color
      linewidth = 0.3
    ) +
    annotate("text", 
             x = (n_samples + 1) / 2,  # 居中放置
             y = max(plot_df$Value) * 1.15,  # 在最高条形图上方15%处
             label = sig_symbol,  # 显著性符号
             size = 4,  # 文本大小
             color = "#4A2C1A",  # 文本颜色
             fontface = "bold") +  # 加粗
    scale_fill_manual(values = colors) +
    labs(x = "", y = metric_name, fill = NULL) +
    theme_classic(base_family = "Arial") +  # Use Arial font
    theme(
      axis.title = element_text(size = 9, color = "black"),
      axis.text = element_text(size = 8, color = "black"),
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "top",
      plot.margin = margin(5, 5, 5, 5, "mm"),
      text = element_text(family = "Arial", color = "black")
    ) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.15)))  # Add 15% space on top
  
  return(bar_plot)
}

# ============================================================================
# 2. BATCH PROCESSING CONFIGURATION
# ============================================================================
# Define the base directory for input files
input_base_dir <- "/media/desk16/tjn050/LSC/deconvolution result/PLOTS"
# Define the output directory
output_base_dir <- "/media/desk16/tjn050/LSC/LSC_deconvolution/LSC/plot"

# Define the three dataset identifiers (crossPlatformConsistency, MultiSourceReferenceIntegration, GeneralizationFromSingleReference)
dataset_ids <- c("crossPlatformConsistency", "MultiSourceReferenceIntegration", "GeneralizationFromSingleReference")
# Define the two metrics
metrics <- c("PEARSON", "RMSE")

# Define color palettes for each metric
# Colors for Pearson correlation plot: MuSic (cool blue), ESMF (warm orange)
pearson_colors <- c("ESMF" = "#E67E22", "MuSic" = "#3A5F8E")
# Colors for RMSE plot: ESMF (green), MuSic (red)
rmse_colors <- c("ESMF" = "#2A5D3D", "MuSic" = "#9A3B3B")

# Plot dimensions (in mm) for publication
golden_width <- 86   # 8.6 cm
golden_height <- 65  # 6.5 cm

# ============================================================================
# 3. BATCH PROCESSING LOOP
# ============================================================================
# This loop iterates over all combinations of dataset identifiers and metrics,
# generating and saving a bar plot for each.

for (dset in dataset_ids) {
  for (met in metrics) {
    
    # Construct the input file path
    input_file <- file.path(input_base_dir, paste0(dset, "-", met, ".csv"))
    
    # Check if the file exists before processing
    if (!file.exists(input_file)) {
      warning(paste("File not found, skipping:", input_file))
      next
    }
    
    # Determine metric-specific parameters
    if (met == "PEARSON") {
      y_label <- "Pearson Correlation"
      colors <- pearson_colors
    } else {  # RMSE
      y_label <- "RMSE Value"
      colors <- rmse_colors
    }
    
    # Generate the plot
    plot_object <- generate_bar_plot(
      file_path = input_file,
      metric_name = y_label,
      colors = colors,
      output_suffix = paste0(dset, "_", met)
    )
    
    # Construct the output file name
    output_file_name <- paste0("GOLDEN_BAR_", dset, "_", met, ".png")
    output_path <- file.path(output_base_dir, output_file_name)
    
    # Save the plot with high-resolution settings
    ggsave(
      filename = output_path,
      plot = plot_object,
      device = png,
      type = "cairo",        # Anti-aliasing renderer
      width = golden_width,
      height = golden_height,
      units = "mm",          # Millimeter precision
      dpi = 600,             # Publication resolution
      bg = "white",          # White background
      family = "Arial"       # Embed Arial font
    )
    
    # Print progress message
    message(paste("Plot saved:", output_path))
  }
}

message("Batch processing complete.")

# ============================================================================
# Title: Line Plot Comparison of Deconvolution Method Performance

# Description: This script generates publication-quality line plots comparing
#              two deconvolution methods (ESMF vs. CIBERSORTx) across multiple
#              samples. The plots visualize Pearson correlation coefficients
#              and RMSE values with paired t-test significance markers.
#              Output: High-resolution PNG images formatted for SCI journals.
# Dependencies: ggplot2, ggpubr, reshape2
# ============================================================================

# Load required libraries
library(ggplot2)
library(ggpubr)
library(reshape2)

# ========================== PART 1: Global Settings ==========================
# Define output dimensions for publication
golden_width <- 85   # mm, standard width for single column
golden_height <- 70  # mm, adjusted for clarity

# Define output directory
output_dir <- "/media/desk16/tjn050/LSC/LSC_deconvolution/LSC/plot"

# Define unified plotting theme for SCI publications
sci_theme <- function(base_family = "Arial") {
  theme_classic(base_family = base_family) +
    theme(
      axis.title = element_text(size = 10, color = "black"),
      axis.text = element_text(size = 8, color = "black"),
      axis.text.x = element_text(angle = 35, hjust = 1, margin = margin(t = 5)),
      legend.position = "top",
      legend.title = element_text(size = 9, color = "black"),
      legend.text = element_text(size = 8, color = "black"),
      legend.spacing.x = unit(0.3, "cm"),
      panel.grid.major.y = element_line(color = "#E5D0A3", linewidth = 0.2),
      plot.margin = margin(7, 7, 7, 7, "mm")
    )
}

# ========================= PART 2: Function Definition =======================
# Function to generate comparison plots
generate_comparison_plot <- function(input_path, y_label) {
  # Function to read data, perform statistical test, and generate line plot
  
  # Read and preprocess data
  data_mat <- read.csv(input_path, sep = ",", header = TRUE)
  methods <- data_mat[, 1]
  data_mat <- data_mat[, -1]
  data_mat <- apply(data_mat, 2, as.numeric)
  data_mat <- as.matrix(data_mat)
  rownames(data_mat) <- methods
  
  # Extract method-specific values
  cibersort_vals <- data_mat["CIBERSORTx", ]
  esmf_vals <- data_mat["ESMF", ]
  
  # Prepare comparison data frame
  comparison_df <- rbind(
    data.frame(
      Sample = colnames(data_mat),
      Value = cibersort_vals,
      Method = "CIBERSORTx"
    ),
    data.frame(
      Sample = colnames(data_mat),
      Value = esmf_vals,
      Method = "ESMF"
    )
  )
  comparison_df$Sample <- factor(comparison_df$Sample, 
                                 levels = colnames(data_mat))
  comparison_df$Method <- factor(comparison_df$Method, 
                                 levels = c("CIBERSORTx", "ESMF"))
  
  # Perform paired t-test
  stat_test <- compare_means(
    Value ~ Method, 
    data = comparison_df,
    method = "t.test",
    paired = TRUE
  )
  
  # Prepare color palette
  if (grepl("PEARSON", input_path)) {
    color_palette <- c(
      "ESMF" = "#E67E22",      # Warm orange
      "CIBERSORTx" = "#3A5F8E"  # Deep blue
    )
  } else {
    color_palette <- c(
      "CIBERSORTx" = "#9A3B3B",  # Vermilion red
      "ESMF" = "#2A5D3D"         # Green-blue
    )
  }
  
  # Create significance annotation
  p_val <- stat_test$p
  sig_symbol <- ifelse(p_val < 0.001, "***",
                       ifelse(p_val < 0.01, "**",
                              ifelse(p_val < 0.05, "*", "ns")))
  
  # Calculate position for significance marker
  n_samples <- length(colnames(data_mat))
  x_center <- (n_samples + 1) / 2
  y_max <- max(comparison_df$Value)
  
  # Generate plot
  line_plot <- ggplot(comparison_df, 
                      aes(x = Sample, y = Value, 
                          color = Method, group = Method)) +
    geom_line(linewidth = 0.8, aes(linetype = Method), key_glyph = "path") +
    geom_point(size = 2.5, fill = "white", stroke = 0.8, shape = 21)
  
  # 只在p<0.05时添加显著性标记
  if (p_val < 0.05) {
    line_plot <- line_plot +
      annotate("text", x = x_center, y = y_max * 1.12, 
               label = sig_symbol, size = 4, color = "black")
  }
  
  # 添加其他图形元素
  line_plot <- line_plot +
    scale_color_manual(values = color_palette) +
    scale_linetype_manual(values = c("solid", "dashed")) +
    labs(x = "", y = y_label, color = "Method", linetype = "Method") +
    sci_theme() +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.18)))
  
  return(line_plot)
}

# ========================= PART 3: Plot Generation ==========================
# Define file paths and labels
pearson_files <- c(
  "/media/desk16/tjn050/LSC/deconvolution result/PLOTS/crossPlatformConsistency-PEARSON.csv",
  "/media/desk16/tjn050/LSC/deconvolution result/PLOTS/MultiSourceReferenceIntegration-PEARSON.csv",
  "/media/desk16/tjn050/LSC/deconvolution result/PLOTS/GeneralizationFromSingleReference-PEARSON.csv"
)

rmse_files <- c(
  "/media/desk16/tjn050/LSC/deconvolution result/PLOTS/crossPlatformConsistency-RMSE.csv",
  "/media/desk16/tjn050/LSC/deconvolution result/PLOTS/MultiSourceReferenceIntegration-RMSE.csv",
  "/media/desk16/tjn050/LSC/deconvolution result/PLOTS/GeneralizationFromSingleReference-RMSE.csv"
)

# Generate and save Pearson correlation plots
for (input_path in pearson_files) {
  plot_object <- generate_comparison_plot(input_path, "Pearson Correlation")
  output_png <- file.path(
    output_dir,
    paste0("GOLDEN_LINE_", tools::file_path_sans_ext(basename(input_path)),
           "_", golden_width, "x", golden_height, "mm.png")
  )
  
  ggsave(output_png,
         plot = plot_object,
         device = png,
         type = "cairo",
         width = golden_width,
         height = golden_height,
         units = "mm",
         dpi = 600,
         bg = "white",
         family = "Arial")
}

# Generate and save RMSE plots
for (input_path in rmse_files) {
  plot_object <- generate_comparison_plot(input_path, "RMSE Value")
  output_png <- file.path(
    output_dir,
    paste0("GOLDEN_LINE_", tools::file_path_sans_ext(basename(input_path)),
           "_", golden_width, "x", golden_height, "mm.png")
  )
  
  ggsave(output_png,
         plot = plot_object,
         device = png,
         type = "cairo",
         width = golden_width,
         height = golden_height,
         units = "mm",
         dpi = 600,
         bg = "white",
         family = "Arial")
}

# ============================================================================
# Script: Polar Bar Charts (Radial Bar Plots) for Deconvolution Metrics
# ============================================================================
# Purpose: Generate publication-quality radial bar charts (rose-like plots) 
#          comparing methods (ESMF vs. DSA) across multiple tissue samples.
#          The script processes two metrics:
#          1. Pearson Correlation (GeneralizationFromSingleReference-PEARSON-W.csv)
#          2. Root Mean Square Error (GeneralizationFromSingleReference-RMSE-W.csv)
#          Each metric generates a separate plot with consistent styling.
# Dependencies: tidyverse, ggplot2
# Output: High-resolution PNG (12×12 in, 600 dpi) in square format
# ============================================================================

# Load required libraries
library(tidyverse)
library(ggplot2)

# ============================================================================
# 1. FUNCTION TO GENERATE POLAR BAR CHART
# ============================================================================
# Args:
#   file_path: Path to the input CSV file
#   metric_type: Type of metric ("PEARSON" or "RMSE") for styling
#   output_suffix: Suffix for output file name
#
# Returns: A ggplot2 polar bar chart object

generate_polar_plot <- function(file_path, metric_type, output_suffix) {
  
  # --- Data Loading and Preprocessing ---
  analysis_data <- read_csv(file_path, show_col_types = FALSE) %>%
    pivot_longer(-method, names_to = "tissue", values_to = "value")
  
  # Clean and format tissue names
  analysis_data <- analysis_data %>%
    mutate(
      tissue = tissue %>% 
        str_replace_all("[-_.]+", "_") %>%  # Unify separators
        str_replace("Alipose", "Adipose") %>%  # Fix spelling
        str_replace("(Hm|Mm)_(\\w+)_(\\w+)", "\\1_\\2\n\\3")  # Add line break
    )
  
  # Set factor levels for methods
  analysis_data$method <- factor(analysis_data$method, 
                                 levels = c("ESMF", "DSA"))
  
  # Get unique tissue order and create a reversed order for polar coordinates
  # Keep the first tissue, reverse the rest
  tissue_order <- unique(analysis_data$tissue)
  new_order <- c(tissue_order[1], rev(tissue_order[-1]))
  analysis_data$tissue <- factor(analysis_data$tissue, levels = new_order)
  
  # --- Special processing for RMSE metric ---
  if (metric_type == "RMSE") {
    # Define scaling parameters for visualization enhancement
    method_ratios <- tribble(
      ~method,   ~base_scale, ~high_boost, ~low_boost,
      "ESMF",    1,           1,           1,
      "DSA",     1,           1,           1
    )
    
    # Apply scaling transformations
    analysis_data <- analysis_data %>%
      left_join(method_ratios, by = "method") %>%
      group_by(method) %>%
      mutate(
        base_scaled = value * base_scale,
        value_threshold = quantile(base_scaled, 0.65),
        scaled_value = case_when(
          base_scaled > value_threshold ~ 
            (base_scaled - value_threshold) * high_boost + value_threshold,
          base_scaled <= value_threshold ~ 
            (base_scaled / value_threshold) * (value_threshold * low_boost)
        )
      ) %>%
      ungroup() %>%
      mutate(scaled_value = pmin(scaled_value, max(base_scaled) * 1.2))
    
    # Use scaled value for RMSE plots
    y_var <- "scaled_value"
  } else {
    # For Pearson correlation, use original values
    y_var <- "value"
  }
  
  # --- Define color palettes based on metric type ---
  if (metric_type == "PEARSON") {
    # Warm gold (ESMF) and cool blue (DSA) for Pearson
    method_palette <- c("ESMF" = "#E9A85D", "DSA" = "#3A5E7A")
  } else {  # RMSE
    # Green (ESMF) and red (DSA) for RMSE
    method_palette <- c("ESMF" = "#2E5E4B", "DSA" = "#9A3B3B")
  }
  
  # --- Create plot with consistent styling ---
  polar_plot <- ggplot(analysis_data, 
                       aes(x = tissue, 
                           y = .data[[y_var]],  # Dynamic y variable
                           fill = method)) +
    geom_col(
      position = position_dodge(width = 0.9),
      width = 0.85,
      color = "white",  # White border for bars
      linewidth = 0.3,
      alpha = 0.95
    ) +
    # Add value labels for bars with value > 0.3
    geom_text(
      aes(label = ifelse(.data[[y_var]] > 0.3, round(.data[[y_var]], 2), "")),
      position = position_dodge(width = 0.9),
      angle = 45,  # 45-degree angle for better readability
      vjust = -0.8,
      hjust = 0.5,
      size = 10,  # Font size for labels
      color = "black",
      fontface = "plain"
    ) +
    # Convert to polar coordinates (radial/rose plot)
    coord_polar(theta = "x", start = 0) +
    # Manual fill colors
    scale_fill_manual(
      values = method_palette,
      guide = guide_legend(
        title = NULL,
        reverse = TRUE,  # Reverse legend order
        direction = "horizontal",
        nrow = 1,
        override.aes = list(color = NA, size = 1.5),
        keywidth = unit(1.2, "cm")
      )
    ) +
    # For RMSE: add log-scale y-axis
    {if(metric_type == "RMSE") 
      scale_y_continuous(
        trans = "log10",
        breaks = scales::trans_breaks("log10", function(x) 10^x),
        labels = scales::label_number()
      )
    } +
    # Consistent theme settings
    theme_minimal() +
    theme(
      axis.text.x = element_text(
        size = 25,
        angle = 90,  # Vertical text for radial segments
        hjust = 0.5,
        vjust = 0.5,
        lineheight = 0.8,
        family = "Arial",  # Use Arial font
        face = "bold"      # Make axis text bold (consistent with other plots)
      ),
      axis.title = element_blank(),
      axis.text.y = element_blank(),  # Remove y-axis labels (removes the scale)
      axis.ticks.y = element_blank(),  # Remove y-axis ticks
      # Legend positioning and styling
      legend.position = "top",
      legend.justification = "center",
      legend.background = element_rect(
        fill = alpha("white", 0.9), 
        color = NA
      ),
      legend.key = element_rect(fill = "white", color = NA),
      legend.text = element_text(
        family = "Arial",
        size = 30,
        color = "black",
        face = "bold"      # Make legend text bold (consistent with other plots)
      ),
      # Remove grid lines
      panel.grid = element_blank(),
      # Plot background
      plot.background = element_rect(fill = "white", color = NA),
      # Text elements
      text = element_text(family = "Arial", color = "black", face = "bold")  # Apply bold font globally
    )
  
  return(polar_plot)
}

# ============================================================================
# 2. BATCH PROCESSING CONFIGURATION
# ============================================================================
# Define base directories
input_base_dir <- "/media/desk16/tjn050/LSC/deconvolution result/PLOTS"
output_base_dir <- "/media/desk16/tjn050/LSC/LSC_deconvolution/LSC/plot"

# Square plot dimensions (inches) for consistent aspect ratio
plot_width <- 12
plot_height <- 12  # Square aspect ratio

# ============================================================================
# 3. BATCH PROCESSING LOOP
# ============================================================================
# Define all dataset identifiers to be processed
dataset_ids <- c("crossPlatformConsistency", "MultiSourceReferenceIntegration", "GeneralizationFromSingleReference")
# Define the two metric types
metrics <- c("PEARSON", "RMSE")

# Square plot dimensions (inches) - maintaining consistent aspect ratio
plot_width <- 12
plot_height <- 12  # Square aspect ratio

# Batch process all datasets and all metrics
for (dset in dataset_ids) {
  for (met in metrics) {
    # Construct file names (e.g., crossPlatformConsistency-PEARSON-W.csv, MultiSourceReferenceIntegration-PEARSON-W.csv, GeneralizationFromSingleReference-PEARSON-W.csv)
    file_suffix <- ifelse(met == "PEARSON", "PEARSON-W", "RMSE-W")
    input_file <- file.path(input_base_dir, paste0(dset, "-", file_suffix, ".csv"))
    
    # Check if the file exists
    if (!file.exists(input_file)) {
      warning(paste("File not found, skipping:", input_file))
      next
    }
    
    # Generate the polar bar plot
    plot_object <- generate_polar_plot(
      file_path = input_file,
      metric_type = met,
      output_suffix = paste0(dset, "_", met)
    )
    
    # Construct output file name
    output_file_name <- paste0("POLAR_BAR_", dset, "_", met, ".png")
    output_path <- file.path(output_base_dir, output_file_name)
    
    # Save the plot with high-resolution settings
    ggsave(
      filename = output_path,
      plot = plot_object,
      device = png,
      type = "cairo",        # Anti-aliasing renderer
      width = plot_width,
      height = plot_height,
      units = "in",          # Inches for square dimensions
      dpi = 600,             # Publication resolution
      bg = "white",          # White background
      family = "Arial"       # Embed Arial font
    )
    
    # Print progress message
    message(paste("Polar bar plot saved:", output_path))
  }
}

# Processing completion statistics
processed_count <- length(dataset_ids) * length(metrics)
message(paste("Batch processing completed. Processed", processed_count, "dataset-metric combinations."))
message(paste("Datasets:", paste(dataset_ids, collapse = ", ")))
message(paste("Metrics:", paste(metrics, collapse = ", ")))