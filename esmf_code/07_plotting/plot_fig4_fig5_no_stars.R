# ============================================================================
# plot-fig4-fig5-no-stars.R
#
# Purpose : R1-Major 7 —— 删除 Figure 4、Figure 5 中的显著性星号，重新出图。
#
# Source  : 逻辑取自 ESMF_codes/LSC_plots_20260531.R
#             Figure 4 = 柱状图，ESMF vs MuSic      （原 generate_bar_plot）
#             Figure 5 = 折线图，ESMF vs CIBERSORTx （原 generate_comparison_plot）
#           本脚本不修改、不读取原脚本，仅复刻其绘图参数。
#
# Change  : 唯一改动是移除星号标注，即原脚本中的
#             annotate("text", label = sig_symbol, ...)
#           以及为其服务的 compare_means(paired t-test)。
#           配色、主题、尺寸、dpi、字体均与原脚本保持一致。
#
# Reason  : n = 6，检验效力有限；两图本身已是逐样本量化对比，读者可直接从
#           柱高/折线差异判断优劣，星号属冗余（见沟通文档 R1-Major 7）。
#
# Output  : plot/GOLDEN_BAR_{dataset}_{metric}_nostars.png    (Figure 4)
#           plot/GOLDEN_LINE_{dataset}_{metric}_nostars.png   (Figure 5)
#           原图（无 _nostars 后缀）保持不动，便于对照。
# ============================================================================

library(ggplot2)

INPUT_DIR  <- "/media/desk16/tjn050/LSC/deconvolution result/PLOTS"
OUTPUT_DIR <- "/media/desk16/tjn050/LSC/LSC_deconvolution/LSC/plot"

DATASETS <- c("crossPlatformConsistency",
              "MultiSourceReferenceIntegration",
              "GeneralizationFromSingleReference")
METRICS  <- c("PEARSON", "RMSE")

# ---- 配色（与原脚本一致）----
bar_pearson_colors <- c("ESMF" = "#E67E22", "MuSic" = "#3A5F8E")
bar_rmse_colors    <- c("ESMF" = "#2A5D3D", "MuSic" = "#9A3B3B")

line_pearson_colors <- c("ESMF" = "#E67E22", "CIBERSORTx" = "#3A5F8E")
line_rmse_colors    <- c("CIBERSORTx" = "#9A3B3B", "ESMF" = "#2A5D3D")

# ---- 读取数据：第一列为方法名，其余列为样本 ----
read_matrix <- function(path) {
  d  <- read.csv(path, header = TRUE, check.names = FALSE)
  rn <- d[[1]]
  m  <- as.matrix(apply(d[, -1, drop = FALSE], 2, as.numeric))
  rownames(m) <- rn
  m
}

# ============================================================================
# Figure 4：柱状图，ESMF vs MuSic（无星号）
# ============================================================================
make_bar_plot <- function(data_mat, y_label, colors) {
  plot_df <- rbind(
    data.frame(Sample = colnames(data_mat), Value = data_mat["MuSic", ], Method = "MuSic"),
    data.frame(Sample = colnames(data_mat), Value = data_mat["ESMF",  ], Method = "ESMF")
  )
  plot_df$Sample <- factor(plot_df$Sample, levels = colnames(data_mat))
  plot_df$Method <- factor(plot_df$Method, levels = c("MuSic", "ESMF"))

  ggplot(plot_df, aes(x = Sample, y = Value, fill = Method)) +
    geom_col(position = position_dodge(0.8), width = 0.7,
             color = "#5D4C46", linewidth = 0.3) +
    scale_fill_manual(values = colors) +
    labs(x = "", y = y_label, fill = NULL) +
    theme_classic(base_family = "Arial") +
    theme(
      axis.title   = element_text(size = 9, color = "black"),
      axis.text    = element_text(size = 8, color = "black"),
      axis.text.x  = element_text(angle = 45, hjust = 1),
      legend.position = "top",
      plot.margin  = margin(5, 5, 5, 5, "mm"),
      text         = element_text(family = "Arial", color = "black")
    ) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.15)))
}

# ============================================================================
# Figure 5：折线图，ESMF vs CIBERSORTx（无星号）
# ============================================================================
sci_theme <- function(base_family = "Arial") {
  theme_classic(base_family = base_family) +
    theme(
      axis.title   = element_text(size = 10, color = "black"),
      axis.text    = element_text(size = 8, color = "black"),
      axis.text.x  = element_text(angle = 35, hjust = 1, margin = margin(t = 5)),
      legend.position = "top",
      legend.title = element_text(size = 9, color = "black"),
      legend.text  = element_text(size = 8, color = "black"),
      legend.spacing.x = unit(0.3, "cm"),
      panel.grid.major.y = element_line(color = "#E5D0A3", linewidth = 0.2),
      plot.margin  = margin(7, 7, 7, 7, "mm")
    )
}

make_line_plot <- function(data_mat, y_label, colors) {
  df <- rbind(
    data.frame(Sample = colnames(data_mat), Value = data_mat["CIBERSORTx", ], Method = "CIBERSORTx"),
    data.frame(Sample = colnames(data_mat), Value = data_mat["ESMF",       ], Method = "ESMF")
  )
  df$Sample <- factor(df$Sample, levels = colnames(data_mat))
  df$Method <- factor(df$Method, levels = c("CIBERSORTx", "ESMF"))

  ggplot(df, aes(x = Sample, y = Value, color = Method, group = Method)) +
    geom_line(linewidth = 0.8, aes(linetype = Method), key_glyph = "path") +
    geom_point(size = 2.5, fill = "white", stroke = 0.8, shape = 21) +
    scale_color_manual(values = colors) +
    scale_linetype_manual(values = c("solid", "dashed")) +
    labs(x = "", y = y_label, color = "Method", linetype = "Method") +
    sci_theme() +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.18)))
}

# ============================================================================
# 批量输出
# ============================================================================
if (!dir.exists(OUTPUT_DIR)) dir.create(OUTPUT_DIR, recursive = TRUE)

for (dset in DATASETS) {
  for (met in METRICS) {

    input_file <- file.path(INPUT_DIR, paste0(dset, "-", met, ".csv"))
    if (!file.exists(input_file)) {
      warning("File not found, skipping: ", input_file)
      next
    }

    data_mat <- read_matrix(input_file)

    if (met == "PEARSON") {
      y_label     <- "Pearson Correlation"
      bar_colors  <- bar_pearson_colors
      line_colors <- line_pearson_colors
    } else {
      y_label     <- "RMSE Value"
      bar_colors  <- bar_rmse_colors
      line_colors <- line_rmse_colors
    }

    # Figure 4（86 x 65 mm，与原脚本一致）
    p_bar <- make_bar_plot(data_mat, y_label, bar_colors)
    out_bar <- file.path(OUTPUT_DIR,
                         paste0("GOLDEN_BAR_", dset, "_", met, "_nostars.png"))
    ggsave(filename = out_bar, plot = p_bar, device = png, type = "cairo",
           width = 86, height = 65, units = "mm",
           dpi = 600, bg = "white", family = "Arial")

    # Figure 5（85 x 70 mm，与原脚本一致）
    p_line <- make_line_plot(data_mat, y_label, line_colors)
    out_line <- file.path(OUTPUT_DIR,
                          paste0("GOLDEN_LINE_", dset, "_", met, "_nostars.png"))
    ggsave(filename = out_line, plot = p_line, device = png, type = "cairo",
           width = 85, height = 70, units = "mm",
           dpi = 600, bg = "white", family = "Arial")

    message("Saved: ", basename(out_bar), "  |  ", basename(out_line))
  }
}

message("Done — no significance annotations are drawn.")
