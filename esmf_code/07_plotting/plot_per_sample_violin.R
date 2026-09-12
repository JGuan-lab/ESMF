# ============================================================================
# plot-per-sample-violin.R
# 逐样本反卷积结果小提琴图（回应 R1-Major 9 / R2 逐样本解释）。
# 每个数据集一张图：上 = per-sample Pearson，下 = per-sample RMSE。共三张图。
# 配色与顺序对齐 LSC_plots_20260531.R 小提琴图（method_colors，ESMF 深绯红排第一）。
# ============================================================================

library(ggplot2)
library(dplyr)
library(patchwork)
library(showtext)

font_add("Arial",
         regular    = "/usr/share/fonts/truetype/msttcorefonts/Arial.ttf",
         bold       = "/usr/share/fonts/truetype/msttcorefonts/arialbd.ttf",
         italic     = "/usr/share/fonts/truetype/msttcorefonts/Arial_Italic.ttf",
         bolditalic = "/usr/share/fonts/truetype/msttcorefonts/arialbi.ttf")
showtext_auto()
showtext_opts(dpi = 600)

method_colors <- c(
  "ESMF"        = "#9E2A2B",
  "CIBERSORTx"  = "#2D4C7A",
  "MuSiC"       = "#C17B3C",
  "DCQ"         = "#006C5B",
  "DSA"         = "#4A5D34",
  "ElasticNet"  = "#3E606F",
  "EPIC"        = "#B3570B",
  "FARDEEP"     = "#6D3D14",
  "Lasso"       = "#7A1B1D",
  "NNLS"        = "#005C73",
  "OLS"         = "#2E5C4A",
  "Ridge"       = "#556B2F",
  "RLR"         = "#4B3D33")

method_display <- c(
  ESMF = "ESMF", CIBERSORTx = "CIBERSORTx", MuSiC = "MuSiC",
  DCQ = "DCQ", DSA = "DSA", elastic_net = "ElasticNet",
  EPIC = "EPIC", FARDEEP = "FARDEEP", lasso = "Lasso",
  nnls = "NNLS", OLS = "OLS", ridge = "Ridge", RLR = "RLR")

OUT_DIR <- "/media/desk16/tjn050/LSC/LSC_deconvolution/LSC/plot"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

datasets <- list(
  Brain = list(
    csv = "/media/desk16/tjn050/LSC/LSC_deconvolution/LSC/results/crossPlatformConsistency_Mm_Brain/per_sample_results/per_sample_long.csv",
    label = "Mouse Brain"),
  Pancreas = list(
    csv = "/media/desk16/tjn050/LSC/LSC_deconvolution/LSC/results/GeneralizationFromSingleReference_Hm_Pancreas/per_sample_results/per_sample_long.csv",
    label = "Human Pancreas"),
  Adipose = list(
    csv = "/media/desk16/tjn050/LSC/LSC_deconvolution/LSC/results/MultiSourceReferenceIntegration_Hm_Adipose/per_sample_results/per_sample_long.csv",
    label = "Human Adipose"))

make_violin <- function(d, ylab, ylim) {
  ggplot(d, aes(x = Method, y = value)) +
    geom_violin(aes(fill = Method), alpha = 0.85, trim = FALSE,
                linewidth = 0.4, color = "gray30") +
    geom_boxplot(width = 0.15, fill = "white", outlier.shape = NA,
                 linewidth = 0.3, color = "gray30") +
    geom_jitter(aes(color = Method), width = 0.15, size = 1.0, alpha = 0.5) +
    stat_summary(fun = mean, geom = "point", shape = 21, size = 2.5,
                 fill = "white", color = "black", stroke = 0.6) +
    scale_x_discrete(limits = names(method_colors), drop = FALSE) +
    scale_fill_manual(values = method_colors, drop = FALSE) +
    scale_color_manual(values = method_colors, guide = "none") +
    coord_cartesian(ylim = ylim) +
    labs(y = ylab) +
    theme_classic(base_size = 10) +
    theme(text = element_text(family = "Arial", color = "black"),
          axis.title.x = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1,
                                     size = 8, face = "bold"),
          axis.text.y = element_text(size = 9),
          axis.title.y = element_text(size = 10, face = "bold"),
          axis.line = element_line(linewidth = 0.4),
          axis.ticks = element_line(linewidth = 0.3))
}

for (nm in names(datasets)) {
  long_df <- read.csv(datasets[[nm]]$csv, stringsAsFactors = FALSE)
  long_df$Method <- factor(method_display[long_df$method],
                           levels = names(method_colors))

  pear_df <- long_df %>%
    group_by(Method, tissue) %>%
    summarise(value = cor(observed_values, expected_values), .groups = "drop") %>%
    filter(!is.na(value))

  rmse_df <- long_df %>%
    group_by(Method, tissue) %>%
    summarise(value = sqrt(mean((observed_values - expected_values)^2)),
              .groups = "drop")

  p_pear <- make_violin(pear_df, "Per-sample Pearson r", c(-1, 1))
  p_rmse <- make_violin(rmse_df, "Per-sample RMSE", NULL)

  combined <- (p_pear / p_rmse) +
    plot_layout(guides = "collect") &
    theme(legend.position = "right",
          legend.title = element_text(size = 9, face = "bold"),
          legend.text = element_text(size = 8))

  combined <- combined +
    plot_annotation(title = paste0(datasets[[nm]]$label,
                                   " — per-sample deconvolution"),
                    theme = theme(plot.title = element_text(
                      face = "bold", size = 13, family = "Arial", hjust = 0.5)))

  ggsave(file.path(OUT_DIR, paste0("per_sample_violin_", nm, ".pdf")),
         combined, width = 11, height = 9, device = "pdf")
  ggsave(file.path(OUT_DIR, paste0("per_sample_violin_", nm, ".png")),
         combined, width = 11, height = 9, dpi = 600, device = "png")
  cat("已保存：", nm, "\n")
}

cat("完成。输出目录：", OUT_DIR, "\n")
