# ============================================================================
# plot-tune-esmf-heatmap.R —— 全组合网格调参结果热力图（三数据集合一）
#
# 数据来源：三个已跑完的 tune_esmf_paras_*.txt 结果文件
# 图：横坐标 = para 组合（para1/para2/para3），纵坐标 = rate（8 档），
#     颜色 = 留出验证集 H_Pearson；三个数据集分面（facet）。
# 风格：对齐 LSC_plots_20260531.R（敦煌暖色调 + Arial + 600dpi cairo）
# ============================================================================

library(ggplot2)
library(reshape2)
library(dplyr)

# ---- 配置 ----
plot_width <- 180   # mm
plot_height <- 90   # mm（三子图横排）
output_dir <- "/media/desk16/tjn050/LSC/LSC_deconvolution/LSC/plot"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# 敦煌暖色调（对齐参考脚本）
warm_low  <- "#6B4226"
warm_mid  <- "#E5D0A3"
warm_high <- "#BF472C"
grid_col  <- "#8E735B"

# ---- 读三组数据 ----
datasets <- list(
  "Pancreas (single-ref)" = "/media/desk16/tjn050/LSC/LSC_deconvolution/LSC/results/Ablation_Study_GeneralizationFromSingleReference_Hm_Pancreas_ParameterSensitivity_v2/tune_esmf_paras_20260828_070251.txt",
  "Adipose (multi-source)" = "/media/desk16/tjn050/LSC/LSC_deconvolution/LSC/results/Ablation_Study_MultiSourceReferenceIntegration_Hm_Adipose_ParameterSensitivity_v2/tune_esmf_paras_20260828_122345.txt",
  "Brain (cross-platform)" = "/media/desk16/tjn050/LSC/LSC_deconvolution/LSC/results/Ablation_Study_crossPlatformConsistency_Mm_Brain_ParameterSensitivity_v2/tune_esmf_paras_20260828_162805.txt"
)

all_df <- do.call(rbind, lapply(names(datasets), function(ds) {
  d <- read.table(datasets[[ds]], header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  d$dataset <- ds
  d
}))

# ---- 构造 para 组合标签（横坐标）----
all_df$para_combo <- paste(all_df$para1, all_df$para2, all_df$para3, sep = "_")
all_df$para_combo <- factor(
  all_df$para_combo,
  levels = unique(all_df$para_combo[order(all_df$para1, all_df$para2, all_df$para3)])
)
all_df$rate <- factor(all_df$rate, levels = sort(unique(all_df$rate)))
all_df$dataset <- factor(all_df$dataset, levels = names(datasets))

# ---- 画热力图（三子图横排）----
p <- ggplot(all_df, aes(x = para_combo, y = rate, fill = H_Pearson)) +
  geom_tile(color = grid_col, linewidth = 0.2) +
  scale_fill_gradient2(
    low = warm_low, mid = warm_mid, high = warm_high,
    midpoint = 0.9, space = "Lab",
    limits = c(0.6, 1.0),
    name = "H_Pearson\n(held-out)"
  ) +
  facet_wrap(~ dataset, nrow = 1) +
  labs(x = "Parameter combination (para1_para2_para3)", y = "rate") +
  theme_classic(base_family = "Arial") +
  theme(
    text = element_text(family = "Arial", color = "black"),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5,
                               size = 4, face = "bold", family = "Arial"),
    axis.text.y = element_text(size = 8, face = "bold", family = "Arial", color = "#4A2C1A"),
    axis.title = element_text(size = 9, face = "bold", family = "Arial"),
    strip.text = element_text(size = 8, face = "bold", family = "Arial"),
    legend.title = element_text(size = 8, face = "bold", family = "Arial"),
    legend.text = element_text(size = 7, family = "Arial"),
    legend.key.height = unit(1.2, "cm"),
    panel.grid = element_blank(),
    plot.margin = margin(5, 5, 5, 5, "mm")
  )

# ---- 保存 ----
out_png <- file.path(output_dir, "TUNE_ESMF_GRID_HEATMAP_3datasets.png")
ggsave(
  filename = out_png,
  plot = p,
  device = png,
  type = "cairo",
  width = plot_width / 25.4,
  height = plot_height / 25.4,
  units = "in",
  dpi = 600,
  bg = "white"
)
cat("Heatmap saved to:", out_png, "\n")

# ---- 额外输出：三组默认参数 vs 最优参数对照表 ----
cat("\n=== 默认参数(1,1,1,1) vs 最优参数 ===\n")
for (ds in names(datasets)) {
  dd <- all_df[all_df$dataset == ds, ]
  default_val <- dd$H_Pearson[dd$rate == "1" & dd$para1 == 1 & dd$para2 == 1 & dd$para3 == 1]
  best_idx <- which.max(dd$H_Pearson)
  cat(sprintf("%s: 默认=%.4f, 最优=%.4f (rate=%s, para=%s,%s,%s)\n",
              ds, default_val, dd$H_Pearson[best_idx],
              dd$rate[best_idx], dd$para1[best_idx], dd$para2[best_idx], dd$para3[best_idx]))
}
