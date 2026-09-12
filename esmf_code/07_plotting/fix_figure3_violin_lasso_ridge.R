# ============================================================================
# fix-figure3-violin-lasso-ridge.R
# 目的：单独修复 Figure 3（小提琴图）中 Lasso / Ridge 缺失的问题。
#   根因：原脚本 LSC_plots_20260531.R 的 recode 只映射了 elastic_net/nnls/ESMF，
#        而 CSV 里的 lasso、ridge 是小写，与 method_colors 里的大写 Lasso/Ridge
#        不匹配，factor(..., levels=names(method_colors)) 时被转成 NA 丢弃。
#   本脚本不动原脚本，独立重新生成 Figure 3，补全 recode 映射。
# ============================================================================

library(ggplot2)
library(dplyr)
library(tidyr)
library(ggsignif)

input_path <- "/media/desk16/tjn050/LSC/deconvolution result/PLOTS/ClinicalValidationWithFlowCytometry_sample_correlations.csv"
output_dir <- "/media/desk16/tjn050/LSC/LSC_deconvolution/LSC/plot/"

# 与原脚本一致的配色
method_colors <- c(
  "ESMF"         = "#9E2A2B",
  "CIBERSORTx"   = "#2D4C7A",
  "MuSiC"        = "#C17B3C",
  "DCQ"          = "#006C5B",
  "DSA"          = "#4A5D34",
  "ElasticNet"   = "#3E606F",
  "EPIC"         = "#B3570B",
  "FARDEEP"      = "#6D3D14",
  "Lasso"        = "#7A1B1D",
  "NNLS"         = "#005C73",
  "OLS"          = "#2E5C4A",
  "Ridge"        = "#556B2F",
  "RLR"          = "#4B3D33"
)

# 数据预处理：补全 recode 映射（修复点）
data_long <- read.csv(input_path) %>%
  pivot_longer(cols = -1, names_to = "Method", values_to = "Correlation") %>%
  filter(!is.na(Correlation), between(Correlation, -1, 1)) %>%
  mutate(
    Method = recode(
      Method,
      "elastic_net" = "ElasticNet",
      "nnls"        = "NNLS",
      "lasso"       = "Lasso",   # 修复：原脚本漏映射
      "ridge"       = "Ridge",   # 修复：原脚本漏映射
      "ESMF"        = "ESMF"
    ),
    Method = factor(Method, levels = names(method_colors))
  )

# 修复后诊断：列出最终保留的方法（确认 Lasso/Ridge 不再被丢）
cat("== 最终保留的方法（应含 Lasso / Ridge）==\n")
print(table(data_long$Method, useNA = "ifany"))

# 汇总统计
stats <- data_long %>%
  group_by(Method) %>%
  summarise(Mean = mean(Correlation, na.rm = TRUE),
            SD = sd(Correlation, na.rm = TRUE),
            .groups = "drop")

# 与原脚本一致的 SCI 主题
sci_theme <- theme_classic(base_size = 10) +
  theme(
    text = element_text(color = "black"),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 12, margin = margin(b = 8)),
    plot.subtitle = element_text(hjust = 0.5, size = 9, color = "gray40", margin = margin(b = 12)),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 10, face = "bold", margin = margin(r = 8)),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 9, face = "bold"),
    axis.text.y = element_text(size = 9, face = "bold"),
    axis.line = element_line(linewidth = 0.4),
    axis.ticks = element_line(linewidth = 0.3),
    legend.position = "right",
    legend.title = element_text(size = 9, face = "bold"),
    legend.text = element_text(size = 8),
    legend.key.size = unit(0.7, "lines"),
    plot.margin = margin(12, 12, 12, 12)
  )

# 绘制小提琴图（与原脚本一致）
violin_plot <- ggplot(data_long, aes(x = Method)) +
  geom_violin(aes(y = Correlation, fill = Method), alpha = 0.85, trim = FALSE,
              linewidth = 0.4, color = "gray30") +
  geom_boxplot(aes(y = Correlation), width = 0.15, fill = "white",
               outlier.shape = NA, linewidth = 0.3) +
  geom_jitter(aes(y = Correlation, color = Method), width = 0.15, size = 1.2, alpha = 0.5) +
  geom_errorbar(data = stats, aes(x = Method, y = Mean, ymin = Mean - SD, ymax = Mean + SD),
                width = 0.1, color = "black", linewidth = 0.6) +
  geom_point(data = stats, aes(x = Method, y = Mean), size = 3,
             color = "white", fill = "black", shape = 21, stroke = 0.8) +
  geom_signif(aes(y = Correlation), comparisons = list(c("ESMF", "EPIC")),
              map_signif_level = TRUE, y_position = 1.25, tip_length = 0.01,
              textsize = 3.2, vjust = -0.3) +
  scale_x_discrete(limits = names(method_colors)) +
  scale_fill_manual(values = method_colors, name = "Method") +
  scale_color_manual(values = method_colors, name = "Method", guide = "none") +
  coord_cartesian(ylim = c(-0.8, 1.3)) +
  # 注意：y 轴标题为竖排，长度受图高限制。
  #   图高 6 cm ≈ 170 pt，扣除上下 plot.margin(12+12) 后可用约 146 pt；
  #   10pt 加粗字体每字符约 5.5 pt，故标题需 ≤ 26 字符，否则上下两端会被裁切。
  #   原标题 "Pearson Correlation Coefficient (r)" 为 35 字符 → 被裁；
  #   现改为 "Pearson correlation (r)"（23 字符 ≈ 126 pt），可完整显示。
  labs(title = "ESMF Demonstrates Superior Correlation Accuracy",
       subtitle = "Comparison with EPIC Deconvolution Method",
       y = "Pearson correlation (r)") +
  sci_theme

print(violin_plot)

# 保存
timestamp <- format(Sys.time(), "%Y%m%d")
filename <- sprintf("violin_plot_fixed_%s_600dpi.png", timestamp)
full_path <- file.path(output_dir, filename)
ggsave(full_path, plot = violin_plot, device = "png",
       width = 18, height = 6, units = "cm", dpi = 600, bg = "white")
cat("\n== 已保存：", full_path, "==\n")
