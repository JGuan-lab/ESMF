# ============================================================================
# plot-ablation.R —— 消融实验可视化（修复 fit_signature mu/size bug 后）
#
# 图1 A1 正则化消融折线（Pearson + RMSE facet）
# 图2 A2 marker 面板消融分组柱状（Pearson + RMSE facet）
# 只画 3 个"好"数据集（胰腺/脂肪/大脑），骨髓忽略。风格对齐 LSC_plots_20260531.R
# ============================================================================

library(ggplot2)
library(dplyr)
library(tidyr)

output_dir <- "/media/desk16/tjn050/LSC/LSC_deconvolution/LSC/plot"

# 敦煌配色（按数据集，只含 3 个好数据集）
dataset_colors <- c(
  "Hm_Pancreas" = "#8A2D20",  # 砖红
  "Hm_Adipose"  = "#4F695F",  # 青灰绿
  "Mm_Brain"    = "#3A4660"   # 深绀青
)

# 修复后数据（只含 3 个好数据集，骨髓忽略不提）
ablation_data <- data.frame(
  dataset = rep(c("Hm_Pancreas", "Hm_Adipose", "Mm_Brain"), each = 8),
  variant = rep(c("rate0", "rate0.1", "rate1", "rate10",
                  "marker_base", "marker_lit", "marker_DE50", "marker_random"), 3),
  H_Pearson = c(
    0.9254, 0.9378, 0.9429, 0.9366,  0.9429, 0.9429, 0.9268, 0.3649,   # 胰腺
    0.9588, 0.9594, 0.9629, 0.9655,  0.9629, 0.9629, 0.9690, 0.9255,   # 脂肪
    0.9250, 0.9263, 0.9279, 0.9292,  0.9279, 0.9279, 0.8976, 0.8608    # 大脑
  ),
  H_RMSE = c(
    0.101, 0.0919, 0.088, 0.0929,    0.088, 0.088, 0.1003, 0.3109,
    0.0769, 0.0764, 0.0731, 0.0703,  0.0731, 0.0731, 0.0677, 0.1015,
    0.0972, 0.0963, 0.0951, 0.0942,  0.0951, 0.0951, 0.1126, 0.1301
  )
)

ablation_data$dataset <- factor(ablation_data$dataset,
                                levels = c("Hm_Pancreas", "Hm_Adipose", "Mm_Brain"))

# 统一 SCI 主题（对齐 0531）
sci_theme <- theme_classic(base_family = "Arial") +
  theme(
    text = element_text(family = "Arial", color = "black"),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 12, margin = margin(b = 8)),
    axis.title = element_text(size = 10, face = "bold", color = "black"),
    axis.text = element_text(size = 9, face = "bold", color = "black"),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    axis.line = element_line(linewidth = 0.4),
    axis.ticks = element_line(linewidth = 0.3),
    legend.position = "top",
    legend.title = element_text(size = 9, face = "bold"),
    legend.text = element_text(size = 8),
    strip.text = element_text(size = 10, face = "bold"),
    strip.background = element_blank(),
    panel.grid.major.y = element_line(color = "#E5D0A3", linewidth = 0.2),
    plot.margin = margin(8, 8, 8, 8, "mm")
  )

# ============================================================
# 图1 A1 正则化消融折线
# ============================================================
a1_data <- ablation_data %>%
  filter(grepl("rate", variant)) %>%
  mutate(rate = factor(variant, levels = c("rate0","rate0.1","rate1","rate10"),
                       labels = c("0 (pure NMF)","0.1","1 (default)","10")))

a1_long <- a1_data %>%
  pivot_longer(cols = c(H_Pearson, H_RMSE),
               names_to = "metric", values_to = "value") %>%
  mutate(metric = factor(metric, levels = c("H_Pearson","H_RMSE"),
                         labels = c("Pearson correlation","RMSE")))

p1 <- ggplot(a1_long, aes(x = rate, y = value, color = dataset, group = dataset)) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2.5, shape = 21, fill = "white", stroke = 1) +
  facet_wrap(~ metric, ncol = 2, scales = "free_y") +
  scale_color_manual(values = dataset_colors, name = "Dataset") +
  labs(x = "Overall prior strength (rate)", y = NULL,
       title = "A1: Regularization ablation — prior strength effect") +
  sci_theme

ggsave(file.path(output_dir, "ABLATION_A1_regularization.png"), p1,
       device = png, type = "cairo", width = 180/25.4, height = 90/25.4,
       units = "in", dpi = 600, bg = "white")

# ============================================================
# 图2 A2 marker 面板消融分组柱状
# ============================================================
a2_data <- ablation_data %>%
  filter(grepl("marker", variant)) %>%
  mutate(marker = factor(variant,
                         levels = c("marker_base","marker_lit","marker_DE50","marker_random"),
                         labels = c("Baseline\nmarkers","Literature\naugmented","DE top50","Random\ngenes")))

a2_long <- a2_data %>%
  pivot_longer(cols = c(H_Pearson, H_RMSE),
               names_to = "metric", values_to = "value") %>%
  mutate(metric = factor(metric, levels = c("H_Pearson","H_RMSE"),
                         labels = c("Pearson correlation","RMSE")))

p2 <- ggplot(a2_long, aes(x = marker, y = value, fill = dataset)) +
  geom_col(position = position_dodge(0.8), width = 0.7,
           color = "#5D4C46", linewidth = 0.3) +
  facet_grid(metric ~ dataset, scales = "free_y") +
  scale_fill_manual(values = dataset_colors, name = "Dataset") +
  labs(x = "Marker panel", y = NULL,
       title = "A2: Marker panel ablation — marker selection effect") +
  sci_theme +
  theme(strip.text.x = element_blank())

ggsave(file.path(output_dir, "ABLATION_A2_marker.png"), p2,
       device = png, type = "cairo", width = 180/25.4, height = 90/25.4,
       units = "in", dpi = 600, bg = "white")

# ============================================================
# 图2b A2 marker 贡献（Δ 相对 baseline 的下降）
# ============================================================
a2_delta <- a2_data %>%
  group_by(dataset) %>%
  mutate(benchmark = H_Pearson[variant == "marker_base"],
         delta = H_Pearson - benchmark) %>%
  ungroup() %>%
  filter(variant != "marker_base") %>%
  mutate(marker = factor(variant,
                         levels = c("marker_lit","marker_DE50","marker_random"),
                         labels = c("Literature\naugmented","DE top50","Random\ngenes")))

p2b <- ggplot(a2_delta, aes(x = marker, y = delta, fill = dataset)) +
  geom_col(position = position_dodge(0.8), width = 0.7,
           color = "#5D4C46", linewidth = 0.3) +
  # 数值标签：把"零增量"也明确标出来（literature 显示 0.0000）
  geom_text(aes(label = sprintf("%.4f", delta)),
            position = position_dodge(0.8),
            vjust = ifelse(a2_delta$delta < 0, 1.2, -0.3),
            size = 2.4, fontface = "bold", family = "Arial",
            color = "black") +
  geom_hline(yintercept = 0, color = "black", linewidth = 0.4) +
  facet_grid(. ~ dataset, scales = "free_y") +
  scale_fill_manual(values = dataset_colors, name = "Dataset") +
  scale_y_continuous(expand = expansion(mult = c(0.08, 0.15))) +
  labs(x = "Marker panel (relative to baseline)", y = "Δ H_Pearson (vs baseline)",
       title = "A2: Marker ablation — performance loss relative to baseline") +
  sci_theme +
  theme(strip.text.x = element_blank())

ggsave(file.path(output_dir, "ABLATION_A2_marker_delta.png"), p2b,
       device = png, type = "cairo", width = 180/25.4, height = 90/25.4,
       units = "in", dpi = 600, bg = "white")

cat("消融图已保存到", output_dir, "\n")
cat("  ABLATION_A1_regularization.png\n")
cat("  ABLATION_A2_marker.png\n")
cat("  ABLATION_A2_marker_delta.png\n")
