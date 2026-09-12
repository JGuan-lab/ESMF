# ============================================================================
# plot_cell_type_expansion.R
#
# 细胞类型扩展前后对比图（回应 R1-Major 2）
# 老鼠下丘脑 5→6→7 种扩展，ESMF 全程第一、不坍缩。
# 对比 ESMF vs MuSiC / CIBERSORTx / FARDEEP（最强对照方法）。
# 左子图 H_Pearson，右子图 H_RMSE。
#
# 风格对齐 LSC_plots_20260531.R：敦煌配色、Arial 加粗、SCI 主题、
# 600 dpi、cairo 渲染、mm 尺寸、黄金比例。
# ============================================================================

library(ggplot2)

# ---- 输出目录 ----
output_dir <- "/media/desk16/tjn050/LSC/LSC_deconvolution/LSC/plot"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# ---- 数据（老鼠下丘脑 5/6/7 种，ESMF + 3 个最强对照）----
# 5 种=单参考下丘脑；6 种=单参考 6CT；7 种=多源 7CT
expansion_data <- data.frame(
  n_types = rep(c(5, 6, 7), times = 4),
  Method   = rep(c("ESMF", "MuSiC", "CIBERSORTx", "FARDEEP"), each = 3),
  H_Pearson = c(
    0.9502, 0.9451, 0.9202,   # ESMF
    0.9080, 0.9256, 0.8588,   # MuSiC
    0.7619, 0.8644, 0.8767,   # CIBERSORTx
    0.8970, 0.8629, 0.8118    # FARDEEP
  ),
  H_RMSE = c(
    0.0663, 0.0643, 0.0725,   # ESMF
    0.0941, 0.0772, 0.0951,   # MuSiC
    0.1558, 0.1032, 0.0876,   # CIBERSORTx
    0.1005, 0.1041, 0.1049    # FARDEEP
  )
)

# ---- 敦煌配色（对齐 LSC_plots 的 method_colors）----
method_colors <- c(
  "ESMF"       = "#9E2A2B",  # 深绯红（主方法）
  "MuSiC"      = "#C17B3C",  # 暖赭
  "CIBERSORTx" = "#2D4C7A",  # 藏青
  "FARDEEP"    = "#6D3D14"   # 巧克力棕
)
method_levels <- c("ESMF", "MuSiC", "CIBERSORTx", "FARDEEP")
expansion_data$Method <- factor(expansion_data$Method, levels = method_levels)

# ---- SCI 主题（对齐 LSC_plots 的 sci_theme）----
sci_theme <- function() {
  theme_classic(base_size = 10) +
    theme(
      text = element_text(family = "Arial", color = "black"),
      plot.title = element_text(hjust = 0.5, face = "bold", size = 12, margin = margin(b = 8)),
      axis.title = element_text(size = 10, face = "bold", color = "black"),
      axis.text = element_text(size = 9, color = "black", face = "bold"),
      axis.line = element_line(linewidth = 0.4),
      axis.ticks = element_line(linewidth = 0.3),
      panel.grid.major.y = element_line(color = "#E5D0A3", linewidth = 0.2),
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "top",
      legend.title = element_blank(),
      legend.text = element_text(size = 9, face = "bold"),
      legend.key.size = unit(0.7, "lines"),
      legend.spacing.x = unit(0.3, "cm"),
      plot.margin = margin(7, 7, 7, 7, "mm"),
      strip.text = element_text(size = 10, face = "bold", color = "black"),
      strip.background = element_rect(fill = "#F5EFDF", color = NA)
    )
}

# ---- 长表：把 Pearson 和 RMSE 合并为 metric 列，方便 facet ----
library(tidyr)
long_df <- expansion_data %>%
  pivot_longer(cols = c("H_Pearson", "H_RMSE"),
               names_to = "metric", values_to = "value")
long_df$metric <- factor(long_df$metric,
                         levels = c("H_Pearson", "H_RMSE"),
                         labels = c("Pearson correlation", "RMSE"))

# ---- 画图：分面折线 + 点 + 数值标注 ----
plot_obj <- ggplot(long_df, aes(x = n_types, y = value, color = Method, group = Method)) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2.8, fill = "white", stroke = 0.9, shape = 21) +
  geom_text(aes(label = sprintf("%.3f", value)),
            vjust = -1.3, size = 2.5, color = "black", family = "Arial",
            fontface = "bold") +
  scale_color_manual(values = method_colors) +
  scale_x_continuous(breaks = c(5, 6, 7)) +
  facet_wrap(~ metric, scales = "free_y", ncol = 2) +
  labs(x = "Number of cell types", y = NULL,
       title = "ESMF remains stable and top-ranked as cell types expand (mouse hypothalamus)") +
  sci_theme()

# ---- 保存（600 dpi, cairo, mm）----
out_png <- file.path(output_dir, "CellType_Expansion_5_6_7_20260828_600dpi.png")
ggsave(
  filename = out_png,
  plot = plot_obj,
  device = png,
  type = "cairo",
  width = 180,
  height = 80,
  units = "mm",
  dpi = 600,
  bg = "white",
  family = "Arial"
)

cat("已保存：", out_png, "\n")
cat("\n数据核对：\n")
print(expansion_data)
