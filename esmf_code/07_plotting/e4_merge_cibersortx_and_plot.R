# ============================================================================
# e4_merge_cibersortx_and_plot.R
#
# 用途：
#   1. 将补算的 CIBERSORTx 结果（E4_CIBERSORTx_backfill_*.txt）合并进 E4 v3
#      汇总结果（E4_expression_shift_v3_*.txt），生成一张统一汇总表。
#   2. 生成对比图：场景 A / B 分面，H_RMSE 与 H_Pearson 两个面板，
#      全部 12 个方法（10 bulk + MuSiC + CIBERSORTx，另有 ESMF），
#      CIBERSORTx 用高亮颜色 + 加粗标出。
# ============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(scales)
})

RESULTS_DIR <- "/media/desk16/tjn050/LSC/LSC_deconvolution/LSC/results/Ablation_Study_crossPlatformConsistency_Mm_Brain_ParameterSensitivity_v2"

# ---- 1. 定位最新文件 ----
v3_file <- sort(list.files(RESULTS_DIR, pattern = "^E4_expression_shift_v3_[0-9]{8}_[0-9]{6}\\.txt$", full.names = TRUE),
                decreasing = TRUE)[1]
backfill_file <- sort(list.files(RESULTS_DIR, pattern = "^E4_CIBERSORTx_backfill_.*\\.txt$", full.names = TRUE),
                      decreasing = TRUE)[1]
if (is.na(v3_file) || is.na(backfill_file)) stop("未找到 E4 v3 或 CIBERSORTx backfill 结果文件")

cat("v3 结果文件：", v3_file, "\n")
cat("CIBERSORTx backfill 文件：", backfill_file, "\n\n")

# ---- 2. 读取并统一列格式 ----
v3 <- read.table(v3_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
cbx <- read.table(backfill_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# 统一列集合：v3 有 iter/seed/W_RMSE/W_Pearson，cbx 没有（补 NA）
common_cols <- c("label", "scenario", "method", "shift", "H_RMSE", "H_Pearson",
                 "iter", "seed", "W_RMSE", "W_Pearson")
for (cc in setdiff(common_cols, names(cbx))) cbx[[cc]] <- NA
cbx <- cbx[, common_cols]
for (cc in setdiff(common_cols, names(v3))) v3[[cc]] <- NA
v3 <- v3[, common_cols]

# shift 统一为数值
v3$shift <- as.numeric(v3$shift)
cbx$shift <- as.numeric(cbx$shift)

merged <- rbind(v3, cbx)
merged <- merged[order(merged$scenario, merged$method, merged$shift), ]

# ---- 3. 写出合并后的统一汇总表 ----
current_datetime <- format(Sys.time(), "%Y%m%d_%H%M%S")
merged_file <- file.path(RESULTS_DIR, paste0("E4_expression_shift_v3_merged_", current_datetime, ".txt"))
write.table(merged, file = merged_file, sep = "\t", quote = FALSE,
            col.names = TRUE, row.names = FALSE)
cat("合并结果已保存到：", merged_file, "\n\n")

# ---- 4. 绘图 ----
suppressPackageStartupMessages({
  library(ggrepel)
  library(patchwork)
})

# 高对比 13 色调色板（色盲友好），并按方法名固定映射
method_order <- sort(unique(merged$method))
pal13 <- c(
  ESMF        = "#E63946",  # 亮红——主角，最醒目
  CIBERSORTx  = "#1D3557",  # 深蓝
  MuSiC       = "#2A9D8F",  # 青绿
  lasso       = "#F4A261",  # 橙
  OLS         = "#6A4C93",  # 紫
  RLR         = "#00B4D8",  # 天蓝
  nnls        = "#EF476F",  # 粉红
  ridge       = "#06D6A0",  # 薄荷绿
  FARDEEP     = "#8338EC",  # 亮紫
  elastic_net = "#FB8500",  # 深橙
  EPIC        = "#9B5DE5",  # 淡紫
  DCQ         = "#80B918",  # 草绿
  DSA         = "#BB3E03"   # 棕
)
pal <- pal13[method_order]

# 方法名缩写映射（折线图末端标注用，避免长名重叠）
method_abbr <- c(
  ESMF = "ESMF", CIBERSORTx = "CIBERSORTx", MuSiC = "MuSiC",
  lasso = "lasso", OLS = "OLS", RLR = "RLR", nnls = "nnls",
  ridge = "ridge", FARDEEP = "FARDEEP", elastic_net = "elastic_net",
  EPIC = "EPIC", DCQ = "DCQ", DSA = "DSA")

merged$method_f <- factor(merged$method, levels = method_order)
merged$scenario_f <- factor(merged$scenario, levels = c("A", "B"),
                            labels = c("Scenario A: Target noise",
                                       "Scenario B: Reference shift"))

# ==========================================================================
# 图 1：基线（shift=0）优势柱状图 —— 一眼看出 ESMF 断崖领先
# ==========================================================================
base <- merged[merged$shift == 0, ]
# 柱状图用两个 scenario 的基线（两者 shift=0 相同），取 scenario A 代表即可
baseA <- base[base$scenario == "A", ]
# 按 H_Pearson 降序排序
baseA$method_f2 <- factor(baseA$method, levels = baseA$method[order(baseA$H_Pearson, decreasing = TRUE)])

bar_rmse <- ggplot(baseA, aes(x = method_f2, y = H_RMSE, fill = method)) +
  geom_col(width = 0.7) +
  geom_text(aes(label = sprintf("%.3f", H_RMSE)), vjust = -0.4, size = 3) +
  scale_fill_manual(values = pal, guide = "none") +
  labs(x = NULL, y = "H RMSE (lower = better)", title = "Baseline (shift=0) H RMSE") +
  theme_bw(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
        panel.grid.minor = element_blank())

bar_pear <- ggplot(baseA, aes(x = method_f2, y = H_Pearson, fill = method)) +
  geom_col(width = 0.7) +
  geom_text(aes(label = sprintf("%.3f", H_Pearson)), vjust = -0.4, size = 3) +
  scale_fill_manual(values = pal, guide = "none") +
  labs(x = NULL, y = "H Pearson (higher = better)", title = "Baseline (shift=0) H Pearson") +
  theme_bw(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
        panel.grid.minor = element_blank())

# ==========================================================================
# 图 2：鲁棒性趋势折线图 —— 线端直接标注方法名（不靠 legend）
# ==========================================================================
plot_trend <- function(metric) {
  d <- merged
  d$value <- d[[metric]]
  ylab <- if (metric == "H_RMSE") "H RMSE (lower = better)" else "H Pearson (higher = better)"
  ttl <- if (metric == "H_RMSE") "H RMSE vs shift" else "H Pearson vs shift"

  # 线端标注点：每个方法取 shift 最大值（场景 B 的右端点）作为标注锚点
  lab_pts <- d[d$shift == max(d$shift) & d$scenario == "B", ]

  ggplot(d, aes(x = shift, y = value, group = method, color = method)) +
    geom_line(aes(linewidth = method), alpha = 0.9) +
    geom_point(size = 1.6) +
    # 线端直接标注方法名
    geom_text_repel(data = lab_pts,
                    aes(x = shift, y = value, label = method_abbr[method]),
                    nudge_x = 0.18, direction = "y", hjust = 0,
                    size = 3.2, fontface = "bold", segment.color = "grey60",
                    segment.size = 0.3, min.segment.length = 0,
                    max.overlaps = 20) +
    facet_grid(~ scenario_f) +
    scale_color_manual(values = pal, guide = "none") +
    scale_linewidth_manual(values = setNames(ifelse(method_order == "ESMF", 1.4, 0.8),
                                             method_order), guide = "none") +
    scale_x_continuous(breaks = c(0, 0.1, 0.2, 0.5, 1), expand = expansion(mult = c(0.05, 0.3))) +
    labs(x = "Expression shift level", y = ylab, title = ttl) +
    theme_bw(base_size = 12) +
    theme(panel.grid.minor = element_blank(),
          strip.text = element_text(face = "bold", size = 10))
}

p_trend_rmse <- plot_trend("H_RMSE")
p_trend_pear <- plot_trend("H_Pearson")

# ---- 输出 ----
fig_bar_rmse <- file.path(RESULTS_DIR, paste0("E4_bar_H_RMSE_", current_datetime))
ggsave(paste0(fig_bar_rmse, ".pdf"), bar_rmse, width = 8, height = 5, device = "pdf")
ggsave(paste0(fig_bar_rmse, ".png"), bar_rmse, width = 8, height = 5, dpi = 300, device = "png")

fig_bar_pear <- file.path(RESULTS_DIR, paste0("E4_bar_H_Pearson_", current_datetime))
ggsave(paste0(fig_bar_pear, ".pdf"), bar_pear, width = 8, height = 5, device = "pdf")
ggsave(paste0(fig_bar_pear, ".png"), bar_pear, width = 8, height = 5, dpi = 300, device = "png")

fig_trend_rmse <- file.path(RESULTS_DIR, paste0("E4_trend_H_RMSE_", current_datetime))
ggsave(paste0(fig_trend_rmse, ".pdf"), p_trend_rmse, width = 12, height = 5, device = "pdf")
ggsave(paste0(fig_trend_rmse, ".png"), p_trend_rmse, width = 12, height = 5, dpi = 300, device = "png")

fig_trend_pear <- file.path(RESULTS_DIR, paste0("E4_trend_H_Pearson_", current_datetime))
ggsave(paste0(fig_trend_pear, ".pdf"), p_trend_pear, width = 12, height = 5, device = "pdf")
ggsave(paste0(fig_trend_pear, ".png"), p_trend_pear, width = 12, height = 5, dpi = 300, device = "png")

# 组合大图：柱状图（上）+ 折线趋势（下）
p_comb <- (bar_rmse | bar_pear) / (p_trend_rmse / p_trend_pear) +
  plot_annotation(title = "E4 Expression Shift: Method Comparison",
                  theme = theme(plot.title = element_text(face = "bold", size = 15)))
fig_comb <- file.path(RESULTS_DIR, paste0("E4_merged_all_", current_datetime))
ggsave(paste0(fig_comb, ".pdf"), p_comb, width = 14, height = 12, device = "pdf")
ggsave(paste0(fig_comb, ".png"), p_comb, width = 14, height = 12, dpi = 300, device = "png")

cat("对比图已保存：\n")
cat("  柱状图 H_RMSE  :", fig_bar_rmse, "\n")
cat("  柱状图 H_Pearson:", fig_bar_pear, "\n")
cat("  折线趋势 H_RMSE :", fig_trend_rmse, "\n")
cat("  折线趋势 H_Pearson:", fig_trend_pear, "\n")
cat("  组合大图        :", fig_comb, "\n")

cat("\n========== 合并后 CIBERSORTx 行 ==========\n")
print(merged[merged$method == "CIBERSORTx",
             c("scenario", "method", "shift", "H_RMSE", "H_Pearson")])
