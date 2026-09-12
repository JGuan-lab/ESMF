# ============================================================================
# plot-per-sample-violin-pearson-rmse.R
# 每个数据集一张图，同时展示逐样本 Pearson 与 RMSE（上下两个面板，共享 x 轴）。
# 配色 / 方法顺序 / SCI 主题 均参考 LSC_plots_20260531.R 的小提琴图。
#
# 逐样本指标定义（与主分析口径一致）：
#   对每个 (method, tissue/样本)，在该样本内的所有 CT 上计算
#     Pearson = cor(observed_values, expected_values)
#     RMSE    = sqrt(mean((observed_values - expected_values)^2))
#   于是每个方法得到一个"每个样本一个值"的分布，用小提琴图展示。
#
# 输出：每个数据集一张 PNG（10 x 8 in, 300 dpi）
# ============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(rlang)
})

# patchwork 用于上下拼图；若未安装则给出明确提示
if (!requireNamespace("patchwork", quietly = TRUE)) {
  stop("需要 patchwork 包来拼接上下两个面板。请先运行：install.packages('patchwork')")
}
library(patchwork)

# ----------------------------------------------------------------------------
# 1. 配置：三个数据集
# ----------------------------------------------------------------------------
RESULTS_ROOT <- "/media/desk16/tjn050/LSC/LSC_deconvolution/LSC/results"
OUT_DIR      <- "/media/desk16/tjn050/LSC/LSC_deconvolution/LSC/plot"

if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE)

# 数据集命名前缀说明：Hm = Human（人类），Mm = Mouse（小鼠，Mus musculus）
DATASETS <- list(
  list(
    name = "Mm_Brain_crossPlatformConsistency",
    sim  = "Simulation 1",
    label = "Mm_Brain (cross-platform consistency)",
    path = file.path(RESULTS_ROOT,
                     "crossPlatformConsistency_Mm_Brain",
                     "per_sample_results", "per_sample_long.csv")
  ),
  list(
    name = "Hm_Adipose_MultiSourceReferenceIntegration",
    sim  = "Simulation 2",
    label = "Hm_Adipose (multi-source reference integration)",
    path = file.path(RESULTS_ROOT,
                     "MultiSourceReferenceIntegration_Hm_Adipose",
                     "per_sample_results", "per_sample_long.csv")
  ),
  list(
    name = "Hm_Pancreas_GeneralizationFromSingleReference",
    sim  = "Simulation 3",
    label = "Hm_Pancreas (generalization from single reference)",
    path = file.path(RESULTS_ROOT,
                     "GeneralizationFromSingleReference_Hm_Pancreas",
                     "per_sample_results", "per_sample_long.csv")
  )
)

# ----------------------------------------------------------------------------
# 2. 配色与方法顺序（严格照搬 LSC_plots_20260531.R 的 method_colors）
#    ESMF 第一（深绛红），其余按参考脚本顺序
# ----------------------------------------------------------------------------
method_colors <- c(
  "ESMF"       = "#9E2A2B",  # Primary: Deep crimson
  "CIBERSORTx" = "#2D4C7A",  # Rich navy blue
  "MuSiC"      = "#C17B3C",  # Warm terracotta
  "DCQ"        = "#006C5B",  # Deep teal
  "DSA"        = "#4A5D34",  # Forest green
  "ElasticNet" = "#3E606F",  # Slate blue
  "EPIC"       = "#B3570B",  # Burnt orange
  "FARDEEP"    = "#6D3D14",  # Chocolate brown
  "Lasso"      = "#7A1B1D",  # Maroon red
  "NNLS"       = "#005C73",  # Deep cyan
  "OLS"        = "#2E5C4A",  # Jade green
  "Ridge"      = "#556B2F",  # Olive green
  "RLR"        = "#4B3D33"   # Dark taupe
)

# CSV 内部方法名 -> 参考脚本显示名
method_recode <- c(
  "elastic_net" = "ElasticNet",
  "nnls"        = "NNLS",
  "lasso"       = "Lasso",
  "ridge"       = "Ridge"
)

# ----------------------------------------------------------------------------
# 3. SCI 主题（照搬 LSC_plots_20260531.R 的 sci_theme）
# ----------------------------------------------------------------------------
sci_theme <- theme_classic(base_family = "Arial") +
  theme(
    text = element_text(family = "Arial", color = "black"),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 11,
                              margin = margin(t = 10, b = 6)),
    plot.subtitle = element_text(hjust = 0.5, size = 9, color = "gray40",
                                 margin = margin(b = 12)),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 10, face = "bold", color = "black"),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1,
                               size = 9, face = "bold", color = "black"),
    axis.text.y = element_text(size = 9, face = "bold", color = "black"),
    axis.line = element_line(linewidth = 0.4),
    axis.ticks = element_line(linewidth = 0.3),
    legend.position = "none",
    panel.grid.major.y = element_line(color = "#E5D0A3", linewidth = 0.2),
    plot.margin = margin(8, 8, 8, 8, "mm")
  )

# ----------------------------------------------------------------------------
# 4. 工具函数：读长表 -> 计算逐样本指标
# ----------------------------------------------------------------------------
# 安全相关系数：任一方方差为 0 时返回 NA，且不产生 warning
safe_cor <- function(x, y) {
  if (length(x) < 2 || length(y) < 2) return(NA_real_)
  sx <- sd(x, na.rm = TRUE); sy <- sd(y, na.rm = TRUE)
  if (is.na(sx) || is.na(sy) || sx < 1e-12 || sy < 1e-12) return(NA_real_)
  suppressWarnings(cor(x, y, use = "complete.obs"))
}

compute_per_sample_metrics <- function(csv_path) {
  df <- read.csv(csv_path, stringsAsFactors = FALSE)

  # 统一方法名到参考脚本的显示名
  df$method_disp <- ifelse(df$method %in% names(method_recode),
                           unname(method_recode[df$method]),
                           df$method)
  # 颜色表里没有的方法（应为空），warning 提示
  unknown <- setdiff(unique(df$method_disp), names(method_colors))
  if (length(unknown) > 0) {
    warning("以下方法名不在 method_colors 中，将被过滤: ",
            paste(unknown, collapse = ", "))
  }
  missing_m <- setdiff(names(method_colors), unique(df$method_disp))
  if (length(missing_m) > 0) {
    warning("数据集中缺少方法: ", paste(missing_m, collapse = ", "))
  }

  df <- df[df$method_disp %in% names(method_colors), ]
  df$Method <- factor(df$method_disp, levels = names(method_colors))

  # ---- 逐样本（tissue 内跨 CT）指标 + 退化诊断 ----
  metrics_raw <- df %>%
    group_by(Method, tissue) %>%
    summarise(
      exp_sd      = sd(expected_values, na.rm = TRUE),
      obs_sd      = sd(observed_values, na.rm = TRUE),
      pearson_raw = safe_cor(observed_values, expected_values),
      rmse        = sqrt(mean((observed_values - expected_values)^2, na.rm = TRUE)),
      .groups = "drop"
    ) %>%
    mutate(
      # (a) 样本本身无区分度：真实比例在所有 CT 上完全相同，与方法无关
      sample_invalid    = is.na(exp_sd) | exp_sd < 1e-12,
      # (b) 方法退化：真实比例有区分度，但该方法把所有 CT 估成同一值
      method_degenerate = (!sample_invalid) & (is.na(obs_sd) | obs_sd < 1e-12)
    )

  n_sample_invalid <- length(unique(metrics_raw$tissue[metrics_raw$sample_invalid]))

  # (a) 类样本在所有方法中统一排除，保证各方法比较同一批样本、
  #     且 Pearson 与 RMSE 两个指标口径一致
  metrics <- metrics_raw %>%
    filter(!sample_invalid) %>%
    mutate(
      # (b) 类：方法退化 -> pearson 记为 0（完全无预测能力），不再静默丢弃
      pearson = ifelse(method_degenerate, 0, pearson_raw)
    ) %>%
    filter(!is.na(pearson)) %>%
    select(Method, tissue, pearson, rmse, method_degenerate)

  # ---- 诊断输出 ----
  cat("   >> 样本本身无区分度(exp_sd=0)、在所有方法中统一排除的样本数: ",
      n_sample_invalid, "\n", sep = "")
  degen_tab <- metrics_raw %>%
    filter(!sample_invalid) %>%
    group_by(Method) %>%
    summarise(
      n_samples      = n(),
      n_degenerate   = sum(method_degenerate),
      pct_degenerate = round(100 * mean(method_degenerate), 2),
      .groups = "drop"
    ) %>%
    arrange(desc(n_degenerate))
  cat("   >> 各方法退化样本数（把所有 CT 估成同一值，pearson 记为 0）:\n")
  print(degen_tab)

  metrics
}

# 单个指标面板的小提琴图
make_violin_panel <- function(metrics, metric_col, ylab, show_x_text) {
  mcol <- rlang::sym(metric_col)

  # 汇总统计（均值 ± SD），用于误差棒 + 均值点
  stats <- metrics %>%
    group_by(Method) %>%
    summarise(
      Mean = mean(!!mcol, na.rm = TRUE),
      SD   = sd(!!mcol, na.rm = TRUE),
      .groups = "drop"
    )

  p <- ggplot(metrics, aes(x = Method)) +
    geom_violin(aes(y = !!mcol, fill = Method),
                alpha = 0.85, trim = FALSE, linewidth = 0.4, color = "gray30") +
    geom_boxplot(aes(y = !!mcol), width = 0.15, fill = "white",
                 outlier.shape = NA, linewidth = 0.3) +
    geom_jitter(aes(y = !!mcol, color = Method),
                width = 0.15, size = 0.9, alpha = 0.35) +
    geom_errorbar(data = stats,
                  aes(x = Method, y = Mean, ymin = Mean - SD, ymax = Mean + SD),
                  width = 0.1, color = "black", linewidth = 0.6) +
    geom_point(data = stats, aes(x = Method, y = Mean),
               size = 3, color = "white", fill = "black",
               shape = 21, stroke = 0.8) +
    scale_x_discrete(limits = names(method_colors)) +
    scale_fill_manual(values = method_colors, name = "Method") +
    scale_color_manual(values = method_colors, name = "Method", guide = "none") +
    labs(y = ylab) +
    sci_theme

  # 上面板隐藏 x 轴文字（共享 x 轴），下面板显示
  if (!show_x_text) {
    p <- p + theme(axis.text.x = element_blank(),
                   axis.ticks.x = element_blank())
  }

  p
}

# ----------------------------------------------------------------------------
# 5. 主循环：每个数据集一张图
# ----------------------------------------------------------------------------
for (ds in DATASETS) {
  cat("\n=========================================================\n")
  cat(">> 数据集:", ds$name, "\n")
  cat(">> 读文件:", ds$path, "\n")

  if (!file.exists(ds$path)) {
    warning("文件不存在，跳过: ", ds$path)
    next
  }

  metrics <- compute_per_sample_metrics(ds$path)

  n_methods <- length(unique(metrics$Method))
  n_samples <- length(unique(metrics$tissue))
  cat(">> 方法数:", n_methods, " 样本数:", n_samples,
      " 总指标点:", nrow(metrics), "\n")

  # 打印每个方法的均值（便于核对）
  cat("\n-- 逐方法均值 --\n")
  print(metrics %>%
          group_by(Method) %>%
          summarise(
            mean_pearson = round(mean(pearson, na.rm = TRUE), 4),
            mean_rmse    = round(mean(rmse, na.rm = TRUE), 4),
            n_samples    = n(),
            .groups = "drop"
          ) %>%
          arrange(desc(mean_pearson)))

  # 上面板：Pearson（隐藏 x 文字）；下面板：RMSE（显示 x 文字）
  # 标题/副标题放在上面板自身的 ggtitle，这样标题上方才能真正留白
  p_pearson <- make_violin_panel(metrics, "pearson",
                                 "Pearson correlation (r)",
                                 show_x_text = FALSE) +
    ggtitle(paste(ds$sim, ds$label, sep = ": "),
            subtitle = "Per-sample accuracy (each point = one sample)") +
    theme(plot.margin = margin(8, 8, 2, 8, "mm"))
  p_rmse    <- make_violin_panel(metrics, "rmse",
                                 "RMSE",
                                 show_x_text = TRUE) +
    theme(plot.margin = margin(2, 8, 8, 8, "mm"))

  combined <- (p_pearson / p_rmse)

  out_png <- file.path(OUT_DIR, paste0("PerSample_Pearson_RMSE_", ds$name, ".png"))
  ggsave(out_png, combined, device = png, type = "cairo",
         width = 180 / 25.4, height = 140 / 25.4, units = "in",
         dpi = 600, bg = "white")
  cat(">> 已保存:", out_png, "\n")
}

cat("\n全部完成。输出目录:", OUT_DIR, "\n")
