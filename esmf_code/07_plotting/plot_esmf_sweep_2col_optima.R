# ============================================================================
# plot-esmf-sweep-2col-optima.R —— ESMF 遍历散点（全部结果点 + 只标红局部最优）
#
# 需求：遍历全部结果点都画出来；只把"局部最优点"标红；
#       配色 = 暖金 + 冷色，两种颜色，以 1% 为阈值（接近最优=暖金，其余=冷色），
#       不要中间过渡色。
# 数据源：esmf_sweep_paras_*.txt（全部组合）/ esmf_sweep_best_*.txt（局部最优）
# 用法：Rscript plot-esmf-sweep-2col-optima.R [PCT]   # PCT 默认 0.01
# ============================================================================
suppressMessages({library(ggplot2); library(dplyr)})

# ---- 配色：暖金 + 冷色（两色）----
col_gold <- "#E0A32E"  # 暖金（≤PCT：接近最优）
col_cold <- "#5E8C8C"  # 冷色（其余）
col_opt  <- "#CE3B2F"  # 朱砂（局部最优点）
ink      <- "#33302E"

pct <- 0.01
args <- commandArgs(trailingOnly = TRUE)
if (length(args) >= 1 && nchar(args[1]) > 0) pct <- as.numeric(args[1])
lab_near <- sprintf("≤%g%% (near optimum)", pct * 100)
lab_far  <- sprintf(">%g%%", pct * 100)

outdir <- "/media/desk16/tjn050/LSC/LSC_deconvolution/LSC/plot"
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

result_dirs <- c(
  "Mm_Brain"    = "/media/desk16/tjn050/LSC/LSC_deconvolution/LSC/results/Ablation_Study_crossPlatformConsistency_Mm_Brain_ParameterSensitivity_v2",
  "Hm_Adipose"  = "/media/desk16/tjn050/LSC/LSC_deconvolution/LSC/results/Ablation_Study_MultiSourceReferenceIntegration_Hm_Adipose_ParameterSensitivity_v2",
  "Hm_Pancreas" = "/media/desk16/tjn050/LSC/LSC_deconvolution/LSC/results/Ablation_Study_GeneralizationFromSingleReference_Hm_Pancreas_ParameterSensitivity_v2"
)
latest_file <- function(dir, pattern) {
  fs <- list.files(dir, pattern = pattern, full.names = TRUE)
  if (length(fs) == 0) return(NULL)
  fs[order(fs, decreasing = TRUE)][1]
}

all_df <- do.call(rbind, lapply(names(result_dirs), function(ds) {
  f <- latest_file(result_dirs[[ds]], "^esmf_sweep_paras_.*\\.txt$")
  d <- read.table(f, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  d$dataset <- ds; d
}))
opt_pts <- do.call(rbind, lapply(names(result_dirs), function(ds) {
  f <- latest_file(result_dirs[[ds]], "^esmf_sweep_best_.*\\.txt$")
  b <- read.table(f, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  b$dataset <- ds; b
}))

all_df$dataset <- factor(all_df$dataset, levels = names(result_dirs),
                         labels = paste0("Simulation ", seq_along(result_dirs), ": ", names(result_dirs)))
all_df$rate <- factor(all_df$rate, levels = sort(unique(as.numeric(as.character(all_df$rate)))))
opt_pts$dataset <- factor(opt_pts$dataset, levels = names(result_dirs),
                          labels = paste0("Simulation ", seq_along(result_dirs), ": ", names(result_dirs)))
opt_pts$rate <- factor(opt_pts$rate, levels = levels(all_df$rate))

# 阈值分色：每个数据集内，(best - PCC)/max(best-PCC) ≤ pct 记为"近最优"
best_pcc <- setNames(opt_pts$PCC.C, opt_pts$dataset)
all_df$dist_best <- best_pcc[all_df$dataset] - all_df$PCC.C
all_df$dist_norm <- ave(all_df$dist_best, all_df$dataset, FUN = function(x) x / max(x))
all_df$colcat <- factor(ifelse(all_df$dist_norm <= pct, lab_near, lab_far),
                        levels = c(lab_near, lab_far))

# 局部最优点的标注偏移
rate_idx <- as.integer(opt_pts$rate); n_rate <- length(levels(all_df$rate))
opt_pts$nudge_x <- ifelse(rate_idx <= n_rate / 2, 0.35, -0.35)
opt_pts$hjust   <- ifelse(rate_idx <= n_rate / 2, 0, 1)

p <- ggplot(all_df, aes(x = rate, y = PCC.C)) +
  geom_jitter(aes(color = colcat, size = PCC.C),
              width = 0.25, height = 0, shape = 16, alpha = 0.75) +
  geom_point(data = opt_pts, aes(x = rate, y = PCC.C),
             shape = 16, size = 5, color = col_opt, inherit.aes = FALSE) +
  geom_text(data = opt_pts, aes(x = rate, y = PCC.C, label = "Local optimum"),
            nudge_y = 0.02, nudge_x = opt_pts$nudge_x, hjust = opt_pts$hjust,
            size = 2.2, fontface = "bold", color = col_opt, family = "Arial") +
  scale_color_manual(name = "distance to\noptimum\n(ΔPearson)",
                     values = setNames(c(col_gold, col_cold), c(lab_near, lab_far))) +
  scale_size_continuous(range = c(1, 4.5), name = "Pearson") +
  facet_wrap(~ dataset, nrow = 1) +
  labs(x = "rate (overall prior strength)",
       y = "Pearson",
       title = "ESMF hyperparameter search (para1 = 1; grid over rate × para2 × para3)") +
  guides(color = guide_legend(order = 1), size = guide_legend(order = 2)) +
  theme_classic(base_family = "Arial") +
  theme(
    text            = element_text(family = "Arial", color = ink),
    plot.title      = element_text(size = 11, face = "bold", hjust = 0.5),
    axis.text       = element_text(size = 9,  face = "bold"),
    axis.title      = element_text(size = 10, face = "bold"),
    strip.text      = element_text(size = 9,  face = "bold"),
    strip.background = element_blank(),
    legend.title    = element_text(size = 7, face = "bold"),
    legend.text     = element_text(size = 6),
    legend.key.size = grid::unit(0.35, "cm"),
    panel.grid.major.y = element_line(color = "#E8E8E8", linewidth = 0.3, linetype = "dashed"),
    plot.margin     = margin(5, 6, 5, 5, "mm")
  )

out <- file.path(outdir, sprintf("TUNE_ESMF_SWEEP_2col_%gpct.png", pct * 100))
ggsave(out, p, device = png, type = "cairo",
       width = 200 / 25.4, height = 105 / 25.4, units = "in", dpi = 600, bg = "white")
cat("saved:", out, "\n")
