# ============================================================================
# plot-ablation-contribution.R —— 消融"双组件贡献"图（回应 R1-M5 审稿意见）
#
# 审稿人要求分别评估两个组件贡献：
#   (1) marker 筛选（cell type specific marker screening）
#   (2) 弹性签名正则化（elastic signature regularization）
#
# 画法（浮动增量柱 + 灰色基准线，基准与增量同时可见）：
#   - 灰色短线 = 零状态基准（rate=0 纯 NMF / random 基因）
#   - 彩色柱   = 相对基准的增量（柱底 = 基准线，柱长 = 增量）
#   - y 轴起点拉近基准值（截断），微小增量也能看清
# 图 A：正则化贡献 —— 基准 rate=0，增量档 rate 0.1/1/10
# 图 B：marker 贡献 —— 基准 random，增量档 DE top50 / ESMF-selected
# ============================================================================

library(ggplot2)
library(dplyr)
library(tidyr)

output_dir <- "/media/desk16/tjn050/LSC/LSC_deconvolution/LSC/plot"

base_grey <- "#9E9E9E"
dataset_colors <- c(
  "Hm_Pancreas" = "#8A2D20",  # 砖红
  "Hm_Adipose"  = "#4F695F",  # 青灰绿
  "Mm_Brain"    = "#3A4660"   # 深绀青
)

sci_theme <- theme_classic(base_family = "Arial") +
  theme(
    text = element_text(family = "Arial", color = "black"),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 11, margin = margin(b = 4)),
    plot.subtitle = element_text(hjust = 0.5, face = "plain", size = 9, color = "black", margin = margin(b = 6)),
    axis.title = element_text(size = 10, face = "bold", color = "black"),
    axis.text = element_text(size = 9, face = "bold", color = "black"),
    axis.line = element_line(linewidth = 0.4),
    axis.ticks = element_line(linewidth = 0.3),
    legend.position = "top",
    legend.title = element_text(size = 9, face = "bold"),
    legend.text = element_text(size = 8),
    strip.text = element_text(size = 10, face = "bold"),
    strip.background = element_blank(),
    panel.grid.major.y = element_line(color = "#E5D0A3", linewidth = 0.2),
    plot.margin = margin(12, 8, 8, 8, "mm")
  )

w <- 0.24
offset_map <- c("Hm_Pancreas" = -0.28, "Hm_Adipose" = 0, "Mm_Brain" = 0.28)

# ============================================================
# 数据读取：全部从 v3 消融结果文件读，不写死数值
#   每个数据集取其 v3 目录下「最新」一份 ablation_results_*.txt；
#   该文件同时含 A1（rate 消融）与 A2（marker 消融）。
#   数据集命名前缀：Hm = Human，Mm = Mouse（Mus musculus）
# ============================================================
RESULTS_ROOT <- "/media/desk16/tjn050/LSC/LSC_deconvolution/LSC/results"

DS_DIR <- c(
  "Hm_Pancreas" = "Ablation_Study_GeneralizationFromSingleReference_Hm_Pancreas/v3",
  "Hm_Adipose"  = "Ablation_Study_MultiSourceReferenceIntegration_Hm_Adipose/v3",
  "Mm_Brain"    = "Ablation_Study_crossPlatformConsistency_Mm_Brain/v3"
)

read_latest_ablation <- function(ds) {
  dp <- file.path(RESULTS_ROOT, DS_DIR[[ds]])
  fl <- list.files(dp, pattern = "^ablation_results_.*\\.txt$", full.names = TRUE)
  if (length(fl) == 0) stop("未找到 ablation_results 文件：", dp)
  fl <- fl[which.max(file.mtime(fl))]          # 取最新一次运行
  cat(sprintf("[%s] 读取: %s\n", ds, basename(fl)))
  abl_src[[ds]] <<- fl                          # 记录来源文件，供溯源清单使用
  read.delim(fl, stringsAsFactors = FALSE, check.names = FALSE)
}

abl <- list()
abl_src <- list()                                # 每个数据集实际读取的结果文件路径
for (ds in names(DS_DIR)) abl[[ds]] <- read_latest_ablation(ds)

get_val <- function(ds, label) {
  d <- abl[[ds]]
  v <- d$H_Pearson[d$label == label]
  if (length(v) != 1) stop(sprintf("[%s] 未唯一匹配标签：%s", ds, label))
  as.numeric(v)
}

# ============================================================
# 图 A：正则化贡献（灰色基准 = rate=0 纯 NMF）
# ============================================================
rate_levels <- c("0 (pure NMF)", "0.1", "1", "10")
rate_label_map <- c("0 (pure NMF)" = "A1_pure_NMF_rate0",
                    "0.1"          = "A1_weak_rate0.1",
                    "1"            = "A1_benchmark_rate1",
                    "10"           = "A1_strong_rate10")

reg_data <- expand.grid(rate = rate_levels, dataset = names(DS_DIR),
                        stringsAsFactors = FALSE)
reg_data$value <- mapply(function(rl, ds) get_val(ds, rate_label_map[[rl]]),
                         reg_data$rate, reg_data$dataset)
reg_data$base <- sapply(reg_data$dataset, function(ds) get_val(ds, "A1_pure_NMF_rate0"))
reg_data$delta <- reg_data$value - reg_data$base
reg_data$is_base <- reg_data$rate == "0 (pure NMF)"
reg_data$x_pos <- as.numeric(factor(reg_data$rate, levels = c("0 (pure NMF)", "0.1", "1", "10"))) +
  offset_map[reg_data$dataset]
reg_data$rate_lab <- factor(reg_data$rate, levels = c("0 (pure NMF)", "0.1", "1", "10"))

pA <- ggplot(reg_data, aes(x = x_pos)) +
  # 灰色基准条：所有柱子统一画（0 -> base）
  geom_rect(data = reg_data,
            aes(xmin = x_pos - w/2, xmax = x_pos + w/2, ymin = 0, ymax = base),
            fill = base_grey, color = NA, alpha = 0.35) +
  # 彩色增量柱（base -> value）
  geom_rect(data = subset(reg_data, delta > 1e-6),
            aes(xmin = x_pos - w/2, xmax = x_pos + w/2, ymin = base, ymax = value,
                fill = dataset), color = NA) +
  # 数值标注
  geom_text(data = subset(reg_data, !is_base),
            aes(x = x_pos, y = value, label = sprintf("%.4f", value)),
            vjust = -0.6, size = 2.4, fontface = "bold", family = "Arial", color = "black") +
  geom_text(data = subset(reg_data, is_base),
            aes(x = x_pos, y = base, label = sprintf("%.4f", base)),
            vjust = 1.6, size = 2.2, fontface = "bold", family = "Arial", color = base_grey) +
  scale_fill_manual(values = dataset_colors, name = "Dataset",
                  breaks = c("Mm_Brain", "Hm_Adipose", "Hm_Pancreas"),
                  labels = c("Simulation 1: Mm_Brain", "Simulation 2: Hm_Adipose", "Simulation 3: Hm_Pancreas")) +
  coord_cartesian(ylim = c(0.916, 0.972)) +
  scale_x_continuous(breaks = 1:4,
                     labels = c("0 (pure\nNMF)", "0.1", "1", "10"),
                     expand = expansion(mult = 0.08)) +
  labs(x = "Prior strength (rate)", y = "Pearson",
       title = "Component 1: elastic regularization contribution",
       subtitle = "grey bar = pure NMF baseline; colored bar height = gain") +
  sci_theme

ggsave(file.path(output_dir, "ABLATION_contribution_regularization.png"), pA,
       device = png, type = "cairo", width = 180/25.4, height = 90/25.4,
       units = "in", dpi = 600, bg = "white")

# ============================================================
# 图 B：marker 贡献（灰色基准 = random 基因）
# ============================================================
marker_levels <- c("Random genes", "DE top50", "ESMF-selected markers")
marker_label_map <- c("Random genes"          = "A2_random_genes",
                      "DE top50"              = "A2_DE_top50",
                      "ESMF-selected markers" = "A2_benchmark_marker")

marker_data <- expand.grid(marker = marker_levels, dataset = names(DS_DIR),
                           stringsAsFactors = FALSE)
marker_data$value <- mapply(function(ml, ds) get_val(ds, marker_label_map[[ml]]),
                            marker_data$marker, marker_data$dataset)
marker_data$base <- sapply(marker_data$dataset, function(ds) get_val(ds, "A2_random_genes"))
marker_data$delta <- marker_data$value - marker_data$base
marker_data$is_base <- marker_data$marker == "Random genes"
marker_data$x_pos <- as.numeric(factor(marker_data$marker,
                                       levels = c("Random genes", "DE top50", "ESMF-selected markers"))) +
  offset_map[marker_data$dataset]
marker_data$marker_lab <- factor(marker_data$marker,
                                 levels = c("Random genes", "DE top50", "ESMF-selected markers"))

pB <- ggplot(marker_data, aes(x = x_pos)) +
  # 灰色基准条：所有柱子统一画（0 -> base）
  geom_rect(data = marker_data,
            aes(xmin = x_pos - w/2, xmax = x_pos + w/2, ymin = 0, ymax = base),
            fill = base_grey, color = NA, alpha = 0.35) +
  geom_rect(data = subset(marker_data, delta > 1e-6),
            aes(xmin = x_pos - w/2, xmax = x_pos + w/2, ymin = base, ymax = value,
                fill = dataset), color = NA) +
  geom_text(data = subset(marker_data, !is_base),
            aes(x = x_pos, y = value, label = sprintf("%.4f", value)),
            vjust = -0.3, size = 2.4, fontface = "bold", family = "Arial", color = "black") +
  geom_text(data = subset(marker_data, is_base),
            aes(x = x_pos, y = base, label = sprintf("%.4f", base)),
            vjust = 1.6, size = 2.2, fontface = "bold", family = "Arial", color = base_grey) +
  scale_fill_manual(values = dataset_colors, name = "Dataset",
                  breaks = c("Mm_Brain", "Hm_Adipose", "Hm_Pancreas"),
                  labels = c("Simulation 1: Mm_Brain", "Simulation 2: Hm_Adipose", "Simulation 3: Hm_Pancreas")) +
  scale_y_continuous(expand = expansion(mult = c(0.08, 0.18))) +
  facet_wrap(~ factor(dataset,
                      levels = c("Mm_Brain", "Hm_Adipose", "Hm_Pancreas"),
                      labels = c("Simulation 1: Mm_Brain", "Simulation 2: Hm_Adipose", "Simulation 3: Hm_Pancreas")),
             scales = "free_y", nrow = 1) +
  scale_x_continuous(breaks = 1:3,
                     labels = c("Random\ngenes", "DE\ntop50", "ESMF-selected\nmarkers"),
                     expand = expansion(mult = 0.08)) +
  labs(x = "Marker panel", y = "Pearson",
       title = "Component 2: marker screening contribution",
       subtitle = "grey bar = random-gene baseline; colored bar height = gain") +
  sci_theme +
  theme(strip.text.x = element_blank())

ggsave(file.path(output_dir, "ABLATION_contribution_marker.png"), pB,
       device = png, type = "cairo", width = 200/25.4, height = 90/25.4,
       units = "in", dpi = 600, bg = "white")

# ============================================================
# 来源清单：与图片同目录保存，记录图上每个数值来自哪个结果文件
# ============================================================
prov_file <- file.path(output_dir, "ABLATION_contribution_provenance.txt")
con <- file(prov_file, open = "wt")
writeLines(c(
  "ABLATION contribution figures - data provenance",
  paste("Generated  :", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
  "Plot script: plot-ablation-contribution.R",
  "Figures    : ABLATION_contribution_regularization.png, ABLATION_contribution_marker.png",
  "Note       : all values are read directly from the source files listed below;",
  "             no number in the figures is hard-coded.",
  ""
), con)

for (ds in names(DS_DIR)) {
  writeLines(c(paste0("===== ", ds, " ====="),
               paste0("source: ", abl_src[[ds]]), ""), con)
  d <- abl[[ds]]
  keep_lb <- c(unname(rate_label_map), unname(marker_label_map))
  d <- d[d$label %in% keep_lb, c("label", "H_RMSE", "H_Pearson")]
  write.table(d, file = con, sep = "\t", quote = FALSE, row.names = FALSE)
  writeLines("", con)
}
close(con)

cat("两张贡献图已保存到", output_dir, "\n")
cat("  ABLATION_contribution_regularization.png\n")
cat("  ABLATION_contribution_marker.png\n")
cat("  ABLATION_contribution_provenance.txt\n")
