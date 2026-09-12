# ============================================================================
# sweep-esmf-mm-brain.R —— 小鼠脑数据 ESMF 超参数遍历（sweep，真实测试集评估）
#
# 做法：
#   1. source 主脚本到 SECTION 9 之前：
#      - 训练集构建参考对象（sig_matrix / sd_matrix / md / prio_markers）；
#      - 真实测试集（SRA866994_SRS4545964）单细胞经 Generator 合成伪 bulk，
#        得到评估用的 V（测试伪 bulk）与 P（测试真比例）。
#   2. 调用 esmf_sweep()：para1 固定=1，遍历 rate × para2 × para3，
#      在真实测试集 V/P 上用 PCC.C 选参。
# ============================================================================

cat("========================================\n")
cat("脚本版本：sweep-esmf-mm-brain.R (2026-09-09 sweep)\n")
cat("esmf_sweep.R 版本：2026-09-09 sweep（para1 固定=1，真实测试集 PCC.C 选参）\n")
cat("运行开始时间：", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("========================================\n\n")

setwd("/media/desk16/tjn050/LSC/LSC_deconvolution/LSC/")
library(NMF); library(CellMix); library(MASS); library(fitdistrplus)
library(scater); library(Matrix); library(Seurat); library(descend)
library(data.table); library(cowplot); library(rhdf5); library(clusterProfiler)
library(clue); library(gtools); library(reshape2); library(dplyr)
source("esmf_deconvolution.R")
source("esmf_helper_functions.R")
source("bulk_deconvolution.R")
source("normalize_col_in_matrix.R")
source("basic_functions.R")
source("cibersort.R")
source("esmf_sweep.R")

MAIN_SCRIPT <- "/media/desk16/tjn050/LSC/LSC_deconvolution/LSC/Ablation_Parameter_Codes/Ablation_Study_crossPlatformConsistency_Mm_Brain_ParameterSensitivity_v2.R"

# ---- 1. source 主脚本到 ESMF 输入就绪处（首个未注释的 md <- md[rownames(sig_matrix), ] 行）----
#     该截断点之前已完成：训练集参考 sig_matrix/sd_matrix 构建 + 真实测试集 V(伪 bulk)/P(真比例) 生成，
#     且 V/sig_matrix/sd_matrix/md 已按共有基因对齐（与稿件 ESMF 评估所用输入完全一致）。
main_lines <- readLines(MAIN_SCRIPT, warn = FALSE)
cut_lines <- grep("^md <- md\\[rownames\\(sig_matrix\\)", main_lines)
if (length(cut_lines) == 0) stop("未找到 ESMF 输入就绪截断点 (md <- md[rownames(sig_matrix), ])")
cut <- cut_lines[1]
preproc_lines <- main_lines[1:cut]
tmp_preproc <- tempfile(fileext = ".R")
writeLines(preproc_lines, tmp_preproc)
cat(">> Sourcing 主脚本（前", cut, "行，训练集参考 sig_matrix + 真实测试集 V/P 就绪）...\n")
sys.source(tmp_preproc, envir = .GlobalEnv)
cat(">> 完成。\n")
cat("   参考 sig_matrix:", nrow(sig_matrix), "基因 x", ncol(sig_matrix), "细胞类型（训练集 SRA866994_SRS4545961 构建）\n")
cat("   真实测试集 V:", nrow(V), "基因 x", ncol(V), "样本（测试集 SRA866994_SRS4545964 合成伪 bulk）\n")
cat("   真实测试集 P:", nrow(P), "细胞类型 x", ncol(P), "样本（真比例）\n\n")

# ---- 2. 调用 esmf_sweep（真实测试集评估）----
cat("\n>> 开始遍历（para1=1，rate x para2 x para3）...\n")
result <- esmf_sweep(
  V = V, P = P,
  sig_matrix = sig_matrix, sd_matrix = sd_matrix,
  md = md, prio_markers = prio_markers,
  rate_grid = c(0, 0.01, 0.1, 0.5, 1, 2, 5, 10),
  para_grid = c(0.1, 1, 10, 100),
  para1 = 1, iter_times = 3000, err_tol = 1e-6
)

# ---- 3. 输出结果 ----
results_dir <- "/media/desk16/tjn050/LSC/LSC_deconvolution/LSC/results/Ablation_Study_crossPlatformConsistency_Mm_Brain_ParameterSensitivity_v2"
current_datetime <- format(Sys.time(), "%Y%m%d_%H%M%S")

cat("\n\n========== 选参依据：", result$select_by, "（组合总数:", result$n_combos, "）==========\n")
print(result$best)

out_file <- file.path(results_dir, paste0("esmf_sweep_paras_", current_datetime, ".txt"))
write.table(result$paras, file = out_file, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
cat("\n全组合结果已保存到:", out_file, "\n")

best_file <- file.path(results_dir, paste0("esmf_sweep_best_", current_datetime, ".txt"))
write.table(result$best, file = best_file, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
cat("最优参数已保存到:", best_file, "\n")
