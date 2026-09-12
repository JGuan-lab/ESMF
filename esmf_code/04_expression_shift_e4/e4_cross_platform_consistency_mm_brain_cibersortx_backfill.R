# ============================================================================
# Ablation_Study_crossPlatformConsistency_Mm_Brain_E4_CIBERSORTx_backfill.R
#
# 用途：补算 E4 表达偏移实验（Mm_Brain）中 CIBERSORTx 的响应指标。
#
# 背景：
#   - E4 v3 脚本（Ablation_Study_crossPlatformConsistency_Mm_Brain_E4_ExpressionShift_v3.R）
#     已经把 CIBERSORTx 的输入（signature + mixture，场景 A 的 5 档 mixture、场景 B 的 5 档 signature）
#     导出到 results/Ablation_Study_crossPlatformConsistency_Mm_Brain_ParameterSensitivity_v2/CIBERSORTx_input_E4/，
#     但其「结果读取」是注释占位（第 273-286 行）。
#   - 现在 CIBERSORTx 的结果 CSV（E4_Brain_scenario{A,B}_shift<档>_Results.csv）已经跑出来，
#     本脚本按原始实验 cross-platform-consistency-mm-brain.R 的 CIBERSORTx 结算逻辑
#     （Section 10，第 645-692 行），逐档补算 H_RMSE / H_Pearson，
#     输出与 E4 v3 汇总表（out_df）相同格式的行，便于合并到 E4 结果中。
#
# 结算口径（严格对齐原始脚本 cross-platform-consistency-mm-brain.R）：
#   results <- read.table(csv, header=TRUE, sep=",")
#   rownums <- dim(results)[2] - 3                 # 细胞类型列数
#   H <- as.matrix(results[, 2:rownums])           # 取细胞类型比例列（跳过 Mixture / P-value / Correlation / RMSE）
#   rownames(H) <- results$Mixture; H <- t(H)      # 转置为 CT × mixture
#   H <- H[gtools::mixedsort(rownames(H)), ]       # 按细胞类型字典序排序
#   H <- reshape2::melt(H); colnames(H) <- c("CT","tissue","observed_values")
#   P <- as.matrix(generator[["P"]]); P <- P[mixedsort(rownames(P)), ]; P <- melt(P)
#   colnames(P) <- c("CT","tissue","expected_values")
#   INDICATOR <- merge(H, P)
#   RMSE / Pearson 计算同上（四舍五入 4 位）
# ============================================================================

setwd("/media/desk16/tjn050/LSC/LSC_deconvolution/LSC/")

# 依赖库（与原始脚本一致）
suppressPackageStartupMessages({
  library(dplyr)
  library(data.table)
  library(gtools)
  library(reshape2)
})

V2_SCRIPT <- "/media/desk16/tjn050/LSC/LSC_deconvolution/LSC/Ablation_Parameter_Codes/Ablation_Study_crossPlatformConsistency_Mm_Brain_ParameterSensitivity_v2.R"
RESULTS_DIR <- "/media/desk16/tjn050/LSC/LSC_deconvolution/LSC/results/Ablation_Study_crossPlatformConsistency_Mm_Brain_ParameterSensitivity_v2"
CIBX_DIR <- file.path(RESULTS_DIR, "CIBERSORTx_input_E4")

# ---- 1. 复用 v2 预处理，得到 ground truth P（与 E4 v3 完全一致）----
# 复用范围对齐 E4 v3：1:1015（P 已构造）+ 1034:1194（辅助函数），
# 保证 generator / P 与 E4 v3 跑 ESMF/bulk/MuSiC 时用的是同一份 ground truth。
cat(">> Sourcing v2 预处理（生成 ground truth P）...\n")
v2_lines <- readLines(V2_SCRIPT, warn = FALSE)
preproc_lines <- v2_lines[c(1:1015, 1034:1194)]
tmp_preproc <- tempfile(fileext = ".R")
writeLines(preproc_lines, tmp_preproc)
sys.source(tmp_preproc, envir = .GlobalEnv)
cat(">> 预处理完成，ground truth P 就绪。\n\n")

P <- generator[["P"]]  # 宽格式：细胞类型 × mixture（列名 mix1..mix1000）

# ---- 2. 档位 ----
shift_grid <- c(0, 0.1, 0.2, 0.5, 1)
sh_names <- gsub("\\.", "_", as.character(shift_grid))

# ---- 3. 结算函数（严格对齐原始脚本 CIBERSORTx 段）----
evaluate_cibersortx <- function(csv_file, P) {
  if (!file.exists(csv_file)) {
    stop("结果文件不存在: ", csv_file)
  }
  results <- read.table(csv_file, header = TRUE, sep = ",")
  rownums <- dim(results)[2] - 3
  H <- as.matrix(results[, 2:rownums])
  rownames(H) <- results$Mixture
  H <- t(H)
  H <- H[gtools::mixedsort(rownames(H)), , drop = FALSE]
  H <- reshape2::melt(H)
  colnames(H) <- c("CT", "tissue", "observed_values")

  Pm <- as.matrix(P)
  Pm <- Pm[gtools::mixedsort(rownames(Pm)), , drop = FALSE]
  Pm <- reshape2::melt(Pm)
  colnames(Pm) <- c("CT", "tissue", "expected_values")

  INDICATOR <- merge(H, Pm)
  INDICATOR$expected_values <- round(INDICATOR$expected_values, 2)
  INDICATOR$observed_values <- round(INDICATOR$observed_values, 2)

  INDICATOR <- INDICATOR %>% dplyr::summarise(
    RMSE = sqrt(mean((observed_values - expected_values)^2)) %>% round(., 4),
    Pearson = cor(observed_values, expected_values) %>% round(., 4)
  )
  c(H_RMSE = as.numeric(INDICATOR$RMSE), H_Pearson = as.numeric(INDICATOR$Pearson))
}

# ---- 4. 逐档补算 ----
rows <- list()
for (i in seq_along(shift_grid)) {
  sh <- shift_grid[i]
  sn <- sh_names[i]
  lab_sh <- as.character(sh)  # 用于 label（与 E4 v3 一致，label 里 shift 用原始数值文本）

  # 场景 A：目标加噪（signature 不变、mixture 5 档）
  fA <- file.path(CIBX_DIR, paste0("E4_Brain_scenarioA_shift", sn, "_Results.csv"))
  if (file.exists(fA)) {
    rA <- evaluate_cibersortx(fA, P)
    rows[[paste0("A_shift_", lab_sh)]] <- c(
      label = paste0("A_CIBERSORTx_shift_", lab_sh),
      scenario = "A", method = "CIBERSORTx", shift = sh,
      H_RMSE = unname(rA["H_RMSE"]), H_Pearson = unname(rA["H_Pearson"]))
    cat(sprintf("[A] shift=%s  CIBERSORTx  H_RMSE=%s  H_Pearson=%s\n",
                lab_sh, unname(rA["H_RMSE"]), unname(rA["H_Pearson"])))
  } else {
    cat(sprintf("[A] shift=%s  结果文件缺失，跳过：%s\n", lab_sh, fA))
  }

  # 场景 B：参考加偏（mixture 不变、signature 5 档）
  fB <- file.path(CIBX_DIR, paste0("E4_Brain_scenarioB_shift", sn, "_Results.csv"))
  if (file.exists(fB)) {
    rB <- evaluate_cibersortx(fB, P)
    rows[[paste0("B_shift_", lab_sh)]] <- c(
      label = paste0("B_CIBERSORTx_shift_", lab_sh),
      scenario = "B", method = "CIBERSORTx", shift = sh,
      H_RMSE = unname(rB["H_RMSE"]), H_Pearson = unname(rB["H_Pearson"]))
    cat(sprintf("[B] shift=%s  CIBERSORTx  H_RMSE=%s  H_Pearson=%s\n",
                lab_sh, unname(rB["H_RMSE"]), unname(rB["H_Pearson"])))
  } else {
    cat(sprintf("[B] shift=%s  结果文件缺失，跳过：%s\n", lab_sh, fB))
  }
}

# ---- 5. 汇总输出（格式与 E4 v3 的 out_df 对齐）----
out_df <- as.data.frame(do.call(rbind, rows), stringsAsFactors = FALSE)
for (cc in c("shift", "H_RMSE", "H_Pearson")) {
  out_df[[cc]] <- as.numeric(out_df[[cc]])
}
out_df <- out_df[order(out_df$scenario, out_df$shift), ]

current_datetime <- format(Sys.time(), "%Y%m%d_%H%M%S")
out_file <- file.path(RESULTS_DIR, paste0("E4_CIBERSORTx_backfill_", current_datetime, ".txt"))
write.table(out_df, file = out_file, sep = "\t", quote = FALSE,
            col.names = TRUE, row.names = FALSE)

cat("\n\n========== CIBERSORTx 补算汇总 ==========\n")
print(out_df)
cat("\n结果已保存到：", out_file, "\n")

# 与 E4 v3 汇总表列对齐，便于直接拼接：
#   E4 v3 输出列：label, iter, seed, H_RMSE, H_Pearson, W_RMSE, W_Pearson, shift, scenario, method
#   本脚本补齐 CIBERSORTx 后，可将 label/H_RMSE/H_Pearson/shift/scenario/method 合并进同一张表，
#   iter/seed 对 CIBERSORTx 无意义，可填 NA 或留空。
cat("\n提示：本表 label/scenario/method/shift/H_RMSE/H_Pearson 与 E4 v3 的 out_df 同名列一致，可直接 rbind 合并。\n")
