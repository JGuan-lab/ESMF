# ============================================================================
# Ablation_Study_MultiSourceReferenceIntegration_Hm_Adipose_E4_CIBERSORTx_backfill.R
#
# 用途：补算 E4 表达偏移实验（Hm_Adipose，人脂肪）中 CIBERSORTx 的响应指标。
#
# 背景：
#   - E4 v3 脚本（Ablation_Study_MultiSourceReferenceIntegration_Hm_Adipose_E4_ExpressionShift_v3.R）
#     已把 CIBERSORTx 输入导出到 CIBERSORTx_input_E4/，结果读取是注释占位。
#   - 现在 CIBERSORTx 结果 CSV 已跑出，本脚本按原始实验的 CIBERSORTx 结算逻辑
#     逐档补算 H_RMSE / H_Pearson。
#
# 细胞类型：Fibroblast, Mono_Macro, CD4T（3 种）
#
# 结算口径（严格对齐原始脚本 CIBERSORTx 段，同 Mm_Brain backfill）：
#   rownums <- dim[2]-3；H <- results[,2:rownums] 转置；mixedsort；melt；merge(P)；算 RMSE/Pearson
# ============================================================================

setwd("/media/desk16/tjn050/LSC/LSC_deconvolution/LSC/")

suppressPackageStartupMessages({
  library(dplyr)
  library(data.table)
  library(gtools)
  library(reshape2)
})

V2_SCRIPT <- "/media/desk16/tjn050/LSC/LSC_deconvolution/LSC/Ablation_Parameter_Codes/Ablation_Study_MultiSourceReferenceIntegration_Hm_Adipose_ParameterSensitivity_v2.R"
RESULTS_DIR <- "/media/desk16/tjn050/LSC/LSC_deconvolution/LSC/results/Ablation_Study_MultiSourceReferenceIntegration_Hm_Adipose_ParameterSensitivity_v2"
CIBX_DIR <- file.path(RESULTS_DIR, "CIBERSORTx_input_E4")

# ---- 1. 复用 v2 预处理（复用范围严格对齐 E4 v3 脚本），得到 ground truth P ----
cat(">> Sourcing v2 预处理（生成 ground truth P）...\n")
v2_lines <- readLines(V2_SCRIPT, warn = FALSE)
preproc_lines <- v2_lines[c(1:1075, 1093:1247)]
tmp_preproc <- tempfile(fileext = ".R")
writeLines(preproc_lines, tmp_preproc)
sys.source(tmp_preproc, envir = .GlobalEnv)
cat(">> 预处理完成，ground truth P 就绪。\n\n")

P <- generator[["P"]]

# ---- 2. 档位 ----
shift_grid <- c(0, 0.1, 0.2, 0.5, 1)
sh_names <- gsub("\\.", "_", as.character(shift_grid))

# ---- 3. 结算函数 ----
evaluate_cibersortx <- function(csv_file, P) {
  if (!file.exists(csv_file)) stop("结果文件不存在: ", csv_file)
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
  lab_sh <- as.character(sh)

  fA <- file.path(CIBX_DIR, paste0("E4_Adipose_scenarioA_shift", sn, "_Results.csv"))
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

  fB <- file.path(CIBX_DIR, paste0("E4_Adipose_scenarioB_shift", sn, "_Results.csv"))
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

# ---- 5. 汇总输出 ----
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
