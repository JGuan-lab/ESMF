# ============================================================================
# Ablation_Study_GeneralizationFromSingleReference_Hm_Pancreas_E4_ExpressionShift_v3.R
#
# E4（v3）：表达偏移模拟——双场景（回应审稿人 R2 M2③）
# v3 相对 v2 的改动：
#   1. 统一偏移口径：所有方法（ESMF/bulk/MuSiC）的偏移矩阵均做
#      非负截断 + 四舍五入：pmax(0, round(x + shift*scale*eps, 3))，
#      解决"加性偏移产生负表达值/15 位长小数"的问题（CIBERSORTx 输入要求非负）。
#   2. 新增 MuSiC（独立循环，严格沿用原脚本 CPM 处理链 + melt/merge 结算）：
#      场景 A 对 bulk.mtx(T) 加噪；场景 B 对单细胞参考 C 加偏。
#   3. CIBERSORTx：输入文件导出逻辑直接内嵌；结果读取留注释占位，
#      结果文件保存到输入文件夹后，重跑本脚本即自动读入并结算。
#
# 场景 A（目标加噪）：V/T 加噪，参考 Σ/C 不变。
# 场景 B（参考加偏）：sig_matrix/C（及 MuSiC 单细胞参考）加偏，目标 V/T 不变。
# ============================================================================

V2_SCRIPT <- "/media/desk16/tjn050/LSC/LSC_deconvolution/LSC/Ablation_Parameter_Codes/Ablation_Study_GeneralizationFromSingleReference_Hm_Pancreas_ParameterSensitivity_v2.R"
RESULTS_DIR <- "/media/desk16/tjn050/LSC/LSC_deconvolution/LSC/results/Ablation_Study_GeneralizationFromSingleReference_Hm_Pancreas_ParameterSensitivity_v2"

# ---- 1. 复用 v2 预处理 + 校验 + 函数 ----
v2_lines <- readLines(V2_SCRIPT, warn = FALSE)
preproc_lines <- v2_lines[c(1:1151, 1153:1240, 1596:1716)]
tmp_preproc <- tempfile(fileext = ".R")
writeLines(preproc_lines, tmp_preproc)
cat(">> Sourcing v2 预处理 ...\n")
sys.source(tmp_preproc, envir = .GlobalEnv)
cat(">> 预处理完成，输入复现校验已通过。\n\n")

# ---- 2. 保存原始对象 ----
V_orig <- V
T_orig <- as.matrix(T)
sig_orig <- sig_matrix
C_orig <- as.matrix(C)

# MuSiC 单细胞参考（沿用原脚本处理链：Transformation(none) + Scaling(TMM) → CPM）
C_sc_orig <- Scaling(Transformation(train_cellID, transformation), normalization_scC)
C_sc_orig <- as.matrix(C_sc_orig)
# MuSiC 的 bulk.mtx（沿用原脚本：Transformation(none) + Scaling(TMM) → CPM）
T_music_orig <- Scaling(Transformation(generator[["T"]], transformation), normalization_scT)
T_music_orig <- as.matrix(T_music_orig)
# MuSiC 的 cellType 标注（沿用原脚本第 757-767 行 phenoDataC 处理：改名 + rownames = cellID）
pDataC_music <- pDataC
if (length(grep("[N-n]ame", colnames(pDataC_music))) > 0) {
  sample_column <- grep("[N-n]ame", colnames(pDataC_music))
} else {
  sample_column <- grep("[S-s]ample|[S-s]ubject", colnames(pDataC_music))
}
colnames(pDataC_music)[sample_column] <- "SubjectName"
rownames(pDataC_music) <- pDataC_music$cellID

# ---- 3. 档位与对照 ----
shift_grid <- c(0, 0.1, 0.2, 0.5, 1)
bulk_ctl_methods <- c("nnls", "FARDEEP", "RLR", "DCQ", "elastic_net", "lasso", "ridge", "OLS", "EPIC", "DSA")

# 噪声方向（基因特异性，固定 seed，跨档复用）
set.seed(20260818)
eps_g_full <- rnorm(nrow(T_orig), mean = 0, sd = 1)
names(eps_g_full) <- rownames(T_orig)

# 噪声尺度（各自矩阵的 SD，使 shift 可解释）
scale_V <- sd(as.numeric(as.matrix(V_orig)))
scale_T <- sd(as.numeric(as.matrix(T_orig)))
scale_sig <- sd(as.numeric(as.matrix(sig_orig)))
scale_C <- sd(as.numeric(as.matrix(C_orig)))
scale_C_sc <- sd(as.numeric(as.matrix(C_sc_orig)))
scale_T_music <- sd(as.numeric(as.matrix(T_music_orig)))

# 统一偏移函数：只对非零位置加噪 + 非负截断 + 四舍五入（保留矩阵维度）
# 说明：0 位置保持 0（不被噪声污染），非零位置加噪并截断为 2 位小数；
#       这样既保证全非负（CIBERSORTx 友好）、又避免全 0 行（NMF/DSA 报错）、
#       且不会像全矩阵加噪那样把 0 填成噪声导致文件暴增。
shift_nonneg <- function(x, sh, scale, eps) {
  # eps 是命名向量（按基因）；按 x 的行名对齐后广播成与 x 同维的矩阵，
  # 避免直接用 eps[nz] 因矩阵逻辑索引与向量长度不匹配产生 NA。
  eps_row <- eps[rownames(x)]
  eps_mat <- matrix(eps_row, nrow = nrow(x), ncol = ncol(x))
  y <- x
  nz <- x > 0
  y[nz] <- pmax(0, round(x[nz] + sh * scale * eps_mat[nz], 1))
  # 全 0 行保护：强负噪声可能把某低表达基因的所有非零值都清零，导致整行全 0，
  # 进而使 DSA（CellMix::ged/NMF）报 "null row" 错误。这里恢复该行原值（原值非负）。
  zero_rows <- which(rowSums(y) == 0)
  if (length(zero_rows) > 0) {
    y[zero_rows, ] <- x[zero_rows, ]
  }
  dim(y) <- dim(x)
  dimnames(y) <- dimnames(x)
  y
}

cat("========== E4 v3: EXPRESSION SHIFT (additive noise, non-neg truncated, 2 scenarios) ==========\n")

results <- list()

# ---- 4. 场景 A：目标加性噪声 ----
cat("\n########## 场景 A：目标数据加性噪声 ##########\n")
for (sh in shift_grid) {
  cat(sprintf("\n-- A shift = %s --\n", sh))

  # A.1 ESMF
  eps_v <- eps_g_full[rownames(V_orig)]
  V <<- shift_nonneg(V_orig, sh, scale_V, eps_v)
  r <- evaluate_esmf(rate_mat = "coef", md_obj = md, para = c(1, 1, 1),
                     iter_times = 3000, seed_val = 1233322,
                     label = paste0("A_ESMF_shift_", sh))
  V <<- V_orig
  results[[paste0("A_ESMF_shift_", sh)]] <- r

  # A.2 bulk 对照（counts 尺度，T 加噪、C 不变）
  T_shifted <- shift_nonneg(T_orig, sh, scale_T, eps_g_full)
  for (mth in bulk_ctl_methods) {
    set.seed(123456)
    RES <- Deconvolution(T = T_shifted, C = C_orig, method = mth, P = P,
                         elem = to_remove, marker_distrib = marker_distrib,
                         refProfiles.var = refProfiles.var)
    RES <- RES %>% dplyr::summarise(
      RMSE = sqrt(mean((observed_values - expected_values)^2)) %>% round(., 4),
      Pearson = cor(observed_values, expected_values) %>% round(., 4))
    results[[paste0("A_", mth, "_shift_", sh)]] <- c(
      label = paste0("A_", mth, "_shift_", sh), iter = 3000, seed = 1233322,
      H_RMSE = as.numeric(RES$RMSE), H_Pearson = as.numeric(RES$Pearson),
      W_RMSE = NA, W_Pearson = NA)
  }
}

# ---- 5. 场景 B：参考加偏移 ----
cat("\n########## 场景 B：参考数据加偏移 ##########\n")
for (sh in shift_grid) {
  cat(sprintf("\n-- B shift = %s --\n", sh))

  # B.1 ESMF
  eps_sig <- eps_g_full[rownames(sig_orig)]
  sig_matrix <<- shift_nonneg(sig_orig, sh, scale_sig, eps_sig)
  r <- evaluate_esmf(rate_mat = "coef", md_obj = md, para = c(1, 1, 1),
                     iter_times = 3000, seed_val = 1233322,
                     label = paste0("B_ESMF_shift_", sh))
  sig_matrix <<- sig_orig
  results[[paste0("B_ESMF_shift_", sh)]] <- r

  # B.2 bulk 对照（counts 尺度，C 加偏、T 不变）
  eps_c <- eps_g_full[rownames(C_orig)]
  C_shifted <- shift_nonneg(C_orig, sh, scale_C, eps_c)
  for (mth in bulk_ctl_methods) {
    set.seed(123456)
    RES <- Deconvolution(T = T_orig, C = C_shifted, method = mth, P = P,
                         elem = to_remove, marker_distrib = marker_distrib,
                         refProfiles.var = refProfiles.var)
    RES <- RES %>% dplyr::summarise(
      RMSE = sqrt(mean((observed_values - expected_values)^2)) %>% round(., 4),
      Pearson = cor(observed_values, expected_values) %>% round(., 4))
    results[[paste0("B_", mth, "_shift_", sh)]] <- c(
      label = paste0("B_", mth, "_shift_", sh), iter = 3000, seed = 1233322,
      H_RMSE = as.numeric(RES$RMSE), H_Pearson = as.numeric(RES$Pearson),
      W_RMSE = NA, W_Pearson = NA)
  }
}

# ---- 6. MuSiC（独立循环，严格沿用原脚本处理链 + melt/merge 结算）----
cat("\n########## MuSiC（独立循环）##########\n")
# 重新取宽格式 P（原脚本 MuSiC 段第 739 行 P <- generator[["P"]]，宽格式）
P_music <- generator[["P"]]

for (sh in shift_grid) {

  # ---- 场景 A：bulk T（CPM）加噪，单细胞参考不变 ----
  eps_tm <- eps_g_full[rownames(T_music_orig)]
  T_music_A <- shift_nonneg(T_music_orig, sh, scale_T_music, eps_tm)
  C.sce_A <- SingleCellExperiment(assays = list(counts = C_sc_orig), colData = pDataC_music)
  set.seed(123456)
  RES_TEMP_A <- t(MuSiC::music_prop(bulk.mtx = T_music_A, sc.sce = C.sce_A,
                                    clusters = 'cellType', samples = 'cellID',
                                    markers = NULL, normalize = FALSE,
                                    verbose = FALSE)$Est.prop.weighted)

  # 沿用原脚本结算（melt + merge）
  RESULTS_A <- RES_TEMP_A[gtools::mixedsort(rownames(RES_TEMP_A)), , drop = FALSE]
  RESULTS_A <- RESULTS_A[gtools::mixedsort(rownames(RESULTS_A)), , drop = FALSE]
  RESULTS_A <- as.data.table(RESULTS_A, keep.rownames = TRUE)
  RESULTS_A <- data.table::melt(RESULTS_A)
  colnames(RESULTS_A) <- c("CT", "tissue", "observed_values")

  P_A <- P_music[gtools::mixedsort(rownames(P_music)), , drop = FALSE]
  P_A$CT <- rownames(P_A)
  setDT(P_A)
  P_A <- data.table::melt(P_A, id.vars = "CT")
  colnames(P_A) <- c("CT", "tissue", "expected_values")

  RESULTS_A <- merge(RESULTS_A, P_A)
  RESULTS_A$expected_values <- round(RESULTS_A$expected_values, 3)
  RESULTS_A$observed_values <- round(RESULTS_A$observed_values, 3)
  IND_A <- RESULTS_A %>% dplyr::summarise(
    RMSE = sqrt(mean((observed_values - expected_values)^2)) %>% round(., 4),
    Pearson = cor(observed_values, expected_values) %>% round(., 4))
  results[[paste0("A_MuSiC_shift_", sh)]] <- c(
    label = paste0("A_MuSiC_shift_", sh), iter = 3000, seed = 1233322,
    H_RMSE = as.numeric(IND_A$RMSE), H_Pearson = as.numeric(IND_A$Pearson),
    W_RMSE = NA, W_Pearson = NA)
  cat(sprintf("  A MuSiC shift=%s 完成 (H_Pearson=%.4f)\n", sh, as.numeric(IND_A$Pearson)))

  # ---- 场景 B：单细胞参考 C（CPM）加偏，bulk T 不变 ----
  eps_csc <- eps_g_full[rownames(C_sc_orig)]
  C_sc_shifted <- shift_nonneg(C_sc_orig, sh, scale_C_sc, eps_csc)
  C.sce_B <- SingleCellExperiment(assays = list(counts = C_sc_shifted), colData = pDataC_music)
  set.seed(123456)
  RES_TEMP_B <- t(MuSiC::music_prop(bulk.mtx = T_music_orig, sc.sce = C.sce_B,
                                    clusters = 'cellType', samples = 'cellID',
                                    markers = NULL, normalize = FALSE,
                                    verbose = FALSE)$Est.prop.weighted)

  RESULTS_B <- RES_TEMP_B[gtools::mixedsort(rownames(RES_TEMP_B)), , drop = FALSE]
  RESULTS_B <- RESULTS_B[gtools::mixedsort(rownames(RESULTS_B)), , drop = FALSE]
  RESULTS_B <- as.data.table(RESULTS_B, keep.rownames = TRUE)
  RESULTS_B <- data.table::melt(RESULTS_B)
  colnames(RESULTS_B) <- c("CT", "tissue", "observed_values")

  P_B <- P_music[gtools::mixedsort(rownames(P_music)), , drop = FALSE]
  P_B$CT <- rownames(P_B)
  setDT(P_B)
  P_B <- data.table::melt(P_B, id.vars = "CT")
  colnames(P_B) <- c("CT", "tissue", "expected_values")

  RESULTS_B <- merge(RESULTS_B, P_B)
  RESULTS_B$expected_values <- round(RESULTS_B$expected_values, 3)
  RESULTS_B$observed_values <- round(RESULTS_B$observed_values, 3)
  IND_B <- RESULTS_B %>% dplyr::summarise(
    RMSE = sqrt(mean((observed_values - expected_values)^2)) %>% round(., 4),
    Pearson = cor(observed_values, expected_values) %>% round(., 4))
  results[[paste0("B_MuSiC_shift_", sh)]] <- c(
    label = paste0("B_MuSiC_shift_", sh), iter = 3000, seed = 1233322,
    H_RMSE = as.numeric(IND_B$RMSE), H_Pearson = as.numeric(IND_B$Pearson),
    W_RMSE = NA, W_Pearson = NA)
  cat(sprintf("  B MuSiC shift=%s 完成 (H_Pearson=%.4f)\n", sh, as.numeric(IND_B$Pearson)))
}

# ---- 7. CIBERSORTx（输入导出 + 结果读取占位）----
# 7.1 导出 CIBERSORTx 输入文件（signature + mixture，含各 shift 档）
cibx_dir <- file.path(RESULTS_DIR, "CIBERSORTx_input_E4")
dir.create(cibx_dir, showWarnings = FALSE, recursive = TRUE)

write_ciber <- function(mat, file) {
  df <- data.frame(GeneSymbol = rownames(mat), as.matrix(mat), check.names = FALSE)
  write.table(df, file = file, sep = "\t", quote = FALSE,
              row.names = FALSE, col.names = TRUE)
}

# 场景 A：signature 不变（C_orig），mixture 5 档（T 加噪）
write_ciber(C_orig, file.path(cibx_dir, "E4_Pancreas_scenarioA_signature_matrix.txt"))
# 场景 B：mixture 不变（T_orig），signature 5 档（C 加偏）
write_ciber(T_orig, file.path(cibx_dir, "E4_Pancreas_scenarioB_mixture_file.txt"))
eps_c_cibx <- eps_g_full[rownames(C_orig)]

for (sh in shift_grid) {
  sh_name <- gsub("\\.", "_", as.character(sh))
  # 场景 A mixture（T 加噪）
  T_A <- shift_nonneg(T_orig, sh, scale_T, eps_g_full)
  write_ciber(T_A, file.path(cibx_dir, paste0("E4_Pancreas_scenarioA_mixture_shift", sh_name, ".txt")))
  # 场景 B signature（C 加偏）
  C_B <- shift_nonneg(C_orig, sh, scale_C, eps_c_cibx)
  write_ciber(C_B, file.path(cibx_dir, paste0("E4_Pancreas_scenarioB_signature_shift", sh_name, ".txt")))
}
cat("CIBERSORTx 输入已导出到:", cibx_dir, "\n")

# 7.2 结果读取占位
# 说明：手动上传上述文件到 https://cibersortx.stanford.edu 运行后，
#       把下载的 Results 文件保存到 cibx_dir（本输入文件夹），
#       命名约定：
#         E4_Pancreas_scenarioA_shift<档>_Results.txt
#         E4_Pancreas_scenarioB_shift<档>_Results.txt
#       保存后重跑本脚本，下面的逻辑会自动读入并结算、加入 results。
cat("\n[CIBERSORTx 占位] 结果文件就绪后，请取消下方注释并补全结算逻辑。\n")
cat("结果文件应保存到:", cibx_dir, "\n")

# 结果读取（占位，待结果文件就位后启用）：
# for (sh in shift_grid) {
#   sh_name <- gsub("\\.", "_", as.character(sh))
#   # 场景 A
#   fA <- file.path(cibx_dir, paste0("E4_Pancreas_scenarioA_shift", sh_name, "_Results.txt"))
#   if (file.exists(fA)) {
#     resA <- read.table(fA, header = TRUE, sep = "\t")
#     rn <- dim(resA)[2] - 3
#     H_A <- as.matrix(resA[, 2:rn]); rownames(H_A) <- resA$Mixture; H_A <- t(H_A)
#     H_A <- H_A[gtools::mixedsort(rownames(H_A)), , drop = FALSE]
#     # ... 结算同 MuSiC（melt + merge），写入 results[[paste0("A_CIBERSORTx_shift_", sh)]]
#   }
#   # 场景 B（同 A，读 B 文件）
# }

# ---- 8. 汇总写出 ----
out_df <- as.data.frame(do.call(rbind, results), stringsAsFactors = FALSE)
for (cc in c("H_RMSE", "H_Pearson", "W_RMSE", "W_Pearson")) {
  out_df[[cc]] <- as.numeric(out_df[[cc]])
}
out_df$shift <- as.numeric(sub(".*_shift_", "", out_df$label))
out_df$scenario <- sub("_.*", "", out_df$label)
out_df$method <- sub("^[AB]_", "", out_df$label)
out_df$method <- sub("_shift_.*", "", out_df$method)
out_df <- out_df[order(out_df$scenario, out_df$method, out_df$shift), ]

current_datetime <- format(Sys.time(), "%Y%m%d_%H%M%S")
out_file <- file.path(RESULTS_DIR, paste0("E4_expression_shift_v3_", current_datetime, ".txt"))
write.table(out_df, file = out_file, sep = "\t", quote = FALSE,
            col.names = TRUE, row.names = FALSE)

cat("\n\n========== E4 v3 汇总 ==========\n")
print(out_df[, c("scenario", "method", "shift", "H_RMSE", "H_Pearson")])
cat("\n结果已保存到：", out_file, "\n")
