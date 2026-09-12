# ============================================================================
# Export CIBERSORTx Input for E4 Expression Shift (v2) —— Adipose
#
# 目的：导出 E4 偏移实验（场景 A 目标加噪 / 场景 B 参考加偏）所需的
#       CIBERSORTx 输入文件（signature matrix + mixture file），供手动网页运行。
# 逻辑：复用 Ablation_Study_MultiSourceReferenceIntegration_Hm_Adipose_ParameterSensitivity_v2.R 的预处理，拿到 T_orig / C_orig / eps_g / scale_* 后，
#       仅导出各 shift 档的输入文件即退出，不跑 10 个方法的偏移反卷积。
#
# 场景 A（目标加噪）：signature 不变 = C_orig；
#                     mixture 加噪 = T_orig + shift * scale_T * eps_g（5 档）
# 场景 B（参考加偏）：mixture 不变 = T_orig；
#                     signature 加偏 = C_orig + shift * scale_C * eps_c（5 档）
# ============================================================================

V2_SCRIPT <- "/media/desk16/tjn050/LSC/LSC_deconvolution/LSC/Ablation_Parameter_Codes/Ablation_Study_MultiSourceReferenceIntegration_Hm_Adipose_ParameterSensitivity_v2.R"
RESULTS_DIR <- "/media/desk16/tjn050/LSC/LSC_deconvolution/LSC/results/Ablation_Study_MultiSourceReferenceIntegration_Hm_Adipose_ParameterSensitivity_v2"

# ---- 1. 复用 v2 预处理 + 校验 + 函数（照抄原 E4 v2 脚本）----
v2_lines <- readLines(V2_SCRIPT, warn = FALSE)
preproc_lines <- v2_lines[c(1:1075, 1093:1247)]
tmp_preproc <- tempfile(fileext = ".R")
writeLines(preproc_lines, tmp_preproc)
cat(">> Sourcing v2 预处理 ...\n")
sys.source(tmp_preproc, envir = .GlobalEnv)
cat(">> 预处理完成，输入复现校验已通过。\n\n")

# ---- 2. 保存原始对象（照抄原 E4 v2 脚本）----
V_orig <- V
T_orig <- as.matrix(T)
sig_orig <- sig_matrix
C_orig <- as.matrix(C)

shift_grid <- c(0, 0.1, 0.2, 0.5, 1)

set.seed(20260818)
eps_g_full <- rnorm(nrow(T_orig), mean = 0, sd = 1)
names(eps_g_full) <- rownames(T_orig)

scale_V <- sd(as.numeric(as.matrix(V_orig)))
scale_T <- sd(as.numeric(as.matrix(T_orig)))
scale_sig <- sd(as.numeric(as.matrix(sig_orig)))
scale_C <- sd(as.numeric(as.matrix(C_orig)))

# ---- 3. 导出 CIBERSORTx 输入 ----
out_dir <- file.path(RESULTS_DIR, "CIBERSORTx_input_E4")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

write_ciber <- function(mat, file) {
  df <- data.frame(GeneSymbol = rownames(mat), as.matrix(mat), check.names = FALSE)
  write.table(df, file = file, sep = "\t", quote = FALSE,
              row.names = FALSE, col.names = TRUE)
}

# ---- 3.1 场景 A：目标加噪（signature 不变，mixture 5 档）----
sigA_file <- file.path(out_dir, "E4_Adipose_scenarioA_signature_matrix.txt")
write_ciber(C_orig, sigA_file)
cat("场景 A signature 已导出: ", basename(sigA_file),
    " (", nrow(C_orig), " genes x ", ncol(C_orig), " cell types)\n")

for (sh in shift_grid) {
  T_sh <- T_orig + sh * scale_T * eps_g_full
  f <- file.path(out_dir, sprintf("E4_Adipose_scenarioA_mixture_shift%s.txt",
                                   sub("\\.", "_", as.character(sh))))
  write_ciber(T_sh, f)
  cat(sprintf("  场景 A mixture shift=%s 已导出 (", sh),
      nrow(T_sh), " genes x ", ncol(T_sh), " samples)\n")
}

# ---- 3.2 场景 B：参考加偏（mixture 不变，signature 5 档）----
mixB_file <- file.path(out_dir, "E4_Adipose_scenarioB_mixture_file.txt")
write_ciber(T_orig, mixB_file)
cat("场景 B mixture 已导出: ", basename(mixB_file),
    " (", nrow(T_orig), " genes x ", ncol(T_orig), " samples)\n")

eps_c <- eps_g_full[rownames(C_orig)]
for (sh in shift_grid) {
  C_sh <- C_orig + sh * scale_C * eps_c
  f <- file.path(out_dir, sprintf("E4_Adipose_scenarioB_signature_shift%s.txt",
                                   sub("\\.", "_", as.character(sh))))
  write_ciber(C_sh, f)
  cat(sprintf("  场景 B signature shift=%s 已导出 (", sh),
      nrow(C_sh), " genes x ", ncol(C_sh), " cell types)\n")
}

cat("\n========================================\n")
cat("全部 CIBERSORTx 输入已保存到:\n  ", out_dir, "\n")
cat("========================================\n")
cat("用法说明：\n")
cat("  场景 A：用同一份 signature + 5 份不同 shift 的 mixture，跑 5 次；\n")
cat("  场景 B：用同一份 mixture + 5 份不同 shift 的 signature，跑 5 次。\n")
cat("  上传到 https://cibersortx.stanford.edu 手动运行。\n")
