# ============================================================================
# list-package-versions.R
#
# 输出 Supplementary Table 5 中各去卷积方法所用 R 包 / 软件的版本号，
# 用于核对与补全表中 "Package (Version)" 列。
#
# 用法: Rscript list-package-versions.R
# ============================================================================

# Table 5 中登记的包名（按表内方法顺序）。
#   NA 表示该行并非 CRAN/Bioconductor 包：
#     ESMF       —— 自定义实现
#     CIBERSORTx —— 外部网页 / Docker 平台，无 R 包版本号，需另行记录
pkg_table <- c(
  ESMF        = NA,
  nnls        = "nnls",
  FARDEEP     = "FARDEEP",
  RLR         = "MASS",
  DCQ         = "ComICS",
  elastic_net = "glmnet",
  lasso       = "glmnet",
  ridge       = "glmnet",
  OLS         = "stats",
  EPIC        = "EPIC",
  DSA         = "CellMix",
  MuSiC       = "MuSiC",
  CIBERSORTx  = NA
)

# 支撑包：非去卷积方法本身，但影响结果可复现性，一并记录
support_pkgs <- c(
  "SingleCellExperiment", "Biobase", "SummarizedExperiment",
  "Matrix", "matrixStats", "limma", "edgeR", "DESeq2",
  "fitdistrplus", "AUCell", "GSEABase",
  "dplyr", "tidyr", "reshape2", "data.table", "ggplot2", "cowplot", "patchwork"
)

ver <- function(p) {
  if (is.na(p)) return("-")
  if (requireNamespace(p, quietly = TRUE)) {
    as.character(packageVersion(p))
  } else {
    "NOT INSTALLED"
  }
}

pad <- function(x, n) formatC(x, width = n, flag = "-")

cat("========================================\n")
cat("R / 平台\n")
cat("========================================\n")
cat(R.version.string, "\n")
cat("Platform :", R.version$platform, "\n")
cat("Locale   :", Sys.getlocale("LC_CTYPE"), "\n\n")

cat("========================================\n")
cat("Supplementary Table 5 · Package (Version)\n")
cat("========================================\n")
cat(pad("Method", 13), pad("Package", 16), "Version\n")
cat(strrep("-", 48), "\n")
for (m in names(pkg_table)) {
  pkg <- pkg_table[[m]]
  label <- if (is.na(pkg)) "(custom / external)" else pkg
  cat(pad(m, 13), pad(label, 16), ver(pkg), "\n")
}

cat("\n========================================\n")
cat("支撑包版本\n")
cat("========================================\n")
for (p in support_pkgs) cat(pad(p, 24), ver(p), "\n")

cat("\n========================================\n")
cat("说明\n")
cat("========================================\n")
cat("ESMF        : 自定义实现，无独立版本号；随 LSC 代码库发布。\n")
cat("CIBERSORTx  : 经 https://cibersortx.stanford.edu 运行，\n")
cat("              无 R 包版本；投稿时请记录平台访问日期，\n")
cat("              若用 Docker 镜像则记录镜像版本。\n")
cat("OLS         : stats 为基础包，版本与 R 一致。\n")
