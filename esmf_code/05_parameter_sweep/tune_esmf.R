# ============================================================================
# tune-esmf.R
#
# ESMF 默认参数网格搜索（tune_esmf）—— 与 Manuscript
#   "20260912-ESMF-SM-LSC.docx" 第二节 "2. Parameter Tuning Guidelines"
# 严格对齐的实现。
#
# 算法要点（来自文档）：
#   * 弹性正则系数矩阵 Theta (genes x cell-types)：
#       - 基线全为 1；
#       - 在每个细胞类型的 "calculated markers" 处被 P1 覆盖；
#       - 在每个细胞类型的 "prior markers"      处被 P2 覆盖
#         （两套 marker 重叠时 P2 优先）；
#       - 全局标量 rate 再整体缩放：effective_coeff = rate * Theta。
#         rate = 0 时系数为 0，退化为普通 NMF。
#   * 网格：global rate 取 8 个值；P1、P2 各取 4 个值（脚本传入的
#     para_grid 同时用于 P1 与 P2），笛卡尔积 8 x 4 x 4 = 128 组。
#   * 验证：5 折交叉验证，每折随机遮住 20% 的条目（genes x samples
#     矩阵中的条目），在剩余 80% 上拟合，重建被遮条目。
#   * 两个分数（均只在被遮条目上计算）：
#       分数1 = 被遮条目上 重建值 vs 真值 的 RMSE
#               （不使用参考谱的重建精度）；
#       分数2 = 推断细胞类型 signature (W) 与 参考谱 (sig_matrix)
#               的各细胞类型 Pearson 相关之均值
#               （参考谱恢复程度）。
#   * 选择：取 5 折平均分数2 最高（Pearson 最大）的参数组为最优。
#   * 另报告 "估计比例 vs 真实比例" 一致性（仅评估、不参与选择）。
#
# 调用约定（与 Tune_ESMF_*.R v2 驱动脚本一致）：
#   result <- tune_esmf(V, P, sig_matrix, sd_matrix, md, prio_markers,
#                       rate_grid, para_grid, iter_grid,
#                       iter_times = 1000, val_frac = 0.3,
#                       n_folds = 5, seed = 1233322, err_tol = 1e-6)
#   result$paras      : 128 行全组合结果（含各项分数）
#   result$best       : 1 行，平均 Pearson 最高的参数组
#   result$iter_scan  : 在最优参数下对 iter_times 单独扫描的收敛曲线
# ============================================================================

suppressMessages({
  require(CellMix, quietly = TRUE)   # 提供 rmatrix() 随机初始化
})

# ---------------------------------------------------------------------------
# 安全 Pearson：方差为零或样本不足时返回 NA，绝不中断网格搜索
# ---------------------------------------------------------------------------
safe_pearson <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 3) return(NA_real_)
  sx <- x[ok]; sy <- y[ok]
  if (sd(sx) == 0 || sd(sy) == 0) return(NA_real_)
  cor(sx, sy, method = "pearson")
}

# ---------------------------------------------------------------------------
# 构建 Theta 层矩阵（基线 1 / P1 于 calculated markers / P2 于 prior markers）
#   sig_matrix : genes x cell-types 参考谱（决定 Theta 维度与行列名）
#   md         : data.frame(gene, CT, ...)  calculated markers
#   prio_markers: list，元素 i 为第 i 个细胞类型的先验 marker 基因向量
#   p1, p2     : 两个层权重标量
# 返回 genes x cell-types 矩阵
# ---------------------------------------------------------------------------
build_theta <- function(sig_matrix, md, prio_markers, p1, p2) {
  genes <- rownames(sig_matrix)
  cts   <- colnames(sig_matrix)
  Theta <- matrix(1, nrow = length(genes), ncol = length(cts),
                  dimnames = list(genes, cts))

  # —— calculated markers -> P1 ——
  if (!is.null(md) && nrow(md) > 0) {
    mgenes <- as.character(md$gene)
    mct    <- as.character(md$CT)
    for (k in seq_along(mgenes)) {
      g <- mgenes[k]; c <- mct[k]
      if (!is.na(g) && g %in% genes && !is.na(c) && c %in% cts) {
        Theta[g, c] <- p1
      }
    }
  }

  # —— prior markers -> P2（重叠处覆盖 P1）——
  if (!is.null(prio_markers) && is.list(prio_markers) && length(prio_markers) > 0) {
    ct_names <- if (!is.null(names(prio_markers)) && all(nzchar(names(prio_markers)))) {
      names(prio_markers)
    } else if (!is.null(md) && nrow(md) > 0) {
      unique(as.character(md$CT))
    } else {
      cts
    }
    for (i in seq_along(prio_markers)) {
      c <- ct_names[min(i, length(ct_names))]
      if (is.null(c) || !(c %in% cts)) next
      for (g in as.character(prio_markers[[i]])) {
        if (!is.na(g) && g %in% genes) Theta[g, c] <- p2
      }
    }
  }

  Theta
}

# ---------------------------------------------------------------------------
# 单次参数组合的交叉验证：返回 5 折平均的各项分数
# ---------------------------------------------------------------------------
cv_one_combo <- function(rate_scalar, p1, p2,
                         V, P, sig_matrix, sd_matrix, md, prio_markers,
                         r, iter_times, err_tol, n_folds, seed, Theta_cache = NULL) {
  m <- nrow(V); n <- ncol(V)

  # 系数矩阵 = 全局 rate 缩放后的 Theta（rate=0 -> 全 0 -> 普通 NMF）
  Theta    <- build_theta(sig_matrix, md, prio_markers, p1, p2)
  rate_mat <- rate_scalar * Theta

  # —— 5 折条目级划分（每折 1/n_folds 的条目）——
  set.seed(seed)
  entry_order <- sample(m * n)
  fold_id     <- cut(seq_len(m * n),
                     breaks = seq(0, m * n, length.out = n_folds + 1),
                     labels = FALSE, include.lowest = TRUE)
  fold_of     <- integer(m * n)
  fold_of[entry_order] <- fold_id

  rmse_v  <- numeric(n_folds)
  pear_v  <- numeric(n_folds)
  ppear_v <- numeric(n_folds)

  for (f in seq_len(n_folds)) {
    pos   <- which(fold_of == f)
    ij    <- arrayInd(pos, dim = c(m, n))
    mask  <- matrix(FALSE, m, n)
    mask[ij] <- TRUE

    Vfit        <- V
    Vfit[mask]  <- 0                     # 遮住的条目按 0 处理（缺失）

    res <- suppressMessages(
      ESMF_deconvolution(V = Vfit, ml = NULL, r = r,
                        err_cutoff = err_tol, rate = rate_mat,
                        sig_matrix = sig_matrix, sd_matrix = sd_matrix,
                        iter_times = iter_times)
    )
    W <- res$W; H <- res$H
    recon <- W %*% H

    # 分数1：被遮条目 RMSE（不使用参考谱）
    rmse_v[f] <- sqrt(mean((V[mask] - recon[mask])^2, na.rm = TRUE))

    # 分数2：推断 signature 与参考谱 的各细胞类型 Pearson 均值
    pear_v[f] <- mean(vapply(seq_len(r), function(k)
      safe_pearson(W[, k], sig_matrix[, k]), numeric(1)), na.rm = TRUE)

    # 报告项：估计比例 vs 真实比例（仅评估，不参与选择）
    if (!is.null(P) && is.matrix(P) && nrow(P) == n && ncol(P) == r) {
      Hn <- sweep(H, 2, colSums(H), "/")   # 每列（样本）比例和为 1
      ppear_v[f] <- mean(vapply(seq_len(r), function(k)
        safe_pearson(Hn[k, ], P[, k]), numeric(1)), na.rm = TRUE)
    } else {
      ppear_v[f] <- NA_real_
    }
  }

  list(rmse       = mean(rmse_v,  na.rm = TRUE),
       pearson_sig = mean(pear_v,  na.rm = TRUE),
       prop_pearson = mean(ppear_v, na.rm = TRUE))
}

# ===========================================================================
# tune_esmf —— 主函数
# ===========================================================================
tune_esmf <- function(V, P = NULL,
                      sig_matrix, sd_matrix,
                      md = NULL, prio_markers = NULL,
                      rate_grid, para_grid, iter_grid,
                      iter_times = 1000,
                      val_frac = 0.3,        # 保留为接口兼容；官方 CV 用 n_folds 等折
                      n_folds = 5,           # 文档：5 折，每折 20% 条目
                      seed = 1233322, err_tol = 1e-6,
                      verbose = TRUE) {

  stopifnot(is.matrix(V) || is(V, "Matrix"))
  stopifnot(nrow(V) == nrow(sig_matrix))
  stopifnot(ncol(sig_matrix) == ncol(sd_matrix))
  r <- ncol(sig_matrix)

  # para_grid 同时用于 P1 与 P2；与 rate_grid 组成 128 组
  grid <- expand.grid(rate = rate_grid,
                      p1   = para_grid,
                      p2   = para_grid,
                      stringsAsFactors = FALSE)
  if (verbose) {
    cat(sprintf(">> 网格规模：%d 组 = %d rate x %d P1 x %d P2；5 折 CV；固定 iter_times=%d\n",
                nrow(grid), length(rate_grid), length(para_grid),
                length(para_grid), iter_times))
  }

  rows <- vector("list", nrow(grid))
  for (i in seq_len(nrow(grid))) {
    g  <- grid[i, ]
    cv <- cv_one_combo(rate_scalar = g$rate, p1 = g$p1, p2 = g$p2,
                       V = V, P = P, sig_matrix = sig_matrix, sd_matrix = sd_matrix,
                       md = md, prio_markers = prio_markers,
                       r = r, iter_times = iter_times, err_tol = err_tol,
                       n_folds = n_folds, seed = seed)
    rows[[i]] <- data.frame(rate = g$rate, p1 = g$p1, p2 = g$p2,
                            rmse = cv$rmse,
                            pearson_sig = cv$pearson_sig,
                            prop_pearson = cv$prop_pearson,
                            stringsAsFactors = FALSE)
    if (verbose && (i %% 16 == 0 || i == nrow(grid))) {
      cat(sprintf("   进度 %d/%d  当前(rate=%.3g,p1=%.3g,p2=%.3g): RMSE=%.4g Pearson=%.4f\n",
                  i, nrow(grid), g$rate, g$p1, g$p2, cv$rmse,
                  ifelse(is.na(cv$pearson_sig), NaN, cv$pearson_sig)))
    }
  }
  paras <- do.call(rbind, rows)

  # 选择：5 折平均 Pearson(signature, reference) 最高者
  best_idx <- which.max(paras$pearson_sig)
  best     <- paras[best_idx, , drop = FALSE]

  if (verbose) {
    cat("\n>> 最优参数（平均 Pearson 最高）：\n")
    print(best)
  }

  # —— iter_times 单独扫描（在最优 rate/P1/P2 下，单次全量拟合，看收敛）——
  if (verbose) cat("\n>> 扫描 iter_times（最优参数下）...\n")
  iter_rows <- vector("list", length(iter_grid))
  for (j in seq_along(iter_grid)) {
    it      <- iter_grid[j]
    Theta   <- build_theta(sig_matrix, md, prio_markers,
                           best$p1, best$p2)
    rate_mat <- best$rate * Theta
    res <- suppressMessages(
      ESMF_deconvolution(V = V, ml = NULL, r = r,
                        err_cutoff = err_tol, rate = rate_mat,
                        sig_matrix = sig_matrix, sd_matrix = sd_matrix,
                        iter_times = it)
    )
    W <- res$W
    pear <- mean(vapply(seq_len(r), function(k)
      safe_pearson(W[, k], sig_matrix[, k]), numeric(1)), na.rm = TRUE)
    iter_rows[[j]] <- data.frame(iter = it, err = res$err,
                                 pearson_sig = pear,
                                 stringsAsFactors = FALSE)
  }
  iter_scan <- do.call(rbind, iter_rows)

  list(paras = paras, best = best, iter_scan = iter_scan)
}
