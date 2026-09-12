ESMF_deconvolution <- function(V, ml, r, err_cutoff, rate, sig_matrix, sd_matrix, iter_times) {
  m <- dim(V)[1]
  n <- dim(V)[2]
  W <- rmatrix(m, r)
  H <- rmatrix(r, n)
  err_vc <- c()
  err <- 0.0
  max_dis_pre <- 0
  max_dis <- 0
  
  for (iteration in 1:iter_times) {
    V_pre <- W %*% H
    E <- V - V_pre
    err_pre <- err
    err <- 0.0
    
    for (i in seq(m)) {
      for (j in seq(n)) {
        err <- err + E[i, j] * E[i, j]
      }
    }
    
    err_vc <- c(err_vc, err)
    message(paste0("Iteration ", iteration, ": reconstruction error = ", format(err, scientific = TRUE, digits = 4)))
    
    dis <- (W - sig_matrix) / (sd_matrix + .Machine$double.eps)
    max_dis_pre <- max_dis
    max_dis <- max(abs(dis))
    message(paste0("Iteration ", iteration, ": maximum deviation = ", format(max_dis, scientific = TRUE, digits = 4)))
    message(paste0("Iteration ", iteration, ": number of significant deviations = ", length(which(abs(dis) > rate * sd_matrix))))
    
    if (err < err_cutoff) {
      message("Termination condition met (reached minimum error threshold).")
      RES <- list(W = W, H = H, err = err)
      return(RES)
    }
    
    if (err == err_pre) {
      message("Termination condition met (error stabilized).")
      RES <- list(W = W, H = H, err = err)
      return(RES)
    }
    
    if (err > err_pre) {
      message("Warning: reconstruction error increased from previous iteration.")
      
      if (iteration > 100) {
        message("Proceeding with termination after exceeding 100 iterations.")
        RES <- list(W = W, H = H, err = err)
        return(RES)
      }
    }
    
    if (err < err_pre) {
      message("Reconstruction error decreased from previous iteration.")
    }
    
    if (max_dis == max_dis_pre) {
      message("Maximum deviation stabilized.")
    }
    
    if (max_dis > max_dis_pre) {
      message("Maximum deviation increased from previous iteration.")
    }
    
    if (max_dis < max_dis_pre) {
      message("Maximum deviation decreased from previous iteration.")
    }
    
    VH <- V %*% t(H)
    WHH <- W %*% H %*% t(H)
    
    for (i_2 in seq(m)) {
      for (j_2 in seq(r)) {
        A <- VH[i_2, j_2] + rate[i_2, j_2] * sig_matrix[i_2, j_2]
        B <- WHH[i_2, j_2] + rate[i_2, j_2] * W[i_2, j_2]
        W[i_2, j_2] <- W[i_2, j_2] * (A / B)
      }
    }
    W <<- W
    
    WV <- t(W) %*% V
    WWH <- t(W) %*% W %*% H
    
    for (i_1 in seq(r)) {
      for (j_1 in seq(n)) {
        if (WWH[i_1, j_1] != 0) {
          H[i_1, j_1] <- H[i_1, j_1] %*% WV[i_1, j_1] / WWH[i_1, j_1]
        }
      }
    }
    H <<- H
    
    if (iter_times %% 10 == 0) {
      H[H < .Machine$double.eps] <- .Machine$double.eps
      W[W < .Machine$double.eps] <- .Machine$double.eps
    }
  }
  
  RES <- list(W = W, H = H, err = err)
  return(RES)
}