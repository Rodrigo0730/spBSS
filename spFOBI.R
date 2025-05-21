library(SpatialBSS)
library(JADE)
library(spGARCH)
library(spdep)
library(sp)
library(dplyr)
library(sf)
library(moments)

set.seed(123)

sim_FOBI <- function(d, reps=100) {
  p <- 3
  n <- d^2
  md_values <- numeric(reps)
  for (r in 1:reps){
    coords <- matrix(runif(n * 2, 0, n), ncol = 2)
    dist_matrix <- as.matrix(dist(coords))
    
    W <- exp(-dist_matrix)
    diag(W) <- 0
    W <- W / rowSums(W)
    
    #generate data from log-spARCH model for heteroscedastic
    rho1 <- 0.5
    rho2 <- 0.2
    alpha1 <- 1
    alpha2 <- 0.5
    sink(tempfile())
    source1 <- sim.spARCH(n = n, rho = rho1, alpha = alpha1, W = W, type = "log-spARCH")
    source2 <- sim.spARCH(n = n, rho = rho2, alpha = alpha2, W = W, type = "log-spARCH")
    source3 <- sim.spARCH(n = n, rho = rho1, alpha = alpha2, W = W, type = "log-spARCH")
    sink()
    noise <- matrix(rt(n * p, df = 5), nrow = n, ncol = p)
    sources <- cbind(source1, source2, source3) 
    
    #True mixing matrix
    A <- matrix(runif(p^2, 0, 1), p, p)
    
    #generate the observed data, whiten the data
    X <- sources %*% t(A)
    X_centered <- scale(X, center = TRUE, scale=FALSE)
    cov_X <- cov(X_centered)
    E <- eigen(cov_X)
    whitening_mat <- E$vectors %*% diag(1 / sqrt(E$values)) %*% t(E$vectors)
    X_white <- X_centered %*% whitening_mat
    
    kernel_bandwidths <- c(0.05, 0.1, 0.25) * n
    kernels <- spatial_kernel_matrix(coords, kernel_type = "gauss", kernel_parameters = kernel_bandwidths)
    
    B_f_list <- array(0, dim = c(p, p, length(kernel_bandwidths)))
    for (i in seq_along(kernels)) {
      norms_sq <- rowSums(X_white^2)   
      Kmat <- kernels[[i]]
      W1 <- Kmat * matrix(norms_sq, nrow = n, ncol = n, byrow = TRUE) / n
      Bf <- t(X_white) %*% W1 %*% X_white
      B_f_list[, , i] <- Bf
    }
    tryCatch({
      jd_res <- rjd(B_f_list, maxiter=10000)
    }, error = function(e) {
      cat("Error in frjd: ", e$message, "\n")
      browser()
    }, finally = {
    })
    U <- jd_res$V
    
    W_est <- t(U) %*% whitening_mat
    md_values[r] <- MD(W_est, A)
    if (r %% 5 == 0) {
      cat("completed", r, "out of", reps, "simulations\n")
    }
  }
  
  return(list(mean_md = mean(md_values), d = d, reps = reps, n = n))
  
}

res <- sim_FOBI(d = 80, reps = 50)
message("mean MD over 20 simulations: ", res$mean_md, "\n")

  
  
