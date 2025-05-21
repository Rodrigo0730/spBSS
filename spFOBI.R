library(SpatialBSS)
library(JADE)
library(spGARCH)
library(spdep)
library(sp)
library(dplyr)
library(sf)
library(moments)

set.seed(1234)
p <- 3
d <- 100
n <- d^2
coords <- matrix(rnorm(n * 2), ncol = 2)
dist_matrix <- as.matrix(dist(coords))

inv_dist <- 1 / dist_matrix
diag(inv_dist) <- 0

W <- inv_dist / rowSums(inv_dist)

#generate data from log-spARCH model for heteroscedastic
rho1 <- 0.5
rho2 <- 0.9
alpha1 <- 1
alpha2 <- 0.5

source1 <- sim.spARCH(n = n, rho = rho1, alpha = alpha1, W = W, type = "log-spARCH", control = list(seed = 224893301))
source2 <- sim.spARCH(n = n, rho = rho2, alpha = alpha2, W = W, type = "log-spARCH", control = list(seed = 575639611))
source3 <- sim.spARCH(n = n, rho = rho1, alpha = alpha2, W = W, type = "log-spARCH", control = list(seed = 907977532))

sources <- cbind(source1, source2, source3)

#True mixing matrix
set.seed(1234)
A <- matrix(runif(p^2, 0, 1), p, p)

#generate the observed data, whiten the data
X <- sources %*% t(A)
X_centered <- scale(X, center = TRUE, scale=FALSE)
cov_X <- cov(X_centered)
E <- eigen(cov_X)
whitening_mat <- E$vectors %*% diag(1 / sqrt(E$values)) %*% t(E$vectors)
X_white <- X_centered %*% whitening_mat

B_f_list <- list()
kernel_bandwidths <- c(0.1, 0.2, 0.3)
kernels <- spatial_kernel_matrix(coords, kernel_type = "gauss", kernel_parameters = kernel_bandwidths)

for (k in 1:p) {
  kernel_mat <- kernels[[k]]
  x_stk <- X_white[, k]
  x_squared <- x_stk^2
  
  mixed_moments_matrix <- tcrossprod(x_squared)
  
  weighted_sum <- sum(kernel_mat * (mixed_moments_matrix + (p - 1)))
  
  B_f <- matrix(0, p, p)
  B_f[k, k] <- weighted_sum / n
  
  B_f_list[[k + 1]] <- B_f
}

#jade
K <- length(B_f_list)
B_array <- array(0, dim = c(p, p, K))

for (i in 1:K) {
  B_array[, , i] <- B_f_list[[i]]
}

jd_res <- frjd(B_array)
U <- jd_res$V

W_est <- t(U) %*% whitening_mat

# sources
S_est <- X %*% t(W_est)

md <- MD(W_est, A)
  
  
  
