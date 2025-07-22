#auxiliary file with core functions used in the data generation process
#



requirements <- c("SpatialBSS","JADE","spGARCH","spdep","sp","dplyr","sf","moments", "Matrix", "FNN")
for (pkg in requirements) {
 if (!require(pkg, character.only = TRUE)) {
   library(pkg, character.only = TRUE)
 }
}

#generates the data points for given sample size n = d^2 and given seed,
#the data points are uniformly distributed in the square of area n

gen_field <- function(d, seed=NULL) {
  if (!is.null(seed)) set.seed(seed)
  n <- d^2
  coords <- matrix(runif(n*2, 0, d), ncol=2)
  return(coords)
}

#pairwise bandwidhts for the kernels from ring, that is: c(0, 0.1, 0.1, 0.5)
#for example

gen_rings <- function(field, bandwidths) {
  n <- nrow(field)
  kernels <- spatial_kernel_matrix(field, "ring", bandwidths)
  kernels_sparse <- lapply(kernels, function(W) {
    W_sp <- as(W, "dgCMatrix")
    W_sp
  })
  kernels <- c(list(diag(1, n, n)), kernels)
  kernels_sparse <- c(list(Diagonal(n, x=1)), kernels_sparse)
  return(list(kernels_sparse=kernels_sparse, kernels=kernels))
}

#generate a symmetric spatial weight matrix based on KNN
#fixed  to k = 15
#estimates sigma as median of the distances to the k nearest neighbors for 
#robustness
gen_weights <- function(field, k = 15) {
  n <- nrow(field)
  knn <- get.knn(field, k = k)
  
  i <- rep(1:n, each = k)
  j <- as.integer(knn$nn.index)
  d <- as.numeric(knn$nn.dist)
  sigma <- median(d)
  
  w <- exp(-(d^2) / (2 * sigma^2))
  
  W0 <- sparseMatrix(i = i, j = j, x = w, dims = c(n, n))
  W_sym <- (W0 + t(W0)) / 2
  rsums <- rowSums(W_sym)
  inv_rsums <- 1 / rsums
  inv_rsums[!is.finite(inv_rsums)] <- 0
  W <- Diagonal(x = inv_rsums) %*% W_sym
  W <- drop0(W, tol = 1e-8)
  

  return(as(W, "generalMatrix"))
}

#produce higher order matrices
make_weighted_sum <- function(W, alphas) {
  n <- nrow(W)
  W1 <- sparseMatrix(i = integer(0), j = integer(0), x = numeric(0), dims = c(n, n))
  W_power <- W
  for (i in seq_along(alphas)) {
    W1 <- W1 + alphas[i] * W_power
    if (i < length(alphas)) {
      W_power <- drop0(W_power %*% W, tol = 1e-8)
    }
  }
  return(as(W1, "generalMatrix"))
}

#generate sources for the datagen process
gen_sources <- function(d, W0, W3, W4, seed=NULL) {
  if (!is.null(seed)) set.seed(seed)
  
  params <- list(
    list(rho=0.1, lambda=0.7, alpha=0.5, W1=W0, W2=W0),
    list(rho=0.9, lambda=0.0, alpha=1, W1=W0, W2=W0),
    list(rho=0.1, lambda=0.8, alpha=0.5, W1=W0, W2=W0),
    list(rho=0.1, lambda=0.5, alpha=1, W1=W3, W2=W0),
    list(rho=0.5, lambda=0.0, alpha=0.5, W1=W3, W2=W0)
  )
  
  sources_list <- lapply(params, function(par) {
    sim.spGARCH(n = d^2, rho = par$rho, lambda = par$lambda, alpha = par$alpha, W1 = par$W1, W2 = par$W2, type = "log-spGARCH")
  })  
  
  sources <- do.call(cbind, sources_list)
  return(sources)
}

#test function
test <- function(d = 20, i) {
  field <- gen_field(d, seed = 123 + d)
  W0 <- gen_weights(field)
  
  alphas3 <- c(0.5714, 0.2857, 0.1429)
  alphas4 <- c(0.5333, 0.2667, 0.1333, 0.0667)
  
  W3  <- make_weighted_sum(W0, alphas3)
  W4  <- make_weighted_sum(W0, alphas4)
  
  sources <- gen_sources(d, W0 = W0, W3 = W3, W4 = W4, seed = i + 123)
  
  rm(W3, W4, field); gc()
  
  return(sources)
}





