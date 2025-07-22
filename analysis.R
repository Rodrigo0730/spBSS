requirements <- c("SpatialBSS","tsBSS", "JADE","spGARCH","spdep","sp","dplyr","sf","moments", "BSSprep", "Matrix", "FNN", "data.table")
for (pkg in requirements) {
  if (!require(pkg, character.only = TRUE)) {
    library(pkg, character.only = TRUE)
  }
}
source("datagen.R")
# need kernels and field, same field reproducible with seed = 123 + d


#whitening exactly as in from JADE:::FOBI
#ensure all 5 methods undergo the same whitening
whiten_data <- function(X, na.action = na.fail) {
  X <- na.action(X)
  X <- as.matrix(X)
  
  X_centered <- sweep(X, 2, colMeans(X), "-")
  
  COV <- crossprod(X_centered) / nrow(X_centered)
  
  EVD <- eigen(COV, symmetric = TRUE)
  U <- EVD$vectors
  Lambda <- EVD$values
  
  COV_inv_sqrt <- U %*% diag(1 / sqrt(Lambda)) %*% t(U)
  
  X_whitened <- X_centered %*% COV_inv_sqrt
  
  return(list(
    Y = X_whitened,
    W0 = COV_inv_sqrt
  ))
}

#computes bkl matrix
#I have just thought that normalizing by sum(kmat) instead of n gives way better results,
#need to look at the mathematics behind it but it makes sense and should work, moreover
#consider doing the same in Lcov
compute_bkl <- function(k, l, kmat, Y) {
  n <- nrow(Y)
  vkl <- Y[, k] * Y[, l]
  w <- kmat %*% vkl
  #return(as.matrix(crossprod(Y * w, Y)/sum(kmat)))
  return(as.matrix(crossprod(Y * w, Y)/n))
}

#calculates Ckl matrices
compute_ckl <- function(k, l, kmat, lcov, Y) {
  n <- nrow(Y)
  p <- ncol(Y)
  Bkl <- compute_bkl(k, l, kmat, Y)
  
  Ekl <- matrix(0, p, p)
  if (k == l) {
    Ekl[k, l] <- 2
  } else {
    Ekl[k, l] <- 1
    Ekl[l, k] <- 1
  }
  L <- lcov %*% Ekl %*% t(lcov)
  
  trEkl <- if (k==l) 1 else 0
  
  Sf <- sum(kmat) /n
  
  Ckl <- Bkl - L  - trEkl * Sf * diag(p)
  
  return(as.matrix(Ckl))
}

#fspice method for the simulation
fspice <- function(field, sources, kernels_sparse, eps= 1e-06, maxiter=100, seed=NULL) {
  if (!is.numeric(field))
    stop("invalid field")
  
  K <- length(kernels_sparse)
  n <- nrow(sources)
  p <- ncol(sources)
  
  set.seed(seed)
  A <- matrix(rnorm(p*p), p, p)
  X <- tcrossprod(sources, A)

  whitening <- whiten_data(X)
  Y <- whitening$Y #whitened data
  W0 <- whitening$W0 #whitening matrix

  B_array <- array(0, dim = c(p, p, K))
  #k is small (number of kernels), iterate through kernels
  for (i in 1:K) {
    kmat <- kernels_sparse[[i]]
    B_sum <- matrix(0, p, p)
    for (k in 1:p) {
      B <- compute_bkl(k, k, kmat, Y)
      B_sum <- B_sum + B
    }
    
    B_array[, , i] <- B_sum
  }
  
  jd_res <- tryCatch(
    suppressWarnings(frjd(B_array, maxiter = maxiter, eps = eps)),
    error = function(e) {
      message("frjd failed: ", e$message)
      return(NULL)
    }
  )
  
  if (is.null(jd_res)) {
    md <- NA
  } else {
    W_est <- crossprod(jd_res$V, W0)
    md <- MD(W_est, A)
  }
  return(list(
    md = md,
    B_array = B_array
  ))
  
}

#jspice method for the simulations, takes original kernels (not sparse) to compute the lcov matrices
jspice <- function(field, sources, kernels, kernels_sparse, eps= 1e-06, maxiter=100, seed=NULL) {
  if (!is.numeric(field))
    stop("invalid field")
  
  K <- length(kernels_sparse)
  n <- nrow(sources)
  p <- ncol(sources)
  
  set.seed(seed)
  A <- matrix(rnorm(p*p), p, p)
  X <- tcrossprod(sources, A)
  
  whitening <- whiten_data(X)
  Y <- whitening$Y #whitened data
  W0 <- whitening$W0 #whitening matrix
  
  lcovs <- local_covariance_matrix(x = Y, kernel_list = kernels)
  
  C_array <- array(0, dim = c(p,p,p*p*K)) 
  
  idx <- 1
  for (i in 1:K) {
    kmat <- kernels_sparse[[i]]
    lcov <- lcovs[[i]]
    for (k in 1:p) {
      for (l in 1:p) {
        C_array[ , , idx] <- compute_ckl(k, l, kmat, lcov, Y)
        idx <- idx + 1
      }
    }
    
  }
  
  jd_res <- tryCatch(
    suppressWarnings(frjd(C_array, maxiter = maxiter, eps = eps)),
    error = function(e) {
      message("frjd failed: ", e$message)
      return(NULL)
    }
  )
  
  if (is.null(jd_res)) {
    md <- NA
  } else {
    U <- jd_res$V
    W_est <- t(U) %*% W0

    md <- MD(W_est, A)
  }
  return(list(
    md = md,
    C_array = C_array
  ))
  
}

#function to compute the other three methods for the simulation
lcovbss_fobi_jade <- function(sources, kernels, seed=NULL) {
  
  K <- length(kernels)
  n <- nrow(sources)
  p <- ncol(sources)
  
  set.seed(seed)
  A <- matrix(rnorm(p*p), p, p)
  sources <- tcrossprod(sources, A)
  
  W_est_sbss <- tryCatch(
    sbss(x = sources, kernel_list = kernels, lcov = "lcov")$w,
    error = function(e) NA
  )

  W_est_fobi <- tryCatch(
    FOBI(X = sources)$W,
    error = function(e) NA
  )

  W_est_jade <- tryCatch(
    JADE(X = sources, n.comp = p)$W,
    error = function(e) NA
  )

  md_sbss <- if (!is.matrix(W_est_sbss)) NA else MD(W_est_sbss, A)
  md_fobi <- if (!is.matrix(W_est_fobi)) NA else MD(W_est_fobi, A)
  md_jade <- if (!is.matrix(W_est_jade)) NA else MD(W_est_jade, A)
  
  # md_sbss <- NA 
  # md_fobi <- NA 
  # md_jade <- NA
  
  return(list(
    md_sbss=md_sbss,
    md_fobi=md_fobi,
    md_jade=md_jade
  ))
  
}

#test function to evaluate the different methods
test <- function(ds, n_rep = 1000) {
  all_res <- list()
  
  for (d in ds) {
    filename <- paste0("data_", d, ".rds")
    data <- readRDS(filename)
    field <- gen_field(d, seed = 123 + d)

    bd <- c(0, 2.2, 2.2, 4.4)
    rings <- gen_rings(field, bd)
    
    kernels <- rings$kernels
    kernels_sparse <- rings$kernels_sparse
    
    # kernels <- list(diag(1, 1600, 1600)) #f_0
    # kernels_sparse <- list(Diagonal(1600, x=1))

    df <- data.frame(
      fspice = numeric(n_rep),
      jspice = numeric(n_rep),
      sbss = numeric(n_rep),
      fobi = numeric(n_rep),
      jade = numeric(n_rep),
      stringsAsFactors = FALSE
    )
    
    for (r in 1:n_rep) {
      sources <- data[[r]]
      df$fspice[r] <- fspice(field, sources, kernels_sparse, seed = r)$md
      df$jspice[r] <- jspice(field, sources, kernels, kernels_sparse, seed = r)$md
      
      res <- lcovbss_fobi_jade(sources, kernels, seed = r)
      df$sbss[r] <- res$md_sbss
      df$fobi[r] <- res$md_fobi
      df$jade[r] <- res$md_jade
    }

    all_res[[as.character(d)]] <- df
  }
  do.call(rbind, all_res)
}


#d <- 20
#filename <- paste0("data_", d, ".rds")
#data     <- readRDS(filename)
#field    <- gen_field(d, seed = 123 + d)
#bd       <- c(0, 0.25, 0.25, 0.5, 0.5, 0.75, 0.75, 1) * d
#bd <- c(0, 1) * sqrt(2)*d
#results <- data.frame()

#rings          <- gen_rings(field, bd)
#kernels        <- rings$kernels
# kernels_sparse <- rings$kernels_sparse
# 
# r <- 20
# sources <- data[[r]]
# res1 <- fspice(field, sources, kernels_sparse, seed = r)
# res2 <- jspice(field, sources, kernels, kernels_sparse, seed = r)
# 
# set.seed(r)
# p <- 5
# A <- matrix(rnorm(p*p), p, p)
# X <- tcrossprod(sources, A)
# 
# whitening <- whiten_data(X)
# Y <- whitening$Y #whitened data
# 
# res3  <- MD(gFOBI(Y)$W, A)
# res4 <- MD(gJADE(Y)$W , A)
# results <- rbind(results, data.frame(fmd = res1$md, jmd = res2$md, gfobi = res3, gjade=res4))
# 
# print(results)

res <- test(ds=c(40), n_rep=200)



  


  
  



