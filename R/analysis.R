
# -----------------------------------------------------------------------------
# Spatial BSS Helper Functions for MPI Simulation Framework
# -----------------------------------------------------------------------------
#
# Implemented Methods:
#   - spFOBI: Fobi-based Spatial Independent Component Estimation
#   - spJADE: Jade-based Spatial Independent Component Estimation
#   - SBSS: Spatial BSS based on local covariance matrices.
#   - FOBI: Fourth Order Blind Identification.
#   - JADE: Joint Approximate Diagonalization of Eigen-matrices.
#
#
# Requirements:
#   - Relies on helper data generation functions sourced from "datagen.R".
#
# -----------------------------------------------------------------------------

source("R/datagen.R")

whiten_data <- function(X, na.action = na.fail) {
  
  # Whitening exactly as in from JADE:::FOBI to
  # ensure all 5 implemented methods undergo the same whitening
  
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

compute_bkl <- function(k, l, kmat, Y) {
  
  # For n distinct locations, computes the fourth order cross-moment
  # matrices B_f^{kl}. The computations are vectorized and kmat is
  # a matrix such that the ij-th element corresponds to f(u_i - u_j).
  
  n <- nrow(Y)
  vkl <- Y[, k] * Y[, l]
  w <- kmat %*% vkl
  return(as.matrix(crossprod(Y * w, Y)/n))
}

compute_ckl <- function(k, l, kmat, lcov, Y) {
  
  # For n distinct locations, computes the local fourth order
  # cross-cumulant matrices C_f^{kl} of the whitened random field Y.
  # The computations are vectorized, Lcov is a local covariance matrix
  # and kmat is a matrix such that the ij-th element corresponds to 
  # f(u_i - u_j).
  
  n <- nrow(Y)
  p <- ncol(Y)
  Bkl <- compute_bkl(k, l, kmat, Y)
  
  u <- lcov[, k, drop=FALSE]
  v <- lcov[, l, drop=FALSE]
  L <- if (k==l) 2*tcrossprod(u) else tcrossprod(u, v) + tcrossprod(v, u)
  
  Ckl <- if (k==l) Bkl - L - diag(p) else Bkl - L
  
  return(as.matrix(Ckl))
}

spFOBI <- function(field, sources, kernels_sparse, eps= 1e-06, maxiter=100, seed=NULL) {
  
  # Main function to apply spFOBI to simulated sources. The function first creates
  # a mixing matrix A and computes the mixed sources as X = sources * A^T. Then, X is
  # whitened using whiten_data() and proceeds to compute the B_f matrices for joint
  # diagonalization utilizing JADE::frjd. It returns the array of matrices which underwent
  # joint diagonalization and the Minimum Distance Index value from JADE::MD of the
  # estimated unmixing matrix and the mixing matrix A.
  #
  # Here the kernels matrices must be sparse, generated, for example, using gen_rings()
  # from the datagen.R helper file.
  #
  # Possible future implementation: compare performance with weighted diagonalization.
  
  if (!is.numeric(field))
    stop("invalid field")
  
  K <- length(kernels_sparse)
  n <- nrow(sources)
  p <- ncol(sources)
  
  set.seed(seed)
  A <- matrix(rnorm(p*p), p, p)
  X <- tcrossprod(sources, A)

  whitening <- whiten_data(X)
  Y <- whitening$Y 
  W0 <- whitening$W0 
  
  B_array <- array(0, dim = c(p, p, K))
  
  for (i in 1:K) {
    kmat <- kernels_sparse[[i]]
    B_sum <- matrix(0, p, p)
    for (k in 1:p) {
      B <- compute_bkl(k, k, kmat, Y)
      B_sum <- B_sum + B
    }
    
    B_array[, , i] <- B_sum
  }
  
  # Weighted Diagonalization could be implemented as follows
  # w <- c(1, 1, 1, 1)
  # B_array <- sweep(B_array, 3, w, `*`)
  
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

spJADE <- function(field, sources, kernels, kernels_sparse, eps= 1e-06, maxiter=100, seed=NULL) {
  
  # Main function to apply spJADE to simulated sources. The function first creates
  # a mixing matrix A and computes the mixed sources as X = sources * A^T. Then, X is
  # whitened using whiten_data() and computes the Local Covariance matrices using
  # SpatialBSS::local_covariance_matrix. Moreover, computes the C_f^{kl} matrices for 
  # joint diagonalization utilizing JADE::frjd. It returns the array of matrices which 
  # underwent joint diagonalization and the Minimum Distance Index value from JADE::MD 
  # of the estimated unmixing matrix and the mixing matrix A.
  #
  # Here the kernels_sparse matrices must be sparse, generated, for example, using gen_rings()
  # from the datagen.R helper file. But it is also needed to pass the full matrices as
  # SpatialBSS::local_covariance_matrix does not accept sparse matrices.
  #
  # Possible future implementations: compare performance with weighted diagonalization. 
  # Write own local_covariance_matrix function to take sparse matrices.
  
  if (!is.numeric(field))
    stop("invalid field")
  
  K <- length(kernels_sparse)
  n <- nrow(sources)
  p <- ncol(sources)
  
  set.seed(seed)
  A <- matrix(rnorm(p*p), p, p)
  X <- tcrossprod(sources, A)
  
  whitening <- whiten_data(X)
  Y <- whitening$Y 
  W0 <- whitening$W0 
  
  lcovs <- local_covariance_matrix(x = Y, kernel_list = kernels, lcov="lcov")

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
  
  # w <- c(1, 0.5)
  # C_array <- sweep(C_array, 3, w, `*`)
  
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
    W_est <- crossprod(jd_res$V, W0)

    md <- MD(W_est, A)
  }
  return(list(
    md = md,
    C_array = C_array
  ))
  
}

lcovbss_fobi_jade <- function(sources, kernels, seed=NULL) {
  
  # This function performs the same procedure as for spFOBI and spJADE but for 
  # sbss (SpatialBSS::sbss), classic FOBI and JADE (JADE::FOBI; JADE::JADE). Returns
  # the Minimum Distance Index values (JADE::MD) of the estimated unmixing matrix and 
  # mixing matrix A.
  
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
  

  md_sbss <- tryCatch(
    MD(W_est_sbss, A),
    error = function(e) NA
  )
  md_fobi <- if (!is.matrix(W_est_fobi)) NA else MD(W_est_fobi, A)
  md_jade <- if (!is.matrix(W_est_jade)) NA else MD(W_est_jade, A)
  
  return(list(
    md_sbss=md_sbss,
    md_fobi=md_fobi,
    md_jade=md_jade
  ))
  
}

#test function to evaluate the different methods
test <- function(ds, n_rep = 1, load=FALSE) {
  
  if (load==TRUE) {
    requirements <- c("doRNG","SpatialBSS","JADE","spGARCH","spdep","sp","dplyr","sf",
                      "moments", "Matrix")
    for (pkg in requirements) {
      if (!require(pkg, character.only = TRUE)) {
        library(pkg, character.only = TRUE)
      }
    }
  }
  

  
  all_res <- list()
  
  filename <- paste0("~/Desktop/Research/spBSS/data/setting_3/data_10.rds")
  file <- readRDS(filename)
  
  data <- file$data
  field <- file$coords
  
  bd <- c(0, 1, 1, 2, 2, 3)
  rings <- gen_rings(field, bd)
  
  kernels <- rings$kernels
  kernels_sparse <- rings$kernels_sparse
  
  for (r in 1:2000) {
  
    sources <- data[[r]]
  
    spFOBI_md <- spFOBI(field, sources, kernels_sparse)$md
    spJADE_md <- spJADE(field, sources, kernels, kernels_sparse)$md
    
    res <- lcovbss_fobi_jade(sources, kernels)
  }
  
  data.frame(
    spFOBI = spFOBI_md,
    spJADE = spJADE_md,
    sbss = res$md_sbss,
    fobi = res$md_fobi,
    jade = res$md_jade,
    stringsAsFactors = FALSE
  )
}





  


  
  



