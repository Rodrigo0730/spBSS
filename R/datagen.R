

# -------------------------------------------------------------------------
# Helper functions for data generation MPI framework
#
# This script provides a collection of helper functions for spatial grid
# generation, neighborhood construction, spatial weights matrix handling,
# and source generation (log-spARCH and log-spGARCH processes). These 
# functions are designed to support distributed simulations or computations 
# in an MPI environment.
#
# Functions included:
#   - gen_field: Creates a grid of coordinates for a d x d field.
#   - gen_rings: Builds ring-shaped spatial kernel matrices (dense and sparse).
#   - gen_nb: Generates rook/queen adjacency-based neighbors.
#   - nb2W: Converts neighbors lists to sparse weights matrices.
#   - nth_W: Constructs k-th order spatial weight matrices.
#   - gen_sources_setting_3/4: Simulates spARCH/spGARCH processes with 
#     specified spatial dependencies.
# 
# To-do:
# Write functions to generate sources for settings 1 and 2 once the results
# of settings 3 and 4 are done. Might be useful to incorporate neighbor 
# based rings instead of distanced based like gen_rings() does.

gen_field <- function(d) {
  
  # Generates a d x d grid of coordinates (field) as a matrix with
  # columns x and y. Returns a matrix with d^2 rows, where each row
  # represents a grid point.

  xs <- rep(1:d, each = d)
  ys <- rep(1:d, times = d)
  coords <- cbind(x = xs, y = ys)
  return(coords)
}

gen_rings <- function(field, bandwidths) {
  
  # Generates a series of ring-shaped spatial kernel matrices (both dense and 
  # sparse formats) for a given field and a set of bandwidths. The function 
  # computes ring kernels using spatial_kernel_matrix, converts them to sparse 
  # dgCMatrix format, and prepends an identity matrix to represent the kernel f0  
  # (no spatial information). Returns a list containing both the sparse 
  # and dense versions of the kernels.
  
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

gen_nb <- function(d, type = c("rook","queen")) {
  
  # Generates a neighbors list (nb) for a d x d grid based on the specified 
  # neighborhood definition ("rook" or "queen").
  
  type <- match.arg(type)
  nb <- cell2nb(nrow = d, ncol = d, type = type, torus = FALSE)
  return(nb)
}

nb2W <- function(nb, style = "W") {
  
  # Converts a neighbors list (nb) into a sparse spatial weights
  # matrix of class "dgCMatrix". The style parameter defines how the weights 
  # are standardized (e.g., "W" for row-standardized weights).

  lw  <- nb2listw(nb, style = style)
  W   <- as(listw2mat(lw), "dgCMatrix")
  return(W)
}

nth_W <- function(nb, k = 2, style = "W") {
  
  # Creates a spatial weights matrix of the k-th order neighbors from a given 
  # neighbors list (nb). It computes neighbors up to lag k, extracts the 
  # k-th order neighbor list, and converts it into a spatial weights matrix 
  # using the specified style (default "W" for row-standardized weights).
  
  lagged_l <- suppressWarnings(nblag(nb, maxlag = k))
  nb_k <- lagged_l[[k]]
  Wk <- nb2W(nb_k, style = style)
  return(Wk)
}
gen_sources_setting_1 <- function() {
  # Check notes for specification of the settings and required libraries.
}

gen_sources_setting_2 <- function() {
  # Same as above but with student-t distributed innovations.
}

gen_sources_setting_3 <- function(d, W1, W2, W3, seed = NULL) {
  
  # helper function to generate the sources for the third setting
  # the sources are as following:
    # 1: log-spARCH(1) with rho = 0.9 and alpha = 0.5
    # 2: log-spARCH(1) with rho =0.7 and alpha = 0.5
    # 3: log-spGARCH(1, 2) with rho = 0.4, lambda = 0.9 and alpha = 1
    # 4: log-spGARCH(3, 1) with rho = 0.8, lambda = 0.5 and alpha = 1
    # 5: log-spGARCH(2, 3) with rho=0.4, lambda = 0.3 and alpha = 0.5
  # all with Gaussian innovations. Moreover, rho and lambda are the
  # spatial dependence parameters for W1 (ARCH component) and
  # W2 (GARCH component) respectively, and alpha represents the 
  # baseline volatility or unconditional variance.
  
  if (!is.null(seed)) set.seed(seed)
  
  params <- list(
    list(rho = 0.9, lambda = 0.0, alpha = 0.5, W1 = W1, W2 = W1),
    list(rho = 0.7, lambda = 0.0, alpha = 0.5, W1 = W2, W2 = W1),
    list(rho = 0.4, lambda = 0.9, alpha = 1, W1 = W1, W2 = W2),
    list(rho = 0.8, lambda = 0.5, alpha = 1, W1 = W3, W2 = W1),
    list(rho = 0.4, lambda = 0.3, alpha = 0.5, W1 = W2, W2 = W3)
  )
  
  sources_list <- lapply(seq_along(params), function(i) {
    par <- params[[i]]
    if (i <= 2L) {
      sim.spARCH(
        n       = d^2,
        rho     = par$rho,
        alpha   = par$alpha,
        W       = par$W1,
        type    = "log-spARCH",
        control = list(silent=TRUE)
      )
    } else {
      sim.spGARCH(
        n       = d^2,
        rho     = par$rho,
        alpha   = par$alpha,
        lambda  = par$lambda,
        W1      = par$W1,
        W2      = par$W2,
        type    = "log-spGARCH",
      )
    }
  })
  
  return(sources_list)
}

gen_sources_setting_4 <- function(d, W1, W2, W3, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  
  # This function uses a modified version library spGARCH >2.3.0.
  
  # helper function to generate the sources for the fourth setting
  # same sources as setting 3 but in sources 1, 3, 4, 5 the innovations
  # now follow a student-t distribution. The second source still has
  # Gaussian innovations.
  
  # to-do: make sure the student-t distributions are applied and 
  # modify the code in the spGARCH package to handle the change in
  # innovations properly.
  
  params <- list(
    list(rho = 0.9, lambda = 0.0, alpha = 0.5, W1 = W1, W2 = W1),
    list(rho = 0.7, lambda = 0.0, alpha = 0.5, W1 = W1, W2 = W1),
    list(rho = 0.4, lambda = 0.9, alpha = 1, W1 = W1, W2 = W2),
    list(rho = 0.8, lambda = 0.5, alpha = 1, W1 = W3, W2 = W1),
    list(rho = 0.4, lambda = 0.3, alpha = 0.5, W1 = W2, W2 = W3)
  )
  
  sources_list <- lapply(seq_along(params), function(i) {
    par <- params[[i]]
    if (i <= 2L) {
      sim.spARCH(
        n       = d^2,
        rho     = par$rho,
        alpha   = par$alpha,
        W       = par$W1,
        type    = "log-spARCH",
        control = list(silent=TRUE)
      )
    } else {
      sim.spGARCH(
        n       = d^2,
        rho     = par$rho,
        alpha   = par$alpha,
        lambda  = par$lambda,
        W1      = par$W1,
        W2      = par$W2,
        type    = "log-spGARCH",
      )
    }
  })
  
  return(sources_list)
}























