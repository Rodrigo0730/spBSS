# this file creates the irregular datasets for the analysis which are not included
# in the final paper, master thesis? 

## ADD TO GITIGNORE ##

library(spdep)
library(Matrix)

source("R/datagen.R")

irregular <- function(s, d) {
  full <- readRDS(sprintf("~/Desktop/Research/spBSS/data/setting_%s/data_100.rds", s))
  coords_full <- gen_field(100) 
  set.seed(123)
  idx <- sample(nrow(coords_full), d^2, replace = FALSE)
  coords <- coords_full[idx, ]
  data <- lapply(full, function(mat) mat[idx, , drop = FALSE])
  filename <- sprintf("~/Desktop/Research/spBSS/data/irregular_%s/data_%d.rds", s, d)
  saveRDS(
    list(
      data   = data,    
      coords = coords  
    ),
    filename,
    compress = "gzip"
  )
  
}

ds <- c(5, 10, 20, 40, 60, 80)
for (d in ds) {
  irregular(1, d)
}



