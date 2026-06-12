# -----------------------------------------------------------------------------
# Spatial BSS MPI Data Analysis Framework
# -----------------------------------------------------------------------------
#
# In Slurm choose the setting to run and values of the parameter d. Note that
# the total number of data points is d^2. Add the following after module loads:
#
# SETTING="setting_1/2/3/4"
# DS_VALUES="40,80,120"
# ensure all matrices are bounded. This natively corrects edge effects and forces the optimizer to equally weight th
# Then add to the srun of this R script the following:
#
# --setting=${SETTING} --ds=${DS_VALUES}
#
# See /scripts/example_mpi.sh for and example. This script generates the data
# for a given setting and values of d. The data is saved in 
# /data/setting_X/data_d.rds for setting number "X" and d value "d".
#
# Requirements:
#   - Relies on helper data generation functions sourced from "R/datagen.R" and
#     analysis functions sources from "R/analysis.R".
#
# -----------------------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)

parse_arg <- function(name, default = NULL) {
  val <- grep(paste0("^--", name, "="), args, value = TRUE)
  if (length(val) == 0) return(default)
  sub(paste0("^--", name, "="), "", val)
}

setting <- parse_arg("setting", default = "setting_3")
ds_input <- parse_arg("ds", default = "5")
ds <- as.integer(strsplit(ds_input, ",")[[1]])
ds <- c(10, 20, 40)

clear_raw_results <- function(settings = 1:4) {
  cat("\n--- Commencing Pre-Run Cleanup ---\n")
  
  for (s in settings) {
    # Define the directory path
    #dir_path <- sprintf("results/raw/setting_%d", s)
    dir_path <- sprintf("results/raw/setting_%d", s)
    
    if (dir.exists(dir_path)) {
      files_to_delete <- list.files(path = dir_path, pattern = "\\.rds$", full.names = TRUE)
      
      if (length(files_to_delete) > 0) {
        file.remove(files_to_delete)
        cat(sprintf("Cleared %d obsolete files from %s\n", length(files_to_delete), dir_path))
      } else {
        cat(sprintf("Directory %s is already clean.\n", dir_path))
      }
    } else {
      dir.create(dir_path, recursive = TRUE)
      cat(sprintf("Created new directory: %s\n", dir_path))
    }
  }
  cat("--- Cleanup Complete. Ready for new simulation. ---\n\n")
}


#clear results form previous run?
clear_raw_results(settings = 3)

message("Setting: ", setting)
message("Values of parameter d: ", paste(ds, collapse = ", "))
master <- c("doParallel", "doRNG", "Matrix", "doMPI", "spdep", "dplyr", "SpatialBSS")
workers <- c("SpatialBSS","JADE","spGARCH","spdep","sp",
             "dplyr","sf","moments", "Matrix")
#in puhti
#.libPaths(c("/projappl/project_2012081/project_rpackages_440", .libPaths()))

#load packages for master script
for (pkg in master) {
  if (!require(pkg, character.only = TRUE)) {
    library(pkg, character.only = TRUE, quietly=TRUE)
  }
}

#inherit functions from analysis.R
#source("/scratch/project_2012081/spBSS/R/analysis.R")
source("~/spBSS/R/analysis.R")



message("[", Sys.time(), "] Starting MPI cluster...")

cl <- makeCluster(8)
registerDoParallel(cl)
registerDoRNG(seed=123)

message("[", Sys.time(), "] Cluster started and RNG registered.")

n_reps <- 2000

for (d in ds) {
  
  all_res <- list()
  
  message("[", Sys.time(), "] Starting loop for d = ", d, ".")
  
  # in Puhti supercomputer
  # filename <- sprintf("/scratch/project_2012081/spBSS/data/%s/data_%d.rds", setting, d)
  
  # locally
  filename <- sprintf("~/spBSS/data/%s/data_%d.rds", setting, d)
  
  data <- readRDS(filename)
  field <- gen_field(d)
  
  #baseline
  bd <- c(0, 1, 1, 2, 2, 3)
  
  rings_all <- gen_rings(field, bd, row_standardize = FALSE, symmetrize = TRUE, f0 = TRUE)
  
  kernels_sparse_f0 <- rings_all$kernels_sparse
  kernels_f0 <- rings_all$kernels
  
  kernels_sparse <- kernels_sparse_f0[-1]
  kernels <- kernels_f0[-1]
  
  message("[", Sys.time(), "] Starting dorng...")
  
  df <- foreach(r = 1:n_reps, .combine = rbind, .packages = workers) %dopar% {
    sources <- data[[r]]
    p <- ncol(sources)
    
    A <- matrix(rnorm(p*p), p, p)
    X <- tcrossprod(sources, A)
    
    spFOBI_md <- spFOBI(field, X, A, kernels_sparse)$md
    spJADE_md <- spJADE(field, X, A, kernels, kernels_sparse)$md
    
    spFOBI_md_f0 <- spFOBI(field, X, A, kernels_sparse_f0)$md
    spJADE_md_f0 <- spJADE(field, X, A, kernels_f0, kernels_sparse_f0)$md
    
    res <- lcovbss_fobi_jade(X, A, kernels) 
    
    #compile Results
    data.frame(
      spFOBI    = spFOBI_md,
      spJADE    = spJADE_md,
      spFOBI_f0 = spFOBI_md_f0,
      spJADE_f0 = spJADE_md_f0,
      SBSS      = res$md_sbss,
      FOBI      = res$md_fobi,
      JADE      = res$md_jade,
      stringsAsFactors = FALSE
    )
  }
  
  all_res[[as.character(d)]] <- df
  
  message("[", Sys.time(), "] Saving results to /results/raw...")
  
  output_file <- sprintf("~/spBSS/results_wo_RO/raw/%s/res_%d.rds", setting, d)
  saveRDS(do.call(rbind, all_res), output_file)
  
  message("[", Sys.time(), "] Saved results.")
  message("[", Sys.time(), "] Finished d = ", d)
  
  rm(df, all_res, rings_all, kernels, kernels_sparse, kernels_f0, kernels_sparse_f0)
  gc()
}

message("[", Sys.time(), "] Closing cluster...")

stopCluster(cl)

message("[", Sys.time(), "] Closed sucesfully.")


