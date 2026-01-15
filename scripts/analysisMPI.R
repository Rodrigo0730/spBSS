
# -----------------------------------------------------------------------------
# Spatial BSS MPI Data Analysis Framework
# -----------------------------------------------------------------------------
#
# In Slurm choose the setting to run and values of the parameter d. Note that
# the total number of data points is d^2. Add the following after module loads:
#
# SETTING="setting_1/2/3/4"
# DS_VALUES="40,80,120"
#
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

setting <- parse_arg("setting", default = "irregular_1")
ds_input <- parse_arg("ds", default = "5")
ds <- as.integer(strsplit(ds_input, ",")[[1]])

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
source("~/Desktop/Research/spBSS/R/analysis.R")



message("[", Sys.time(), "] Starting MPI cluster...")

cl <- makeCluster(6)
registerDoParallel(cl)
registerDoRNG(seed=123)

message("[", Sys.time(), "] Cluster started and RNG registered.")

n_reps <- 2000
ds <- c(5, 10, 20, 40, 60, 80)

for (d in ds) {
  
  all_res <- list()
  
  message("[", Sys.time(), "] Starting loop for d = ", d, ".")
  
  # in Puhti supercomputer
  # filename <- sprintf("/scratch/project_2012081/spBSS/data/%s/data_%d.rds", setting, d)
  
  # locally
  filename <- sprintf("~/Desktop/Research/spBSS/data/%s/data_%d.rds", setting, d)
  
  data <- readRDS(filename)
  
  field <- gen_field(d)
  # 
  # field <- data$coords
  # data <- data$data

  bd <- c(0, 1, 1, 2, 2, 3)
  rings <- gen_rings(field, bd)
  kernels_sparse <- rings$kernels_sparse
  kernels <- rings$kernels
  
  message("[", Sys.time(), "] Starting dorng...")
  
  df <- foreach(r = 1:n_reps, .combine = rbind, .packages = workers) %dopar% {
    sources <- data[[r]]
    
    spFOBI_md <- spFOBI(field, sources, kernels_sparse)$md
    spJADE_md <- spJADE(field, sources, kernels, kernels_sparse)$md
    
    res <- lcovbss_fobi_jade(sources, kernels)
    data.frame(
      spFOBI = spFOBI_md,
      spJADE = spJADE_md,
      SBSS = res$md_sbss,
      FOBI = res$md_fobi,
      JADE = res$md_jade,
      stringsAsFactors = FALSE
    )
  }
  
  all_res[[as.character(d)]] <- df
  
  message("[", Sys.time(), "] Saving results to /results/raw...")
  
  output_file <- sprintf("~/Desktop/Research/spBSS/results/raw/%s/res_%d.rds", setting, d)
  saveRDS(do.call(rbind, all_res), output_file)
  
  message("[", Sys.time(), "] Saved results.")
  
  message("[", Sys.time(), "] Finished d = ", d)
  
  rm(df, all_res)
  gc()
}

message("[", Sys.time(), "] Closing cluster...")

stopCluster(cl)

message("[", Sys.time(), "] Closed sucesfully.")


