
# -----------------------------------------------------------------------------
# Spatial BSS MPI Data Analysis Framework
# -----------------------------------------------------------------------------
#
# In slurm choose the setting to run and values of the parameter d. Note that
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

setting <- parse_arg("setting", default = "setting_3")
ds_input <- parse_arg("ds", default = "20")
ds <- as.integer(strsplit(ds_input, ",")[[1]])

message("Setting: ", setting)
message("Values of parameter d: ", paste(ds, collapse = ", "))
master <- c("doRNG", "Matrix", "doMPI", "spdep", "dplyr")
workers <- c("SpatialBSS","JADE","spGARCH","spdep","sp",
             "dplyr","sf","moments", "Matrix")
#in puhti
.libPaths(c("/projappl/project_2012081/project_rpackages_440", .libPaths()))

#load packages for master script
for (pkg in master) {
  if (!require(pkg, character.only = TRUE)) {
    library(pkg, character.only = TRUE, quietly=TRUE)
  }
}

#inherit functions from analysis.R
source("/scratch/project_2012081/spBSS/analysis.R")


message("[", Sys.time(), "] Starting MPI cluster...")

cl <- startMPIcluster()
registerDoMPI(cl)
registerDoRNG(seed=123)

message("[", Sys.time(), "] Cluster started and RNG registered.")

n_reps <- 2000

for (d in ds) {
  
  message("[", Sys.time(), "] Starting loop for d = ", d, ".")
  
  filename <- sprintf("/scratch/project_2012081/spBSS/data/%s/data_%d.rds", setting, d)
  data <- readRDS(filename)
  
  field <- gen_field(d)
  bd <- c(0, 1, 1, 2, 2, 3)
  rings <- gen_rings(field, bd)
  kernels_sparse <- rings$kernels_sparse
  kernels <- rings$kernels
  
  message("[", Sys.time(), "] Starting dorng...")
  
  results_d <- foreach(
    r = 1:n_reps,
    .packages = requirements,
    .combine  = rbind
  ) %dorng% {
    sources <- data[[r]]
    
    #check seed
    
    m_fspice <- fspice(field, sources, kernels_sparse, seed = r)$md
    m_jspice <- jspice(field, sources, kernels, kernels_sparse, seed = r)$md
    tmp      <- lcovbss_fobi_jade(sources, kernels, seed = r)
    m_sbss   <- tmp$md_sbss
    m_fobi   <- tmp$md_fobi
    m_jade   <- tmp$md_jade
    
    data.frame(
      rep    = r,
      fspice = m_fspice,
      jspice = m_jspice,
      sbss   = m_sbss,
      fobi   = m_fobi,
      jade   = m_jade,
      stringsAsFactors = FALSE
    )
  }
  
  df_long <- reshape(
    results_d,
    varying   = list(names(results_d)[2:6]),
    v.names   = "md",
    timevar   = "model",
    times     = names(results_d)[2:6],
    direction = "long"
  )
  df_long$d <- d
  rownames(df_long) <- NULL
  
  output_file <- sprintf("/scratch/project_2012081/spBSS/results/raw/%s/res_%d.rds", setting, d)
  saveRDS(df_long, output_file)
  
  message("[", Sys.time(), "] Finished d = ", d)
  
  rm(results_d, df_long)
  gc()
}

message("[", Sys.time(), "] Closing cluster...")

closeCluster(cl)
mpi.quit()

message("[", Sys.time(), "] Closed sucesfully.")


