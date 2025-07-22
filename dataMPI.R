library(doMPI)
library(doRNG)
library(Matrix)
library(FNN)

requirements <- c("SpatialBSS","JADE","spGARCH","spdep","sp","dplyr","sf","moments", "Matrix")
#in puhti
#.libPaths(c("/projappl/project_2012081/project_rpackages_440", .libPaths()))

source("~/Desktop/Research/spBSS/datagen.R")

message("[", Sys.time(), "] Starting MPI cluster…")
cl <- startMPIcluster()
registerDoMPI(cl)
registerDoRNG(seed=123)
message("[", Sys.time(), "] Cluster started and RNG registered.")

ds <- rev(c(5, 10))
#ds <- rev(c(5, 10, 20, 40 60, 80, 100, 120, 140, 160))
n_reps <- 100

for (d in ds) {
  message("[", Sys.time(), "] Beginning loop for d = ", d)
  
  field <- gen_field(d, seed = 123 + d)
  W <- gen_knn(field)

  message("[", Sys.time(), "] Starting foreach for ", n_reps, " reps…")
  results_d <- foreach(i = 1:n_reps,
                       .combine  = "c",
                       .packages = requirements) %dorng% {
       gen_sources(d, W)
  }
  message("[", Sys.time(), "] Saving results to data_", d, ".rds …")
  
  saveRDS(results_d, file = sprintf("data_%d.rds", d), compress = "gzip")
  
  message("[", Sys.time(), "] Saved results.")
  
  message("[", Sys.time(), "] Clearing memory…")
  
  rm(results_d); gc()
  
  message("[", Sys.time(), "] Cleanup done.")
  
}

message("[", Sys.time(), "] About to closeCluster(cl)…")
closeCluster(cl)
message("[", Sys.time(), "] closeCluster() returned.")

message("[", Sys.time(), "] About to call mpi.finalize()…")
mpi.quit()
message("[", Sys.time(), "] mpi.finalize() returned. Script should exit now.")
