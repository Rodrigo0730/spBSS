library(doMPI)
library(doRNG)
library(Matrix)
library(FNN)
library(dplyr)

requirements <- c("SpatialBSS",
                  "JADE",
                  "spGARCH",
                  "spdep",
                  "sp",
                  "dplyr",
                  "sf",
                  "moments",
                  "Matrix")
#in puhti
#.libPaths(c("/projappl/project_2012081/project_rpackages_440", .libPaths()))

#inherit functions from analysis.R
source("~/Desktop/Research/spBSS/analysis.R")

#in puhti
#source("analysis.R")

message("[", Sys.time(), "] Starting MPI cluster...")

cl <- startMPIcluster()
registerDoMPI(cl)
registerDoRNG(seed=123)

message("[", Sys.time(), "] Cluster started and RNG registered.")

ds <- c(10)
#ds <- rev(c(5, 10, 20, 40 60, 80, 100, 120, 140, 160))

all_res <- list()
for (d in ds) {
  
  message("[", Sys.time(), "] Starting loop for d = ", d, ".")
  
  filename <- file.path(paste0("data_", d, ".rds"))
  data <- readRDS(filename)
  
  field <- gen_field(d, seed = 123 + d)
  bd <- c(0,1,1,2,2,3,3,4,4,5,5,6,6,7,7,8)
  rings <- gen_rings(field, bd)
  kernels_sparse <- rings$kernels_sparse
  kernels <- rings$kernels
  
  message("[", Sys.time(), "] Starting dorng...")
  
  results_d <- foreach(
    r = 1:1000,
    .packages = requirements,
    .combine  = rbind
  ) %dorng% {
    sources <- data[[r]]
    
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
  
  # Reshape to long format
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
  
  # Store
  all_results[[as.character(d)]] <- df_long
  
  message("[", Sys.time(), "] Finished d = ", d)
  
  # Clean up
  rm(results_d, df_long)
  gc()
}

message("[", Sys.time(), "] Closing cluster...")

closeCluster(cl)
mpi.quit()

message("[", Sys.time(), "] Closed sucesfully.")

res_df <- do.call(rbind, all_res)

means <- res_df %>% group_by(d, model) %>% summarise(mean_md = mean(md, na.rm = TRUE))

print(means)


