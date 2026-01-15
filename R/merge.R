# this file is just a merge helper to create the data to be analyzed when
# sometimes some runs have had to be rerun in the supercomputer 
# to regenerate some data without having to generate the complete
# data, just some sources and not all of them

merge <- function (s1, s2, d) {
  
  original <- readRDS(sprintf("~/Desktop/Research/spBSS/data/%s/data_%d.rds", s1, d))
  fix <- readRDS(sprintf("~/Desktop/Research/spBSS/data/%s/data_%d.rds", s2, d))
  
  for (i in (1:2000)) {
    fix[[i]] <- cbind(fix[[i]][, 1, drop=FALSE], original[[i]][, 2], fix[[i]][, 2:3, drop=FALSE])
  }
  saveRDS(fix, sprintf("~/Desktop/Research/spBSS/data/%s/data_%d.rds", s2, d), compress = "gzip")
}

ds <- c(5, 10, 20, 40, 60, 80, 100)

for (d in ds) {
  merge("setting_3", "setting_4", d)
}


