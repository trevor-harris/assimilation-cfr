rm(list = ls()); gc()

# get args
args = commandArgs(TRUE)
cov_scale = as.double(args[1])
simulations = as.integer(args[2])

library(extdepth)
library(dplyr)
library(reshape2)
library(plgp)

source('../../ks_field_functions.R')
source('../../sim_functions.R')


cat("#### Starting Simulation \n")
diffs = matrix(0, 9, simulations)
for (i in 1:simulations) {
  prior.gp = sim_gp(mu = 0, l = 1)
  post.gp = sim_gp(mu = 0, l = 1*cov_scale)
  
  # split em all
  sim.prior.split = vapply(1:100, function(x) matsplitter(prior.gp[,,x], 10, 10),
                           FUN.VALUE = array(0, dim = c(10, 10, 9)))
  
  sim.post.split = vapply(1:100, function(x) matsplitter(post.gp[,,x], 10, 10),
                          FUN.VALUE = array(0, dim = c(10, 10, 9)))
  
  
  # find the observed kst field
  kol.field = kst.field(sim.prior.split, sim.post.split, 200)
  
  # find the permutation distribution
  perm.fields = kst.permute(sim.prior.split, sim.post.split, 200, 200)
  
  # find the central regions
  perm.ed = edepth_set(perm.fields, depth_function = "rank")
  perm.cr = central_region(perm.fields, perm.ed)
  
  diffs[,i] = sapply(1:length(kol.field), function(x) as.integer(kol.field[x] > perm.cr[[2]][x]))
  cat(paste0("sim ", i, "\n"))
}

cat("#### Saving Data \n")
saveRDS(diffs, file = paste0("../../../outdata/cov_scale_", cov_scale,".rds"))

