rm(list = ls()); gc()

# get args
args = commandArgs(TRUE)
mean_shift = as.double(args[1])
simulations = as.integer(args[2])

library(extdepth)
library(dplyr)
library(reshape2)
library(plgp)

source('../../ks_field_functions.R')
source('../../sim_functions.R')

prior_mu = 0

cat("#### Starting Simulation \n")
diffs = matrix(0, 9, simulations)
diffs_pw = matrix(0, 9, simulations)
diffs_bf = matrix(0, 9, simulations)

for (i in 1:simulations) {
  prior.gp = sim_gp(mu = prior_mu, l = 1)
  post.gp = sim_gp(mu = prior_mu + mean_shift, l = 1)
  
  # split em all
  prior.gp.split = vapply(1:100, function(x) matsplitter(prior.gp[,,x], 10, 10),
                           FUN.VALUE = array(0, dim = c(10, 10, 9)))
  
  post.gp.split = vapply(1:100, function(x) matsplitter(post.gp[,,x], 10, 10),
                          FUN.VALUE = array(0, dim = c(10, 10, 9)))
  
  
  # find the observed kst field
  kol.field = kst.field(prior.gp.split, post.gp.split)
  
  # find the permutation distribution
  perm.fields = kst.permute(prior.gp.split, post.gp.split, 100, 100)
  
  # PW central regions
  bf_val = (1-(0.05/9))
  diffs_pw[,i] = sapply(1:length(kol.field), function(r) as.integer(kol.field[r] > quantile(perm.fields[r,], 0.95)))
  diffs_bf[,i] = sapply(1:length(kol.field), function(r) as.integer(kol.field[r] > quantile(perm.fields[r,], bf_val)))
  
  perm.ed = edepth_set(perm.fields)
  perm.cr = central_region(perm.fields, perm.ed)
  
  diffs[,i] = sapply(1:length(kol.field), function(x) as.integer(kol.field[x] > perm.cr[[2]][x]))
  
  cat(paste0("sim ", i, "\n"))
  cat(diffs_bf[,i], "\n")
}

cat("#### Saving Data \n")
saveRDS(diffs, file = paste0("../../../outdata/depth_cr_", mean_shift,".rds"))
saveRDS(diffs_pw, file = paste0("../../../outdata/pw_cr_", mean_shift,".rds"))
saveRDS(diffs_bf, file = paste0("../../../outdata/bf_cr_", mean_shift,".rds"))
