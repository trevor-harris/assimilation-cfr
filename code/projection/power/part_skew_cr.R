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

mean_shift = 0.2
prior_mu = matrix(0, 30, 30)
post_mu = kronecker(diag(1, 3, 3), matrix(mean_shift, 11, 12))
post_mu = post_mu[1:30,1:30]

prior_mu = as.vector(prior_mu)
post_mu = as.vector(post_mu)

cat("#### Starting Simulation \n")
diffs = matrix(0, 9, simulations)
diffs_pw = matrix(0, 9, simulations)
diffs_bf = matrix(0, 9, simulations)

for (i in 1:simulations) {
  prior.gp = sim_gp(mu = prior_mu, l = 1)
  post.gp = sim_gp(mu = post_mu, l = 1)
  
  # split em all
  prior.gp.split = vapply(1:100, function(x) matsplitter(prior.gp[,,x], 10, 10),
                          FUN.VALUE = array(0, dim = c(10, 10, 9)))
  
  post.gp.split = vapply(1:100, function(x) matsplitter(post.gp[,,x], 10, 10),
                         FUN.VALUE = array(0, dim = c(10, 10, 9)))
  
  
  # find the observed kst field
  kol.field = kst.field(prior.gp.split, post.gp.split, 200)
  
  # find the permutation distribution
  perm.fields = kst.permute(prior.gp.split, post.gp.split, 200, 200)
  
  # PW central regions
  bf_val = (1-(0.05/9))
  diffs_pw[,i] = sapply(1:length(kol.field), function(r) as.integer(kol.field[r] > quantile(perm.fields[r,], 0.95)))
  diffs_bf[,i] = sapply(1:length(kol.field), function(r) as.integer(kol.field[r] > quantile(perm.fields[r,], bf_val)))
  
  perm.ed = edepth_set(perm.fields, depth_function = "rank")
  perm.cr = central_region(perm.fields, perm.ed)
  
  diffs[,i] = sapply(1:length(kol.field), function(x) as.integer(kol.field[x] > perm.cr[[2]][x]))
  
  cat(paste0("sim ", i, "\n"))
  cat(diffs_bf[,i], "\n")
}

cat("#### Saving Data \n")
saveRDS(diffs, file = paste0("../../../outdata/pointwise/depth_skew_cr_", mean_shift,".rds"))
saveRDS(diffs_pw, file = paste0("../../../outdata/pointwise/pw_skew_cr_", mean_shift,".rds"))
saveRDS(diffs_bf, file = paste0("../../../outdata/pointwise/bf_skew_cr_", mean_shift,".rds"))
