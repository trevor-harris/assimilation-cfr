rm(list = ls()); gc()

# get args
args = commandArgs(TRUE)
mean_shift = as.double(args[1])
simulations = as.integer(args[2])

library(extdepth)
library(plgp)

source('../../ks_field_functions.R')
source('../../sim_functions.R')

prior_mu = matrix(0, 40, 40)
post_mu = kronecker(diag(1, 4, 4), matrix(mean_shift, 10, 10))

regions = 64

cat("#### Starting Simulation \n")
diffs = matrix(0, regions, simulations)
diffs_pw = matrix(0, regions, simulations)
diffs_bf = matrix(0, regions, simulations)

for (i in 1:simulations) {
  prior.gp = sim_gp(mu = prior_mu, l = 1, pts = 40)
  post.gp = sim_gp(mu = post_mu, l = 1, pts = 40)
  
  # split em all
  prior.gp.split = vapply(1:100, function(x) matsplitter(prior.gp[,,x], 5, 5),
                          FUN.VALUE = array(0, dim = c(5, 5, regions)))
  
  post.gp.split = vapply(1:100, function(x) matsplitter(post.gp[,,x], 5, 5),
                         FUN.VALUE = array(0, dim = c(5, 5, regions)))
  
  
  # find the observed kst field
  kol.field = kst.field(prior.gp.split, post.gp.split, 100)
  
  # find the permutation distribution
  perm.fields = kst.permute(prior.gp.split, post.gp.split, 100, 100)
  
  # PW central regions
  bf_val = (1-(0.05/regions))
  diffs_pw[,i] = sapply(1:length(kol.field), function(r) as.integer(kol.field[r] > quantile(perm.fields[r,], 0.95)))
  diffs_bf[,i] = sapply(1:length(kol.field), function(r) as.integer(kol.field[r] > quantile(perm.fields[r,], bf_val)))
  
  perm.ed = edepth_set(perm.fields, depth_function = "rank")
  perm.cr = central_region(perm.fields, perm.ed)
  
  diffs[,i] = sapply(1:length(kol.field), function(x) as.integer(kol.field[x] > perm.cr[[2]][x]))
  
  cat(paste0("sim ", i, "\n"))
}

cat("#### Saving Data \n")
saveRDS(diffs, file = paste0("../../../outdata/pointwise/depth_part_cr_high_", mean_shift,".rds"))
saveRDS(diffs_pw, file = paste0("../../../outdata/pointwise/pw_part_cr_high_", mean_shift,".rds"))
saveRDS(diffs_bf, file = paste0("../../../outdata/pointwise/bf_part_cr_high_", mean_shift,".rds"))
