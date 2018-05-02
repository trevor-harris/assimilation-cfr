rm(list = ls()); gc()

# get args
args = commandArgs(TRUE)
batch_no = as.double(args[1])
simulations = as.integer(args[2])
run = as.integer(args[3])

library(extdepth)
library(plgp)

source('../ks_field_functions.R')
source('../sim_functions.R')

# marginal number of points in the field (field is pts x pts)
pts = 40

# number of regions to subdivide the fields into
regions = 64

# standard flat prior mean
prior_mu = matrix(0, pts, pts)
post_mu = readRDS(paste0("../../simdata/run", run,"/post_mu.rds"))

prior_mu = as.vector(prior_mu)
post_mu = as.vector(post_mu)

cat("#### Starting Simulation \n")
diffs_de = matrix(0, regions, simulations)
diffs_bf = matrix(0, regions, simulations)

for (i in 1:simulations) {
  t0 = Sys.time()
  
  prior.gp = sim_gp(mu = prior_mu, l = 1, pts = pts)
  post.gp = sim_gp(mu = post_mu, l = 1, pts = pts)
  
  # split em all
  prior.gp.split = vapply(1:100, function(x) matsplitter(prior.gp[,,x], 5, 5),
                          FUN.VALUE = array(0, dim = c(5, 5, regions)))
  
  post.gp.split = vapply(1:100, function(x) matsplitter(post.gp[,,x], 5, 5),
                         FUN.VALUE = array(0, dim = c(5, 5, regions)))
  
  # find the observed kst field
  kol.field = kst.field(prior.gp.split, post.gp.split, 100)
  
  # find the permutation distribution
  perm.fields = kst.permute(prior.gp.split, post.gp.split, 100, 100)
  
  # Bonferroni central regions
  bf_val = (1-(0.05/regions))
  diffs_bf[,i] = sapply(1:length(kol.field), function(r) as.integer(kol.field[r] > quantile(perm.fields[r,], bf_val)))
  
  # Depth central regions
  perm.ed = edepth_set(perm.fields, depth_function = "rank")
  perm.cr = central_region(perm.fields, perm.ed)
  
  diffs_de[,i] = sapply(1:length(kol.field), function(x) as.integer(kol.field[x] > perm.cr[[2]][x]))
  
  cat(paste0("sim ", i, Sys.time()-t0, "\n"))
}

cat("#### Saving Data \n")
saveRDS(diffs_de, file = paste0("../../simdata/run", run,"/Depth", batch_no, ".rds"))
saveRDS(diffs_bf, file = paste0("../../simdata/run", run,"/Bonferroni", batch_no, ".rds"))


