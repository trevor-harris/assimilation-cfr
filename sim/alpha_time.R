rm(list = ls()); gc()

# get args
args = commandArgs(TRUE)
batch_no = as.double(args[1])
simulations = as.integer(args[2])

library(extdepth)

source('/home/trevorh2/assimilation-cfr/code/ks_field_functions.R')
source('/home/trevorh2/assimilation-cfr/code/sim_functions.R')

# marginal number of points in the field (field is pts x pts)
pts = 40

# number of regions to subdivide the fields into
regions = 64

# number of time points
time_points = 10

# standard flat prior mean
prior_mu = matrix(0, pts, pts)
post_mu = matrix(0, pts, pts)

prior_mu = as.vector(prior_mu)
post_mu = as.vector(post_mu)

cat("#### Starting Simulation \n")
upper_de = matrix(0, regions*time_points, simulations)
upper_bf = matrix(0, regions*time_points, simulations)
upper_pw = matrix(0, regions*time_points, simulations)
ks_value = matrix(0, regions*time_points, simulations)
pval_de = rep(0, simulations)

for (i in 1:simulations) {
  t0 = Sys.time()
  
  kst = rep(0, regions*time_points)
  ksp = matrix(0, regions*time_points, 1000)
  
  for (t in 1:time_points) {
    
    prior = sim_gp(100, mu = prior_mu, l = 5, pts = pts)
    post = sim_gp(100, mu = post_mu, l = 5, pts = pts)
    
    # split em all
    prior.split = vapply(1:100, function(x) matsplitter(prior[,,x], 5, 5),
                         FUN.VALUE = array(0, dim = c(5, 5, regions)))
    
    post.split = vapply(1:100, function(x) matsplitter(post[,,x], 5, 5),
                        FUN.VALUE = array(0, dim = c(5, 5, regions)))
    
    # find the observed kst field
    kst[1:regions + regions*(t - 1)] = kst.field(prior.split, post.split, 100)
    ksp[1:regions + regions*(t - 1),] = kst.permute(prior.split, post.split, 1000, 1)
    
  }
  
  # Bonferroni central regions
  bf_val = (1-(0.05/(regions*time_points)))
  
  # Depth central regions
  perm.ed = edepth_set(ksp, depth_function = "rank")
  perm.cr = central_region(ksp, perm.ed)
  
  kst.ed = edepth(kst, ksp, depth_function = "rank")
  
  # Upper regions
  upper_bf[,i] = sapply(1:length(kst), function(r) quantile(ksp[r,], bf_val))
  upper_pw[,i] = sapply(1:length(kst), function(r) quantile(ksp[r,], 0.95))
  upper_de[,i] = perm.cr[[2]]
  
  # save ks values
  ks_value[,i] = kst
  
  # save ks p values
  pval_de[i] = kst.ed
  
  cat(paste0("sim ", i, "\t", Sys.time()-t0, "\n"))
}

cat("#### Saving Data \n")
saveRDS(upper_de, file = paste0("/home/trevorh2/scratch/simdata/alpha/Depth", batch_no, ".rds"))
saveRDS(upper_bf, file = paste0("/home/trevorh2/scratch/simdata/alpha/Bonferroni", batch_no, ".rds"))
saveRDS(upper_pw, file = paste0("/home/trevorh2/scratch/simdata/alpha/Pointwise", batch_no, ".rds"))
saveRDS(ks_value, file = paste0("/home/trevorh2/scratch/simdata/alpha/Values", batch_no, ".rds"))
saveRDS(pval_de, file = paste0("/home/trevorh2/scratch/simdata/alpha/pvals", batch_no, ".rds"))

