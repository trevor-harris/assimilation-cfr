rm(list = ls()); gc()

library(reshape2)
library(tictoc)
library(future)
library(future.apply)

source("/home/trevorh2/assimilation-cfr/code/depth_tests.R")
source("/home/trevorh2/assimilation-cfr/code/depths.R")
source("/home/trevorh2/assimilation-cfr/code/simulation.R")

kolm = function(f, g) {
  ff.xd = xdepth(f, f)
  fg.xd = xdepth(f, g)
  
  gg.xd = xdepth(g, g)
  gf.xd = xdepth(g, f)
  
  ff.cdf = sapply(ff.xd, function(y) mean(ff.xd <= y))
  gf.cdf = sapply(ff.xd, function(y) mean(gf.xd <= y))
  fg.cdf = sapply(gg.xd, function(y) mean(fg.xd <= y))
  gg.cdf = sapply(gg.xd, function(y) mean(gg.xd <= y))
  
  rate = sqrt((ncol(g)*ncol(f)) / (ncol(g) + ncol(f)))
  
  ksf = rate*max(abs(ff.cdf - gf.cdf))
  ksg = rate*max(abs(gg.cdf - fg.cdf))
  
  ks = max(ksf, ksg)
  ks
}
kolm.perm = function(f, g, perms=500) {
  h = cbind(f, g)
  hn = ncol(h)
  fn = ncol(f)
  
  ksd.dist = rep(0, perms)
  
  ksd.dist = future_sapply(1:perms, function(y) {
    hstar = h[,sample(1:hn, hn, replace = F)]
    kolm(hstar[,1:fn], hstar[,-(1:fn)])
  })
}

plan(multiprocess)

#### Convergence
# get args
args = commandArgs(TRUE)
rng = as.double(args[1])
nu = as.double(args[2])
bat = as.integer(args[3])
seed = (as.integer(args[4]) + 100*rng + 100*nu + 100*rng*nu) * bat

set.seed(seed)
perms = 500
sims = 10
n = c(25, 50, 75, 100)

t = seq(0, 2, length.out = 1000)
asym_cdf = sapply(t, function(x) ks_cdf(x))
asym_cval = quantile(asym_cdf, probs = c(0.90, 0.95, 0.99))

conv = matrix(NA, sims*length(n)^2, 6)

k = 1
for(n1 in n) {
    tic(paste0("n1 = ", n1, " n2 = ", n1))
    
    for(j in 1:sims) {
      
      f = gp2d.mat(fields = n1, range = rng, nu = nu)
      g = gp2d.mat(fields = n1, range = rng, nu = nu)
      
      f = flatten(f)
      g = flatten(g)
      
      perm.table = kolm.perm(f, g, perms)
      
      perm_cdf = sapply(t, function(x) mean(perm.table < x))
      perm_cval = quantile(perm_cdf, probs = c(0.90, 0.95, 0.99))
      
      cdf_diff = mean((perm_cdf - asym_cdf)^2)
      cval_diff = perm_cval - asym_cval
      
      conv[k,] = c(cdf_diff, cval_diff, n1, n1)
      
      k = k + 1
    }
    toc()
}

conv = data.frame(n1 = conv[,5],
                  n2 = conv[,6],
                  range = rng,
                  nu = nu,
                  cdf_diff = conv[,1],
                  cval_90 = conv[,2],
                  cval_95 = conv[,3],
                  cval_99 = conv[,4])

cat("Saving...")
# saveRDS(conv, file = "/home/trevorh2/assimilation-cfr/sim/conv/out/conv.RDS")
saveRDS(conv, file = paste0("/home/trevorh2/matern/conv/out/conv_r_",rng,"_nu_",nu,"_bat_",bat,"_diag.RDS"))
cat("Done")
