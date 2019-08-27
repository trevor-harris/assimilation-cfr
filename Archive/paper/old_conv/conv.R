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
  
  max(ksf, ksg)
}
kolm.perm = function(f, g, perms=1000) {
  prog = as.integer(quantile(1:perms, 1:20/20))
  
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
set.seed(072393)

sims = 50
ranges = c(5, 10, 15, 20, 25, 30)
nfuns = c(25, 50, 100, 150)

perms = 1000
t = seq(0, 2, length.out = 1000)
asym_cdf = sapply(t, function(x) ks_cdf(x))
asym_cval = quantile(asym_cdf, probs = c(0.90, 0.95, 0.99))


conv = data.frame(cdf_diff = numeric(0),
                  cval_90 = numeric(0),
                  cval_95 = numeric(0),
                  cval_99 = numeric(0),
                  nfun = numeric(0),
                  range = numeric(0))

k = 1

for(n in nfuns) {
  for(r in ranges) {
    tic(paste0("n = ", n, " r = ", r))
    
    for(j in 1:sims) {
      
      f = flatten(gp2d(n, sd = 1, l = r, pts = 20))
      g = flatten(gp2d(n+1, sd = 1, l = r, pts = 20))
      
      perm.table = kolm.perm(f, g, perms)
      
      perm_cdf = sapply(t, function(x) mean(perm.table < x))
      perm_cval = quantile(perm_cdf, probs = c(0.90, 0.95, 0.99))
      
      cdf_diff = mean((perm_cdf - asym_cdf)^2)
      cval_diff = perm_cval - asym_cval
      
      conv[k,] = c(cdf_diff, cval_diff, n, r)
      
      k = k + 1
    }
    
    toc()
  }
}


saveRDS(conv, file = "/home/trevorh2/assimilation-cfr/sim/conv/out/conv.RDS")
