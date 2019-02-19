
rm(list = ls()); gc()

library(reshape2)
library(tictoc)
library(future)
library(future.apply)

source("../research/assimilation-cfr/code/simulation.R")
source("../research/assimilation-cfr/code/depths.R")
source("../research/assimilation-cfr/code/depth_tests.R")

set.seed(042696)


##### multicore version
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

n = 50
perms = 1000
ranges = c(10, 20, 30, 40)
stdevs = c(0.1, 0.5, 1, 2, 3)
nfuns = c(50, 100, 200)
t = seq(0, 1, length.out = 50)

diff_cdf = data.frame(size = numeric(0),
                      nfun = numeric(0),
                      range = numeric(0),
                      std = numeric(0))

plan(multiprocess)

k = 1

for(n in nfuns) {
  for(r in ranges) {
    for (s in stdevs) {
      tic(paste0("n = ", n, " r = ", r, " s = ", s))
      
      f = gp1d(n, sd = s, l = r)
      g = gp1d(n+1, sd = s, l = r)
      
      perm.table = kolm.perm(f, g, perms)
      
      t = seq(0, 2.1, length.out = 1000)
      perm_cdf = sapply(t, function(x) mean(perm.table < x))
      asym_cdf = sapply(t, function(x) (ks_cdf(x)))
      
      diff_cdf[k,] = c(mean((perm_cdf - asym_cdf)^2), n, r, s)
      toc()
      
      k = k + 1
    }
  }
}


