rm(list = ls()); gc()


library(reshape2)
library(tictoc)
library(future)
library(future.apply)

# set to the top level folder
setwd("/Users/trevorh2/research/assimilation-cfr/submit/")

source("method/depth_tests.R")
source("method/depths.R")
source("method/simulation.R")

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

set.seed(072393)
perms = 500
sims = 10
n = c(25, 50, 75, 100)

t = seq(0, 2, length.out = 1000)
asym_cdf = sapply(t, function(x) ks_cdf(x))
asym_cval = quantile(asym_cdf, probs = c(0.90, 0.95, 0.99))

conv = matrix(NA, sims*length(n)^2, 6)

for (bat in 1:10) {
  for (rng in c(0.2, 0.3, 0.4, 0.5)) {
    for (nu in c(0.5, 1.0, 1.5)) {
      
      k = 1
      
      for(n1 in n) {
        tic(paste0("n1 = ", n1, " n2 = ", n1))
        
        for(j in 1:sims) {
          
          f = tp2d(fields = n1, df = 3, range = rng, nu = nu)
          g = tp2d(fields = n2, df = 3, range = rng, nu = nu)
          
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
      
      saveRDS(conv, file = paste0("out_nongp/conv_rng_",rng,"_nu_",nu,"_bat_",bat,"_diag.RDS"))
    }
  }
}
