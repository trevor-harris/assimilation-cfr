### HERE I WAS TESTING THE CORRELATION BETWEEN THE KSF AND KSG

rm(list = ls()); gc()
library(extdepth)
library(tictoc)
library(ggplot2)

source("research/assimilation-cfr/sim/reference.R")

library(microbenchmark)

#### bootstrap the KS
ksd = function(f, g) {
  ffxd = xdepth(f, f)
  gfxd = xdepth(g, f)
  fgxd = xdepth(f, g)
  ggxd = xdepth(g, g)
  
  tf = seq(0, 1, length.out = max(1000, 3*length(ffxd)))  
  ffr = sapply(tf, function(y) mean(ffxd <= y))
  gfr = sapply(tf, function(y) mean(gfxd <= y))
  
  tg = seq(0, 1, length.out = max(1000, 3*length(ggxd)))
  fgr = sapply(tg, function(y) mean(fgxd <= y))
  ggr = sapply(tg, function(y) mean(ggxd <= y))
  
  rate = sqrt((ncol(g)*ncol(f)) / (ncol(g) + ncol(f)))
  
  ksf = rate*max(abs(ffr - gfr))
  ksg = rate*max(abs(fgr - ggr))
  
  max(ksf, ksg)
}

ksd.perm = function(f, g, perms=500) {
  h = cbind(f, g)
  
  hn = ncol(h)
  fn = ncol(f)
  
  ksd.dist = rep(0, perms)
  
  for(i in 1:perms) {
    hstar = h[,sample(1:hn, hn)]
    ksd.dist[i] = ksd(hstar[,1:ncol(f)], hstar[,-(1:ncol(f))])
  }
  ksd.dist
}

f = gp1d(100)
g = gp1d(100)

perm.table = ksd.perm(f, g, 5000)
hist(perm.table)


ks_cdf = function(x, n = 20) {
  if(x < 0.1) return(0)
  1 - 2*(sum(sapply(1:n, function(k) ((-1)^(k-1)) * exp(-2*(k^2)*(x^2)))))
}

t = seq(0, 2.1, length.out = 1000)
cov_cdf = sapply(t, function(x) mean(perm.table < x))
bro_cdf = sapply(t, function(x) (ks_cdf(x)))

plot(t, cov_cdf, type = "l", main = "Theoretical (red) v.s. Observed (black)")
lines(t, bro_cdf, col = "red")




