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


cvmd = function(f, g) {
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
  
  rate = (ncol(g)*ncol(f)) / (ncol(g) + ncol(f))
  
  ksf = rate*mean((ffr - gfr)^2)
  ksg = rate*mean((fgr - ggr)^2)
  
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

cvmd.perm = function(f, g, perms=500) {
  h = cbind(f, g)
  
  hn = ncol(h)
  fn = ncol(f)
  
  cvmd.dist = rep(0, perms)
  
  for(i in 1:perms) {
    hstar = h[,sample(1:hn, hn)]
    cvmd.dist[i] = cvmd(hstar[,1:ncol(f)], hstar[,-(1:ncol(f))])
  }
  cvmd.dist
}

perm.pval = function(est, table) {
  mean(est > table)
}

f = gp1d(50)
g = gp1d(50)

perm.table = ksd.perm(f, g, 5000)
hist(perm.table)

cvmd.table = cvmd.perm(f, g, 2000)
plot(density(cvmd.table))

xn = 100
yn = 100

#### SHIFT MEANS
par1 = rep(0, 10)
par2 = seq(0, 1, length.out = 10)

set.seed(1)

sims = 500
trevor = matrix(0, length(par2), sims)
regina = matrix(0, length(par2), sims)

for(p in 1:length(par2)) {
  tic("Total")
  cat("Mean shift: ", par2[p], "\n")
  for(s in 1:sims) {
    gp1 = gp1d(xn, mu = par1[p], l = 10)
    gp2 = gp1d(yn, mu = par2[p], l = 10)
    
    trevor[p,s] = perm.pval(ksd(gp1, gp2), ksd.table)
    regina[p,s] = perm.pval(cvmd(gp1, gp2), cvmd.table)
    
  }
  toc()
  cat("\n")
}
location_summary = data.frame(power = c(rowMeans(regina < 0.05),
                                        rowMeans(trevor < 0.05)),
                              method = c(rep("QI", length(par2)), 
                                         rep("KSD", length(par2))),
                              meanshift = rep(round(par2, 1), 2))

ggplot(location_summary, aes(meanshift, power, color = method)) +
  geom_point() +
  geom_path() +
  geom_abline(intercept = 0.05, slope  = 0) +
  theme_classic() +
  ylab("Power") +
  xlab("Mean Shift") + 
  ggtitle("Power against mean changes")







ks_cdf = function(x, n = 20) {
  if(x < 0.1) return(0)
  1 - 2*(sum(sapply(1:n, function(k) ((-1)^(k-1)) * exp(-2*(k^2)*(x^2)))))
}

t = seq(0, 2.1, length.out = 1000)
cov_cdf = sapply(t, function(x) mean(perm.table < x))
bro_cdf = sapply(t, function(x) (ks_cdf(x)))

plot(t, cov_cdf, type = "l", main = "Theoretical (red) v.s. Observed (black)")
lines(t, bro_cdf, col = "red")




