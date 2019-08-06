rm(list = ls()); gc()

library(ggplot2)
library(dplyr)
library(reshape2)
library(future)
library(future.apply)
# 
# source("research/temp/code/depth_tests.R")
# source("research/temp/code/depths.R")
# source("research/temp/code/simulation.R")

source("../research/assimilation-cfr/code/depths.R")
source("../research/assimilation-cfr/code/simulation.R")

ksd = function(f, g) {
  ff.xd = xdepth(f, f)
  fg.xd = xdepth(f, g)
  
  gg.xd = xdepth(g, g)
  gf.xd = xdepth(g, f)
  
  ff.cdf = sapply(sort(ff.xd), function(y) mean(ff.xd <= y))
  gf.cdf = sapply(sort(ff.xd), function(y) mean(gf.xd <= y))
  
  fg.cdf = sapply(sort(gg.xd), function(y) mean(fg.xd <= y))
  gg.cdf = sapply(sort(gg.xd), function(y) mean(gg.xd <= y))
  
  rate = sqrt((ncol(g)*ncol(f)) / (ncol(g) + ncol(f)))
  
  ksf = rate*max(abs(ff.cdf - gf.cdf))
  ksg = rate*max(abs(gg.cdf - fg.cdf))
  
  max(ksf, ksg)
}


ks_cdf = function(x, n = 100) {
  if(x < 0.05) return(0)
  1 - 2*(sum(sapply(1:n, function(k) ((-1)^(k-1)) * exp(-2*(k^2)*(x^2)))))
}

plan(multiprocess)
ksd.perm = function(f, g, perms=500) {
  prog = as.integer(quantile(1:perms, 1:20/20))
  
  h = cbind(f, g)
  hn = ncol(h)
  fn = ncol(f)
  
  ksd.dist = rep(0, perms)
  
  ksd.dist = future_sapply(1:perms, function(y) {
    hstar = h[,sample(1:hn, hn, replace = F)]
    ksd(hstar[,1:fn], hstar[,-(1:fn)])
  })
}


set.seed(1023)

n = 100

t = seq(0, 1, length.out = 50)

f = gp1d(fields = n, mu = sin(t*2*pi))
g = gp1d(fields = n+1, mu = sin(t*2*pi))


perm.table = ksd.perm(f, g, 10000)

t = seq(0, 2.1, length.out = 1000)
perm_cdf = sapply(t, function(x) mean(perm.table < x))
asym_cdf = sapply(t, function(x) (ks_cdf(x)))

cdf.gg = data.frame(Method = rep(c("Permutation", "Kolmogorov"), each=length(t)),
                    Value = c(t, t),
                    CDF = c(perm_cdf, asym_cdf))

ggplot(cdf.gg, aes(Value, CDF, color=Method)) +
  geom_line() +
  geom_vline(xintercept = t[min(which(perm_cdf > 0.95))], color="#00BFC4") +
  geom_vline(xintercept = t[min(which(asym_cdf > 0.95))], color="#F8766D") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5)) +
  ggtitle("Kolmogorov v.s. Permutation (n = 100, p = 10000)") +
  xlab("t") +
  ylab("P(K < t)")
# ggsave(paste0("research/assimilation-cfr/paper/misc/", "perm_dist.png"), width = 5, height = 3.2)





