rm(list = ls()); gc()

library(ggplot2)
library(dplyr)
library(reshape2)

source("../research/assimilation-cfr/code/depth_tests.R")
source("../research/assimilation-cfr/code/depths.R")
source("../research/assimilation-cfr/code/simulation.R")

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


#### Examples
r = 10

f = gp1d(25, sd = 1, l = r)
g = gp1d(25, sd = 1, l = r)

perm.table = kolm.perm(f, g, perms)
perm_cdf = sapply(t, function(x) mean(perm.table < x))

cdf1.gg = data.frame(Method = rep(c("Permutation", "Kolmogorov"), each=length(t)),
                     Value = c(t, t),
                     CDF = c(perm_cdf, asym_cdf))

f = gp1d(50, sd = 1, l = r)
g = gp1d(50, sd = 1, l = r)

perm.table = kolm.perm(f, g, perms)
perm_cdf = sapply(t, function(x) mean(perm.table < x))

cdf2.gg = data.frame(Method = rep(c("Permutation", "Kolmogorov"), each=length(t)),
                     Value = c(t, t),
                     CDF = c(perm_cdf, asym_cdf))

f = gp1d(150, sd = 1, l = r)
g = gp1d(150, sd = 1, l = r)

perm.table = kolm.perm(f, g, perms)
perm_cdf = sapply(t, function(x) mean(perm.table < x))


cdf3.gg = data.frame(Method = rep(c("Permutation", "Kolmogorov"), each=length(t)),
                     Value = c(t, t),
                     CDF = c(perm_cdf, asym_cdf))

cdf.gg = rbind(cdf1.gg, cdf2.gg, cdf3.gg)
cdf.gg[["NFUN"]] = factor(rep(c("N = 25", "N = 50", "N = 100"), each = 2000), 
                          levels = c("N = 25", "N = 50", "N = 100"))

ggplot(cdf.gg, aes(Value, CDF, color=Method)) +
  geom_line() +
  # geom_vline(xintercept = t[min(which(perm_cdf > 0.95))], color="#00BFC4") +
  # geom_vline(xintercept = t[min(which(asym_cdf > 0.95))], color="#F8766D") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5)) +
  ggtitle("Convergence of K") +
  xlab("t") +
  ylab("P(K < t)") +
  # theme(legend.position="none") +
  facet_grid(. ~ NFUN)
ggsave("../research/assimilation-cfr/paper/misc/perm_conv.png", width = 5, height = 3.2)





