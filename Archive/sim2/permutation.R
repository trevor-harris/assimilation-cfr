rm(list = ls()); gc()

library(ggplot2)
library(dplyr)
library(reshape2)

source("research/temp/code/depth_tests.R")
source("research/temp/code/depths.R")
source("research/temp/code/simulation.R")

kolm.wa = function(f, g) {
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

#### permute K
skt.perm = function(f, g, perms=500) {
  prog = as.integer(quantile(1:perms, 1:20/20))
  
  h = cbind(f, g)
  hn = ncol(h)
  fn = ncol(f)
  
  ksd.dist = rep(0, perms)
  
  tic()
  for(i in 1:perms) {
    hstar = h[,sample(1:hn, hn, replace = F)]
    ksd.dist[i] = kolm.wa(hstar[,1:ncol(f)], hstar[,-(1:ncol(f))])
    
    if (i %in% prog)  {
      cat(paste0(100*prog[which(prog == i)]/perms, "% done \n"))
      toc()
      cat("\n")
      tic()
    }
  }
  ksd.dist
}


set.seed(1023)
n = 300
pts = 30
infil = 2
mu = 0
sd = 1
rng = 100
f = gp1d(fields = n, mu = 0, sd = 1, pts = pts, l = rng)
g = gp1d(fields = n+1, mu = mu, sd = sd, pts = pts, l = rng)

f = apply(f, 2, function(x) smooth.spline(x, nknots = 20)$y)
g = apply(g, 2, function(x) smooth.spline(x, nknots = 20)$y)

f = apply(f, 2, function(x) spline(1:pts, x, n = infil*pts)$y)
g = apply(g, 2, function(x) spline(1:pts, x, n = infil*pts)$y)

plt_funs(f, g)

perm.table = skt.perm(f, g, 5000)

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
  ggtitle("Kolmogorov Distribution v.s. Permutation Distribution") +
  xlab("t") +
  ylab("P(K < t)")
ggsave(paste0("research/assimilation-cfr/paper/misc/", "perm_dist.png"), width = 5, height = 3.2)




