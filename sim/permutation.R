### HERE I WAS TESTING THE CORRELATION BETWEEN THE KSF AND KSG

rm(list = ls()); gc()
library(extdepth)
library(tictoc)
library(ggplot2)

source("research/assimilation-cfr/sim/reference.R")
# source("../research/assimilation-cfr/sim/reference.R")

library(microbenchmark)

#### bootstrap the KS
sk.test = function(f, g) {
  ff.xd = xdepth(f, f)
  fg.xd = xdepth(f, g)
  
  gg.xd = xdepth(g, g)
  gf.xd = xdepth(g, f)
  
  tf = seq(0, 1, length.out = max(1000, 3*length(ff.xd)))  
  ff.cdf = sapply(tf, function(y) mean(ff.xd <= y))
  gf.cdf = sapply(tf, function(y) mean(gf.xd <= y))
  
  tg = seq(0, 1, length.out = max(1000, 3*length(gg.xd)))
  fg.cdf = sapply(tg, function(y) mean(fg.xd <= y))
  gg.cdf = sapply(tg, function(y) mean(gg.xd <= y))
  
  rate = sqrt((ncol(g)*ncol(f)) / (ncol(g) + ncol(f)))
  
  # ks1 = max(abs(ff.cdf - fg.cdf))
  # ks2 = max(abs(gf.cdf - gg.cdf))
  # ks3 = rate*max(abs(ff.cdf - gf.cdf))
  # ks4 = rate*max(abs(gg.cdf - fg.cdf))
  
  ks3 = rate*max(abs(ff.cdf - gf.cdf))
  ks4 = rate*max(abs(gg.cdf - fg.cdf))
  
  # ks1 = sqrt(ncol(f))*max(abs(ff.cdf - fg.cdf))
  # ks2 = sqrt(ncol(g))*max(abs(gf.cdf - gg.cdf))
  
  max(ks3, ks4)
}
skt.perm = function(f, g, perms=500) {
  prog = as.integer(quantile(1:perms, 1:20/20))
  
  h = cbind(f, g)
  hn = ncol(h)
  fn = ncol(f)
  
  ksd.dist = rep(0, perms)
  
  tic()
  for(i in 1:perms) {
    hstar = h[,sample(1:hn, hn, replace = F)]
    ksd.dist[i] = sk.test(hstar[,1:ncol(f)], hstar[,-(1:ncol(f))])
    
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
n = 400
pts = 10
infil = 5
mu = 0
sd = 1
rng = 10
f = gp1d(fields = n, mu = 0, sd = 1, pts = pts, l = rng)
g = gp1d(fields = n, mu = mu, sd = sd, pts = pts, l = rng)

f = apply(f, 2, function(x) smooth.spline(x, nknots = 10)$y)
g = apply(g, 2, function(x) smooth.spline(x, nknots = 10)$y)

f = apply(f, 2, function(x) spline(1:pts, x, n = infil*pts)$y)
g = apply(g, 2, function(x) spline(1:pts, x, n = infil*pts)$y)

# plt_funs(f, g)

perm.table = skt.perm(f, g, 1000)

# rate1 = sqrt((ncol(g)*ncol(f)) / (ncol(g) + ncol(f)))
# rate2 = sqrt(ncol(f))
# perm.table2 = rate1*perm.table

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

# ggsave(paste0("../research/assimilation-cfr/paper/misc/", "bootstrap.png"), width = 5, height = 3.2)




