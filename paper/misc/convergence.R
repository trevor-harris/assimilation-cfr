rm(list = ls()); gc()

library(ggplot2)
library(dplyr)
library(reshape2)

source("research/assimilation-cfr/code/depth_tests.R")
source("research/assimilation-cfr/code/depths.R")
source("research/assimilation-cfr/code/simulation.R")

kolm = function(f, g) {
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
    ksd.dist[i] = kolm(hstar[,1:ncol(f)], hstar[,-(1:ncol(f))])
    
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
n = 25
pts = 20
infil = 1
mu = 0
sd = 1
rng = 40
f = gp1d(fields = n, mu = 0, sd = 1, pts = pts, l = rng)
g = gp1d(fields = n+1, mu = mu, sd = sd, pts = pts, l = rng)

perm.table = skt.perm(f, g, 5000)

t = seq(0, 2.1, length.out = 1000)
perm_cdf = sapply(t, function(x) mean(perm.table < x))
asym_cdf = sapply(t, function(x) (ks_cdf(x)))

cdf.gg = data.frame(Method = rep(c("Permutation", "Kolmogorov"), each=length(t)),
                    Value = c(t, t),
                    CDF = c(perm_cdf, asym_cdf))

ggplot(cdf.gg, aes(Value, CDF, color=Method)) +
  geom_line(alpha = 1, size = 0.6) +
  geom_vline(xintercept = t[min(which(asym_cdf > 0.95))], color="#F8766D", alpha = 1, size = 0.5) +
  geom_vline(xintercept = t[min(which(perm_cdf > 0.95))], color="#00BFC4", alpha = 1, size = 0.5) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5)) +
  ggtitle("Kolmogorov v.s. Permutation Distribution (n = 25)") +
  xlab("t") +
  ylab("P(K < t)")
ggsave(paste0("research/assimilation-cfr/paper/misc/", "perm_dist25.png"), width = 5, height = 3.2)

set.seed(1023)
n = 50
pts = 20
infil = 1
mu = 0
sd = 1
rng = 40
f = gp1d(fields = n, mu = 0, sd = 1, pts = pts, l = rng)
g = gp1d(fields = n+1, mu = mu, sd = sd, pts = pts, l = rng)

perm.table = skt.perm(f, g, 5000)

t = seq(0, 2.1, length.out = 1000)
perm_cdf = sapply(t, function(x) mean(perm.table < x))
asym_cdf = sapply(t, function(x) (ks_cdf(x)))

cdf.gg = data.frame(Method = rep(c("Permutation", "Kolmogorov"), each=length(t)),
                    Value = c(t, t),
                    CDF = c(perm_cdf, asym_cdf))

ggplot(cdf.gg, aes(Value, CDF, color=Method)) +
  geom_line(alpha = 1, size = 0.6) +
  geom_vline(xintercept = t[min(which(asym_cdf > 0.95))], color="#F8766D", alpha = 1, size = 0.5) +
  geom_vline(xintercept = t[min(which(perm_cdf > 0.95))], color="#00BFC4", alpha = 1, size = 0.5) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5)) +
  ggtitle("Kolmogorov v.s. Permutation Distribution (n = 50)") +
  xlab("t") +
  ylab("P(K < t)")
ggsave(paste0("research/assimilation-cfr/paper/misc/", "perm_dist50.png"), width = 5, height = 3.2)


set.seed(1023)
n = 100
pts = 20
infil = 1
mu = 0
sd = 1
rng = 40
f = gp1d(fields = n, mu = 0, sd = 1, pts = pts, l = rng)
g = gp1d(fields = n+1, mu = mu, sd = sd, pts = pts, l = rng)

perm.table = skt.perm(f, g, 5000)

t = seq(0, 2.1, length.out = 1000)
perm_cdf = sapply(t, function(x) mean(perm.table < x))
asym_cdf = sapply(t, function(x) (ks_cdf(x)))

cdf.gg = data.frame(Method = rep(c("Permutation", "Kolmogorov"), each=length(t)),
                    Value = c(t, t),
                    CDF = c(perm_cdf, asym_cdf))

ggplot(cdf.gg, aes(Value, CDF, color=Method)) +
  geom_line(alpha = 1, size = 0.6) +
  geom_vline(xintercept = t[min(which(asym_cdf > 0.95))], color="#F8766D", alpha = 1, size = 0.5) +
  geom_vline(xintercept = t[min(which(perm_cdf > 0.95))], color="#00BFC4", alpha = 1, size = 0.5) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5)) +
  ggtitle("Kolmogorov v.s. Permutation Distribution (n = 100)") +
  xlab("t") +
  ylab("P(K < t)")
ggsave(paste0("research/assimilation-cfr/paper/misc/", "perm_dist100.png"), width = 5, height = 3.2)


set.seed(1023)
n = 200
pts = 20
infil = 1
mu = 0
sd = 1
rng = 40
f = gp1d(fields = n, mu = 0, sd = 1, pts = pts, l = rng)
g = gp1d(fields = n+1, mu = mu, sd = sd, pts = pts, l = rng)

perm.table = skt.perm(f, g, 5000)

t = seq(0, 2.1, length.out = 1000)
perm_cdf = sapply(t, function(x) mean(perm.table < x))
asym_cdf = sapply(t, function(x) (ks_cdf(x)))

cdf.gg = data.frame(Method = rep(c("Permutation", "Kolmogorov"), each=length(t)),
                    Value = c(t, t),
                    CDF = c(perm_cdf, asym_cdf))

ggplot(cdf.gg, aes(Value, CDF, color=Method)) +
  geom_line(alpha = 1, size = 0.6) +
  geom_vline(xintercept = t[min(which(asym_cdf > 0.95))], color="#F8766D", alpha = 1, size = 0.5) +
  geom_vline(xintercept = t[min(which(perm_cdf > 0.95))], color="#00BFC4", alpha = 1, size = 0.5) +
  theme_classic() +
  theme(plot.title = element_text(size=10)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  ggtitle("Kolmogorov v.s. Permutation Distribution (n = 200)") +
  xlab("t") +
  ylab("P(K < t)")
ggsave(paste0("research/assimilation-cfr/paper/misc/", "perm_dist200.png"), width = 5, height = 3.2)





