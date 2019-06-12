
rm(list = ls()); gc()

library(reshape2)
library(ggplot2)
library(tictoc)
library(future)
library(future.apply)

source("../research/assimilation-cfr/code/simulation.R")
source("../research/assimilation-cfr/code/depths.R")
source("../research/assimilation-cfr/code/depth_tests.R")

set.seed(042696)


##### multicore version
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

#### settings
perms = 2500
ranges = c(5, 10, 20, 30, 40)
stdevs = c(0.1, 0.5, 1, 2, 3)
nfuns = c(25, 50, 100, 200)
t = seq(0, 1, length.out = 50)

mse_cdf = data.frame(size = numeric(0),
                      nfun = numeric(0),
                      range = numeric(0))

kol_cdf = data.frame(size = numeric(0),
                     nfun = numeric(0),
                     range = numeric(0))

plan(multiprocess)

t = seq(0, 2.1, length.out = 1000)
asym_cdf = sapply(t, function(x) (ks_cdf(x)))

k = 1

for(i in 1:20) {
  for(n in nfuns) {
    for(r in ranges) {
      tic(paste0("n = ", n, " r = ", r))
      
      f = gp1d(n, sd = 1, l = r)
      g = gp1d(n+1, sd = 1, l = r)
      
      perm.table = kolm.perm(f, g, perms)
      perm_cdf = sapply(t, function(x) mean(perm.table < x))
      
      mse_cdf[k,] = c(mean((perm_cdf - asym_cdf)^2), n, r)
      kol_cdf[k,] = c(max(abs(perm_cdf - asym_cdf)), n, r)
      toc()
      
      k = k + 1
    }
  }
}

diffs = kol_cdf
diffs[["samples"]] = as.factor(diffs[["nfun"]])
diffs[["Range"]] = as.factor(diffs[["range"]])
diffs[["Error"]] = diffs[["size"]]

ggplot(diffs, aes(x = samples, y = Error, fill = Range)) + 
  geom_boxplot() +
  theme_classic() +
  facet_grid(. ~ Range) +
  xlab("Sample sizes of X and Y") +
  ylab("Maximum Error")
ggsave("../research/assimilation-cfr/paper/misc/max_error.png", width = 6.3, height = 3.2)

diffs = mse_cdf
diffs[["samples"]] = as.factor(diffs[["nfun"]])
diffs[["Range"]] = as.factor(diffs[["range"]])
diffs[["Error"]] = diffs[["size"]]

ggplot(diffs, aes(x = samples, y = Error, fill = Range)) + 
  geom_boxplot() +
  theme_classic() +
  facet_grid(. ~ Range) +
  xlab("Sample sizes of X and Y") +
  ylab("Mean Squared Error")
ggsave("../research/assimilation-cfr/paper/misc/mse_error.png", width = 6.3, height = 3.2)

# for(n in nfuns) {
#   for(r in ranges) {
#     for (s in stdevs) {
#       tic(paste0("n = ", n, " r = ", r, " s = ", s))
#       
#       f = gp1d(n, sd = s, l = r)
#       g = gp1d(n+1, sd = s, l = r)
#       
#       perm.table = kolm.perm(f, g, perms)
#       perm_cdf = sapply(t, function(x) mean(perm.table < x))
#       
#       mse_cdf[k,] = c(mean((perm_cdf - asym_cdf)^2), n, r, s)
#       kol_cdf[k,] = c(max(abs(perm_cdf - asym_cdf)), n, r, s)
#       toc()
#       
#       k = k + 1
#     }
#   }
# }


perms = 2500
t = seq(0, 2.1, length.out = 1000)
asym_cdf = sapply(t, function(x) (ks_cdf(x)))





#### Examples
r = 30

f = gp1d(25, sd = 1, l = r)
g = gp1d(26, sd = 1, l = r)

perm.table = kolm.perm(f, g, perms)
perm_cdf = sapply(t, function(x) mean(perm.table < x))

cdf.gg = data.frame(Method = rep(c("Permutation", "Kolmogorov"), each=length(t)),
                    Value = c(t, t),
                    CDF = c(perm_cdf, asym_cdf))

ggplot(cdf.gg, aes(Value, CDF, color=Method)) +
  geom_line() +
  geom_vline(xintercept = t[min(which(perm_cdf > 0.95))], color="#00BFC4") +
  geom_vline(xintercept = t[min(which(asym_cdf > 0.95))], color="#F8766D") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5)) +
  ggtitle("N = 25") +
  xlab("t") +
  ylab("P(K < t)") +
  theme(legend.position="none")
# ggsave("../research/assimilation-cfr/paper/misc/perm_dist25.png", width = 5, height = 3.2)



f = gp1d(50, sd = 1, l = r)
g = gp1d(51, sd = 1, l = r)

perm.table = kolm.perm(f, g, perms)
perm_cdf = sapply(t, function(x) mean(perm.table < x))

cdf.gg = data.frame(Method = rep(c("Permutation", "Kolmogorov"), each=length(t)),
                    Value = c(t, t),
                    CDF = c(perm_cdf, asym_cdf))

ggplot(cdf.gg, aes(Value, CDF, color=Method)) +
  geom_line() +
  geom_vline(xintercept = t[min(which(perm_cdf > 0.95))], color="#00BFC4") +
  geom_vline(xintercept = t[min(which(asym_cdf > 0.95))], color="#F8766D") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5)) +
  ggtitle("N = 50") +
  xlab("t") +
  ylab("P(K < t)") +
  theme(legend.position="none")
# ggsave("../research/assimilation-cfr/paper/misc/perm_dist50.png", width = 5, height = 3.2)




f = gp1d(100, sd = s, l = r)
g = gp1d(101, sd = s, l = r)

perm.table = kolm.perm(f, g, perms)
perm_cdf = sapply(t, function(x) mean(perm.table < x))

cdf.gg = data.frame(Method = rep(c("Permutation", "Kolmogorov"), each=length(t)),
                    Value = c(t, t),
                    CDF = c(perm_cdf, asym_cdf))

ggplot(cdf.gg, aes(Value, CDF, color=Method)) +
  geom_line() +
  geom_vline(xintercept = t[min(which(perm_cdf > 0.95))], color="#00BFC4") +
  geom_vline(xintercept = t[min(which(asym_cdf > 0.95))], color="#F8766D") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5)) +
  ggtitle("N = 100") +
  xlab("t") +
  ylab("P(K < t)") +
  theme(legend.position="none")
# ggsave("../research/assimilation-cfr/paper/misc/perm_dist100.png", width = 5, height = 3.2)




f = gp1d(200, sd = s, l = r)
g = gp1d(201, sd = s, l = r)

perm.table = kolm.perm(f, g, perms)
perm_cdf = sapply(t, function(x) mean(perm.table < x))

cdf.gg = data.frame(Method = rep(c("Permutation", "Kolmogorov"), each=length(t)),
                    Value = c(t, t),
                    CDF = c(perm_cdf, asym_cdf))

ggplot(cdf.gg, aes(Value, CDF, color=Method)) +
  geom_line() +
  geom_vline(xintercept = t[min(which(perm_cdf > 0.95))], color="#00BFC4") +
  geom_vline(xintercept = t[min(which(asym_cdf > 0.95))], color="#F8766D") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5)) +
  ggtitle("N = 200") +
  xlab("t") +
  ylab("P(K < t)") +
  theme(legend.position="none")
# ggsave("../research/assimilation-cfr/paper/misc/perm_dist200.png", width = 5, height = 3.2)


f = gp1d(500, sd = s, l = 30)
g = gp1d(501, sd = s, l = 30)

perm.table = kolm.perm(f, g, perms)
perm_cdf = sapply(t, function(x) mean(perm.table < x))

cdf.gg = data.frame(Method = rep(c("Permutation", "Kolmogorov"), each=length(t)),
                    Value = c(t, t),
                    CDF = c(perm_cdf, asym_cdf))

ggplot(cdf.gg, aes(Value, CDF, color=Method)) +
  geom_line() +
  geom_vline(xintercept = t[min(which(perm_cdf > 0.95))], color="#00BFC4") +
  geom_vline(xintercept = t[min(which(asym_cdf > 0.95))], color="#F8766D") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5)) +
  ggtitle("N = 500") +
  xlab("t") +
  ylab("P(K < t)") +
  theme(legend.position="none")
