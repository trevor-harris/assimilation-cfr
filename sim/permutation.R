### HERE I WAS TESTING THE CORRELATION BETWEEN THE KSF AND KSG

rm(list = ls()); gc()
library(extdepth)
library(tictoc)
library(ggplot2)

# source("research/assimilation-cfr/sim/reference.R")
source("../research/assimilation-cfr/sim/reference.R")

library(microbenchmark)

f = gp1d()
g = gp1d()

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
  
  ksf = rate*max(abs(ff.cdf - fg.cdf))
  ksg = rate*max(abs(gf.cdf - gg.cdf))
  
  max(ksf, ksg)
}

skt.boot = function(f, g, perms=500) {
  prog = as.integer(quantile(1:perms, 1:20/20))
  
  h = cbind(f, g)
  hn = ncol(h)
  fn = ncol(f)
  
  ksd.dist = rep(0, perms)
  
  tic()
  for(i in 1:perms) {
    hstar = h[,sample(1:hn, hn, replace = T)]
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

boot.pval = function(est, table) {
  mean(est > table)
}

f = gp1d(1000, pts = 500)
g = gp1d(1000, pts = 500)

boot.table = skt.boot(f, g, 1000)
hist(boot.table, breaks = 20)

t = seq(0, 2.1, length.out = 1000)
cov_cdf = sapply(t, function(x) mean(boot.table < x))
bro_cdf = sapply(t, function(x) (ks_cdf(x)))

cdf.gg = data.frame(Method = rep(c("Bootstrap", "Kolmogorov"), each=length(t)),
                    Value = c(t, t),
                    CDF = c(cov_cdf, bro_cdf))

ggplot(cdf.gg, aes(Value, CDF, color=Method)) +
  geom_line() +
  geom_vline(xintercept = t[min(which(cov_cdf > 0.95))], color="#00BFC4") +
  geom_vline(xintercept = t[min(which(bro_cdf > 0.95))], color="#F8766D") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5)) +
  ggtitle("Kolmogorov Distribution v.s. Bootstrap Distribution") +
  xlab("K Value") +
  ylab("Probability")
ggsave(paste0("../research/assimilation-cfr/paper/misc/", "bootstrap.png"), width = 5, height = 3.2)



