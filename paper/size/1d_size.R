rm(list = ls()); gc()

library(reshape2)

source("research/assimilation-cfr/code/simulation.R")
source("research/assimilation-cfr/code/depths.R")
source("research/assimilation-cfr/code/depth_tests.R")

pts = 25
l = 30
nfun = c(50, 100, 200, 300, 400, 500)
n = length(nfun)

sims = 500
iter = 10

ksize = matrix(0, iter, n)
qsize = matrix(0, iter, n)
for(k in 1:n) {
  klevel = rep(0, iter)
  qlevel = rep(0, iter)
  for(j in 1:iter) {
    kpval = rep(0, sims)
    qpval = rep(0, sims)
    for (i in 1:sims) {
      f = gp1d(nfun[k], pts = 25, l = l)
      g = gp1d(nfun[k]+1, pts = 25, l = l)
      
      kpval[i] = kolm(f, g)[2]
      qpval[i] = quality(f, g)[2]
    }
    klevel[j] = mean(kpval <= 0.05)
    qlevel[j] = mean(qpval <= 0.05)
  }
  ksize[,k] = klevel
  qsize[,k] = qlevel
}

ksize.gg = melt(ksize)
ksize.gg[["Stat"]] = "K"
ksize.gg[["Size"]] = ksize.gg[["value"]]
ksize.gg[["Functions"]] = as.factor(ksize.gg[["Var2"]])
levels(ksize.gg[["Functions"]]) = nfun

qsize.gg = melt(qsize)
qsize.gg[["Stat"]] = "Q"
qsize.gg[["Size"]] = qsize.gg[["value"]]
qsize.gg[["Functions"]] = as.factor(qsize.gg[["Var2"]])
levels(qsize.gg[["Functions"]]) = nfun

size = rbind(ksize.gg, qsize.gg)
size = size[,c("Stat", "Size", "Functions")]

ggplot(size, aes(x = Functions, y = Size, color = Stat)) + 
  geom_boxplot() +
  theme_classic() +
  geom_hline(yintercept = 0.05)

