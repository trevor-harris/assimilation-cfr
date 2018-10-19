### HERE I WAS TESTING THE CORRELATION BETWEEN THE KSF AND KSG


rm(list = ls()); gc()
library(extdepth)
library(tictoc)
library(ggplot2)

source("research/assimilation-cfr/sim/reference.R")

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
  
  ksf = max(abs(ffr - gfr))
  ksg = max(abs(fgr - ggr))
  
  return(list(ksf = ksf, ksg = ksg))
}

ks_cdf = function(x, n = 20) {
  if(x < 0.1) return(0)
  1 - 2*(sum(sapply(1:n, function(k) ((-1)^(k-1)) * exp(-2*(k^2)*(x^2)))))
}


sims = 1000

ksf = rep(0, sims)
ksg = rep(0, sims)

for(s in 1:sims) {
  tic("Total")
  cat("Simulation ", s, "\n")
  
  gp1 = gp1d(1000)
  gp2 = gp1d(1000)
  
  ks = ksd(gp1, gp2)
  
  ksf[s] = ks$ksf
  ksg[s] = ks$ksg
  
  toc()
  cat("\n")
}

ksf1 = ksf
ksg1 = ksg

equal500 = cor.test(ksf, ksg)
plot(ksf, ksg, main = paste0("Corr: ", round(equal500$estimate, 3)))

f = gp1
g = gp2

rate = sqrt((ncol(g)*ncol(f)) / (ncol(g) + ncol(f)))
cover = rate*sapply(1:1000, function(x) max(ksf[x], ksg[x]))

t = seq(0, 2, length.out = 1000)
cov_cdf = sapply(t, function(x) mean(cover <= x))
bro_cdf = sapply(t, function(x) (ks_cdf(x)))
plot(t, cov_cdf, type = "l", main = "Theoretical (red) v.s. Observed (black)")
lines(t, bro_cdf, col = "red")







ksf = rep(0, sims)
ksg = rep(0, sims)

for(s in 1:sims) {
  tic("Total")
  cat("Simulation ", s, "\n")
  
  gp1 = gp1d(1000)
  gp2 = gp1d(1000)
  
  ks = ksd(gp1, gp2)
  
  ksf[s] = ks$ksf
  ksg[s] = ks$ksg
  
  toc()
  cat("\n")
}

ksf2 = ksf
ksg2 = ksg

ksfbig = c(ksf1, ksf2)
ksgbig = c(ksg1, ksg2)



equal500 = cor.test(ksfbig, ksgbig)
plot(ksf, ksg, main = paste0("Corr: ", round(equal500$estimate, 3)))

f = gp1
g = gp2

rate = sqrt((ncol(g)*ncol(f)) / (ncol(g) + ncol(f)))
cover = rate*sapply(1:2000, function(x) max(ksfbig[x], ksgbig[x]))

t = seq(0, 2, length.out = 1000)
cov_cdf = sapply(t, function(x) mean(cover <= x))
bro_cdf = sapply(t, function(x) (ks_cdf(x)))
plot(t, cov_cdf, type = "l", main = "Theoretical (red) v.s. Observed (black)")
lines(t, bro_cdf, col = "red")





sims = 2500

ksf = rep(0, sims)
ksg = rep(0, sims)

for(s in 1:sims) {
  tic("Total")
  cat("Simulation ", s, "\n")
  
  gp1 = gp1d(200)
  gp2 = gp1d(50)
  
  ks = ksd(gp1, gp2)
  
  ksf[s] = ks$ksf
  ksg[s] = ks$ksg
  
  toc()
  cat("\n")
}

plot(ksf, ksg)
nonequal50 = cor.test(ksf, ksg)






sims = 2500

ksf = rep(0, sims)
ksg = rep(0, sims)

for(s in 1:sims) {
  tic("Total")
  cat("Simulation ", s, "\n")
  
  gp1 = flatten(gp2d(100, pts = 20))
  gp2 = flatten(gp2d(100, pts = 20))
  
  ks = ksd(gp1, gp2)
  
  ksf[s] = ks$ksf
  ksg[s] = ks$ksg
  
  toc()
  cat("\n")
}

plot(ksf, ksg)
nonequal50 = cor.test(ksf, ksg)


ks_cdf = function(x, n = 20) {
  if(x < 0.1) return(0)
  1 - 2*(sum(sapply(1:n, function(k) ((-1)^(k-1)) * exp(-2*(k^2)*(x^2)))))
}

t = seq(0, 2, length.out = 1000)

bro_cdf = sapply(t, function(x) (ks_cdf(x)))
plot(bro_cdf, type = "l")


library(e1071)

sims = 25

cover = rep(0, sims)
brown = rep(0, sims)

for(s in 1:sims) {
  tic("Total")
  cat("Simulation ", s, "\n")
  
  cover[s] = coverage(gp1d(), gp1d())
  brown[s] = e1071::r
  
  toc()
  cat("\n")
}





