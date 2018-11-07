### HERE I "INVENTED" THE HEIGHT VERSION OF THE TEST


rm(list = ls()); gc()
library(extdepth)
library(tictoc)
library(reshape2)
library(ggplot2)
library(e1071)

source("research/assimilation-cfr/sim/reference.R")

height = function(g, fmat) {
  
  # Computes the depth values of a function with respect to a set of functions (fmat)
  fn = ncol(fmat)
  depth = rep(0, length(g))
  
  for (row in 1:nrow(fmat)) {
    diff = sum(g[row] < fmat[row,])
    depth[row] = 1 - (diff / fn)
  }
  
  return(depth)
}
xheight = function(gmat, fmat) {
  apply(gmat, 2, function(x) mean(height(x, fmat)))
}
ks_cdf = function(x, n = 20) {
  if(x < 0.1) return(0)
  1 - 2*(sum(sapply(1:n, function(k) ((-1)^(k-1)) * exp(-2*(k^2)*(x^2)))))
}

f = gp1d(1000)
g = gp1d(1000)
fh = xheight(f, f)


fgg = melt(f)
fgg["xh"] = rep(fh, each = 50)

ggplot(fgg, aes(x=Var1, y=value, group=Var2, color=xh)) +
  geom_line() +
  theme_classic()

coverage = function(f, g) {
  ffxd = xheight(f, f)
  gfxd = xheight(g, f)
  
  tf = seq(0, 1, length.out = max(1000, 3*length(ffxd)))  
  ffr = sapply(tf, function(y) mean(ffxd <= y))
  gfr = sapply(tf, function(y) mean(gfxd <= y))
  
  plot(tf, ffr, type = "l")
  lines(tf, gfr, col = "red")
  
  1 - ks_cdf(sqrt(ncol(f))*abs(max(ffr - gfr)))
}


coverage = function(f, g) {
  fxd = sort(xheight(f, f))
  gxd = sort(xheight(g, f))
  
  fr = sapply(fxd, function(y) mean(fxd <= y))
  gr = sapply(fxd, function(y) mean(gxd <= y))
  
  qfr = qnorm(fr)
  qfr = qfr[is.finite(qfr)]
  
  qgr = qnorm(gr)
  qgr = qgr[is.finite(qgr)]
  
  t.test(qfr, qgr)$p.value
}



qqplot(ffr, gfr)

xn = 100
yn = 100


### SIZE
sims = 1000
for(s in 1:sims) {
  tic("Total")
  cat("Iteration: ", s, "\n")
  
  gp1 = gp1d(xn, mu = 0, pts = 200)
  gp2 = gp1d(yn, mu = 0, pts = 200)
  
  trevor[s] = coverage(gp1, gp2)
  regina[s] = quality(gp1, gp2)
  
  cat("Trevor: ", mean(trevor[1:s] < 0.05), "\n")
  cat("Regina: ", mean(regina[1:s] < 0.05), "\n")
  toc()
  cat("\n")
}


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
    gp1 = gp1d(xn, mu = par1[p], l = 30)
    gp2 = gp1d(yn, mu = par2[p], l = 30)
    
    trevor[p,s] = coverage(gp1, gp2)
    regina[p,s] = quality(gp1, gp2)
    
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

library(fdasrvf)

coverage = function(f, g) {
  f = gp1d(xn, sd = 1, l = 30)
  g = gp1d(yn, sd = 0.1, l = 30)
  h = gp1d(yn, sd = 10, l = 30)
  
  fxd = sort(xheight(f, f))
  gxd = sort(xheight(g, f))
  hxd = sort(xheight(h, f))
  
  fr = sapply(fxd, function(y) mean(fxd <= y))
  gr = sapply(fxd, function(y) mean(gxd <= y))
  hr = sapply(fxd, function(y) mean(hxd <= y))
  
  t = seq(0, 1, length.out = max(1000, 3*length(fxd)))
  fr = sapply(t, function(y) mean(fxd <= y))
  gr = sapply(t, function(y) mean(gxd <= y))
  hr = sapply(t, function(y) mean(hxd <= y))
  
  plot(fr, type = "l")
  lines(gr, col = "red")
  lines(hr, col = "blue")
  
  qfr = qnorm(fr)
  qfr = qfr[is.finite(qfr)]
  
  qgr = qnorm(gr)
  qgr = qgr[is.finite(qgr)]
  
  t.test(fr, gr)$p.value
  
  ncol(g)*mean((fr - gr)^2)
  ncol(h)*mean((fr - hr)^2)
  }  
  


coverage = function(f, g) {
  f = gp1d(xn, sd = 1, l = 30)
  g = gp1d(yn, sd = 0.1, l = 30)
  h = gp1d(yn, sd = 10, l = 30)
  
  fxd = sort(xheight(f, f))
  gxd = sort(xheight(g, f))
  hxd = sort(xheight(h, f))
  
  fr = sapply(fxd, function(y) mean(fxd <= y))
  gr = sapply(fxd, function(y) mean(gxd <= y))
  hr = sapply(fxd, function(y) mean(hxd <= y))
  
  t = seq(0, 1, length.out = max(1000, 3*length(fxd)))
  fr = sapply(t, function(y) mean(fxd <= y))
  gr = sapply(t, function(y) mean(gxd <= y))
  hr = sapply(t, function(y) mean(hxd <= y))
  
  plot(fr, type = "l")
  lines(gr, col = "red")
  lines(hr, col = "blue")
  
  qfr = qnorm(fr)
  qfr = qfr[is.finite(qfr)]
  
  qgr = qnorm(gr)
  qgr = qgr[is.finite(qgr)]
  
  t.test(fr, gr)$p.value
  
  ncol(g)*mean((fr - gr)^2)
}

xn = 100
yn = 100

#### SHIFT SCALE
par1 = rep(1, 10)
par2 = seq(0.5, 1.5, length.out = 10)

set.seed(1)

sims = 500
trevor = matrix(0, length(par1), sims)
regina = matrix(0, length(par2), sims)

for(p in 1:length(par2)) {
  tic("Total")
  cat("Scale shift: ", par2[p], "\n")
  for(s in 1:sims) {
    gp1 = gp1d(xn, sd = par1[p], l = 30)
    gp2 = gp1d(yn, sd = par2[p], l = 30)
    
    gp1 = gp1d(xn, sd = 1, l = 30)
    gp2 = gp1d(yn, sd = 0.5, l = 30)
    
    plt_funs(gp1, gp2)
    
    trevor[p,s] = coverage(gp1, gp2)
    regina[p,s] = quality(gp1, gp2)
    
  }
  toc()
  cat("\n")
}
scale_summary = data.frame(power = c(rowMeans(regina < 0.05),
                                     rowMeans(trevor < 0.05)),
                           method = c(rep("QI", length(par2)), 
                                      rep("KSD", length(par2))),
                           meanshift = rep(round(par2, 1), 2))

ggplot(scale_summary, aes(meanshift, power, color = method)) +
  geom_point() +
  geom_path() +
  geom_abline(intercept = 0.05, slope  = 0) +
  theme_classic() +
  ylab("Power") +
  xlab("Scale Shift") +
  ggtitle("Power against scale changes")