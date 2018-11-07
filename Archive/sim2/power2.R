rm(list = ls()); gc()
library(extdepth)
library(tictoc)
library(ggplot2)

source("research/assimilation-cfr/sim/reference.R")
quality = function(f, g) {
  ffxd = xdepth(f, f)
  gfxd = xdepth(g, f)
  rf = sapply(gfxd, function(y) mean(y <= ffxd))
  
  ggxd = xdepth(g, g)
  fgxd = xdepth(f, g)
  rg = sapply(fgxd, function(y) mean(y <= ggxd))
  
  r = max(mean(rf), mean(rg))
  
  1 - pnorm((r - 0.5) / sqrt((1/ncol(f) + 1/ncol(g))/12))^2
}

xn = 100
yn = 100

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
    gp1 = gp1d(xn, mu = par1[p], l = 10)
    gp2 = gp1d(yn, mu = par2[p], l = 10)
    
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
    gp1 = gp1d(xn, sd = par1[p], l = 10)
    gp2 = gp1d(yn, sd = par2[p], l = 10)
    
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





#### SHIFT CORR
par1 = rep(10, 10)
par2 = seq(1, 50, length.out = 10)

set.seed(1)

sims = 500
trevor = matrix(0, length(par1), sims)
regina = matrix(0, length(par2), sims)

for(p in 1:length(par2)) {
  tic("Total")
  cat("Correlation shift: ", par2[p], "\n")
  for(s in 1:sims) {
    gp1 = gp1d(xn, l = par1[p])
    gp2 = gp1d(yn, l = par2[p])
    
    trevor[p,s] = coverage(gp1, gp2)
    regina[p,s] = quality(gp1, gp2)
    
  }
  toc()
  cat("\n")
}
corr_summary = data.frame(power = c(rowMeans(regina < 0.05),
                                    rowMeans(trevor < 0.05)),
                          method = c(rep("QI", length(par2)), 
                                     rep("KSD", length(par2))),
                          meanshift = rep(round(par2, 1), 2))

ggplot(corr_summary, aes(meanshift, power, color = method)) +
  geom_point() +
  geom_path() +
  theme_classic() +
  ylab("Power") +
  xlab("Correlation") +
  ggtitle("Power against correlation changes")




#### SHIFT SHAPE
par1 = rep(0, 10)
par2 = seq(-3, 3, length.out = 10)

shape = sin(1:50 * pi /50)
shape = shape - mean(shape)

set.seed(1)

sims = 500
trevor = matrix(0, length(par1), sims)
regina = matrix(0, length(par2), sims)

for(p in 1:length(par2)) {
  tic("Total")
  cat("amplitude shift: ", par2[p], "\n")
  for(s in 1:sims) {
    gp1 = gp1d(xn, mu = par1[p] * shape)
    gp2 = gp1d(yn, mu = par2[p] * shape)
    
    trevor[p,s] = coverage(gp1, gp2)
    regina[p,s] = quality(gp1, gp2)
    
  }
  toc()
  cat("\n")
}
shape_summary = data.frame(power = c(rowMeans(regina < 0.05),
                                     rowMeans(trevor < 0.05)),
                           method = c(rep("QI", length(par2)), 
                                      rep("KSD", length(par2))),
                           meanshift = rep(round(par2, 1), 2))

ggplot(shape_summary, aes(meanshift, power, color = method)) +
  geom_point() +
  geom_path() +
  theme_classic() +
  ylab("Power") +
  xlab("Shape amplitude") +
  ggtitle("Power against shape changes")




#### SHIFT MODES
par1 = c(seq(3, 0.5, length.out = 5), rep(0, 5))
par2 = c(rep(0, 5), seq(0.5, 3, length.out = 5))

shift = c(-par1[1:5], par2[6:10])

set.seed(1)

sims = 500
trevor = matrix(0, length(par1), sims)
regina = matrix(0, length(par2), sims)

for(p in 1:length(par2)) {
  tic("Total")
  cat("modal shift: ", shift[p], "\n")
  for(s in 1:sims) {
    gp1 = gp1d(xn, l = 10)
    gp2 = gp1d(yn, l = 10)
    
    gp1 = cbind(gp1[,1:(yn/2)] - par1[p]/2, gp1[,(yn/2 + 1):yn] + par1[p]/2)
    gp2 = cbind(gp2[,1:(yn/2)] - par2[p]/2, gp2[,(yn/2 + 1):yn] + par2[p]/2)
    
    trevor[p,s] = coverage(gp1, gp2)
    regina[p,s] = quality(gp1, gp2)
    
  }
  toc()
  cat("\n")
}
bimodal_summary = data.frame(power = c(rowMeans(regina < 0.05),
                                       rowMeans(trevor < 0.05)),
                             method = c(rep("QI", length(par2)), 
                                        rep("KSD", length(par2))),
                             meanshift = rep(round(shift, 1), 2))

ggplot(bimodal_summary, aes(meanshift, power, color = method)) +
  geom_point() +
  geom_path() +
  theme_classic() +
  ylab("Power") +
  xlab("Mode Separation") +
  ggtitle("Power against mode changes")
