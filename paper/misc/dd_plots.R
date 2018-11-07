rm(list = ls()); gc()

library(ggplot2)
library(reshape2)
library(extdepth)
source('../assimilation-cfr/sim/reference.R')
# source("../research/assimilation-cfr/sim/reference.R")

xdepth = function(gmat, fmat) {
  sapply(gmat, function(x) mean(depth(x, fmat)))
}

gmat = gp1d()
fmat = gp1d()
proj = 10

pdepth = function(gmat, fmat, proj = 5000) {
  u = matrix(runif(nrow(gmat) * proj, -1, 1), nrow(gmat), proj)
  u = apply(u, 2, function(x) x / sqrt(sum(x^2)))
  
  Fu = t(u) %*% fmat
  Gu = t(u) %*% gmat
  
  Fu.mu = apply(Fu, 1, median)
  Fu.sig = apply(Fu, 1, mad)
  
  out = t(abs(t(Gu) - Fu.mu) / Fu.sig)

  return(1 / (1 + apply(out, 2, max)))
}

set.seed(1023)
mu = 0
std = 1
rng = 50
pts = 20

f = gp1d(pts = pts)
g = gp1d(mu = mu, sd = std, l = rng, pts = pts)

# ## Expected Depth
# ffxd = xdepth(f, f)
# fgxd = xdepth(f, g)
# 
# ggxd = xdepth(g, g)
# gfxd = xdepth(g, f)
# 
# ## Extremal Depths
# ffxd = edepth_set(f)
# fgxd = apply(f, 2, function(x) edepth(x, g))
# 
# ggxd = edepth_set(g)
# gfxd = apply(g, 2, function(x) edepth(x, f))

set.seed(1023)
pts = 20
infil = 10
f = gp1d(100, sd = 1, pts = pts, l = 10)
g = gp1d(100, sd = 1, pts = pts, l = 10)

f = apply(f, 2, function(x) smooth.spline(x, nknots = 10)$y)
g = apply(g, 2, function(x) smooth.spline(x, nknots = 10)$y)

f = apply(f, 2, function(x) spline(1:pts, x, n = infil*pts)$y)
g = apply(g, 2, function(x) spline(1:pts, x, n = infil*pts)$y)

## Projection Depths
ffxd = xdepth(f, f)
fgxd = xdepth(f, g)

ggxd = xdepth(g, g)
gfxd = xdepth(g, f)

f.dd = data.frame(Reference = ffxd, Proposed = fgxd)
g.dd = data.frame(Reference = ggxd, Proposed = gfxd)
dd = rbind(f.dd, g.dd)
dd["Group"] = rep(c("F", "G"), each=50)

# ggplot(dd, aes(Reference, Proposed, color = Group)) +
#   geom_abline(intercept = 0, slope = 1, alpha = 0.7) +
#   geom_point() +
#   theme_classic() +
#   ylim(c(0, 1)) +
#   xlim(c(0, 1)) +
#   ylab("Proposed Depths") +
#   xlab("Reference Depths") +
#   ggtitle(paste0("DD-plot with G ~ GP(mu=", mu, ", sd=", std, ", l=", rng, ")"))
# 
# 
# mean(fgxd)
# mean(gfxd)
# 
# sd(fgxd)
# sd(gfxd)

ggplot(dd, aes(Reference, Proposed, color = Group)) +
  geom_abline(intercept = 0, slope = 1, alpha = 0.7) +
  geom_point() +
  geom_rug(alpha = 0.5) +
  theme_classic() +
  ylim(c(0, 1)) +
  xlim(c(0, 1)) +
  ylab("Proposed Depths") +
  xlab("Reference Depths") +
  ggtitle(paste0("DD-plot with G ~ GP(mu=", mu, ", sd=", std, ", l=", rng, ")"))







