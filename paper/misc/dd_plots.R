rm(list = ls()); gc()

library(ggplot2)
library(reshape2)
library(extdepth)
source('../assimilation-cfr/sim/reference.R')

set.seed(0426)
mu = 2
std = 2
rng = 2
pts = 100

f = gp1d(pts = pts)
g = gp1d(mu = mu, sd = std, l = rng, pts = pts)

ffxd = xdepth(f, f)
fgxd = xdepth(f, g)

ggxd = xdepth(g, g)
gfxd = xdepth(g, f)

f.dd = data.frame(Reference = ffxd, Proposed = fgxd)
g.dd = data.frame(Reference = ggxd, Proposed = gfxd)
dd = rbind(f.dd, g.dd)
dd["Group"] = rep(c("X", "Y"), each=100)

ggplot(dd, aes(Reference, Proposed, color = Group)) +
  geom_abline(intercept = 0, slope = 1, alpha = 0.7) +
  geom_point() +
  theme_classic() +
  ylim(c(0, 1)) +
  xlim(c(0, 1)) +
  ylab("Proposed Depths") +
  xlab("Reference Depths") +
  ggtitle(paste0("DD-plot with G ~ GP(mu=", mu, ", sd=", std, ", l=", rng, ")"))


# ggplot(dd, aes(Reference, Proposed, color = Group)) +
#   geom_abline(intercept = 0, slope = 1, alpha = 0.7) +
#   geom_point() +
#   geom_rug(alpha = 0.5) +
#   theme_classic() +
#   ylim(c(0, 1)) +
#   xlim(c(0, 1)) +
#   ylab("Proposed Depths") +
#   xlab("Reference Depths") +
#   ggtitle(paste0("DD-plot with G ~ GP(mu=", mu, ", sd=", std, ", l=", rng, ")"))







