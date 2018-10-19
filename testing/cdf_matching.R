### HERE I WAS TESING HOW BAD OF AN ASSUMPTION THE CURRENT INDEPEDENT MODEL IS
### ACTUALLY IT LOOKS LIKE THE MAX FOLLOWS THE SAME AS THE COMPONENTS....


rm(list = ls()); gc()
library(extdepth)
library(tictoc)
library(ggplot2)
library(e1071)

source("../research/assimilation-cfr/sim/reference.R")

coverage = function(f, g) {
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
  
  rate = sqrt((ncol(g)*ncol(f)) / (ncol(g) + ncol(f)))
  
  ksf = rate*max(abs(ffr - gfr))
  ksg = rate*max(abs(fgr - ggr))
  
  max(ksf, ksg)
}

# sims = 500
# 
# cover = rep(0, sims)
# brown = rep(0, sims)
# 
# for(s in 1:sims) {
#   tic("Total")
#   cat("Simulation ", s, "\n")
#   
#   cover[s] = coverage(gp1d(), gp1d())
#   brown[s] = max(abs(rbridge()))^2
#   
#   toc()
#   cat("\n")
# }
# 
# dens = data.frame(freq = c(cover, brown),
#                   method = c(rep("cover", sims),
#                              rep("brown", sims)))
# 
# ggplot(dens, aes(x = freq, color = method, fill = method)) + 
#   geom_density(alpha = 0.2) +
#   theme_classic()


ks_cdf = function(x, n = 20) {
  if(x < 0.1) return(0)
  1 - 2*(sum(sapply(1:n, function(k) ((-1)^(k-1)) * exp(-2*(k^2)*(x^2)))))
}



sims = 1000
cover = rep(0, sims)

for(s in 1:sims) {
  tic("Total")
  cat("Simulation ", s, "\n")
  
  cover[s] = coverage(gp1d(750), gp1d(750))
  
  toc()
  cat("\n")
}

t = seq(0, 2, length.out = 1000)
cov_cdf = sapply(t, function(x) mean(cover <= x))
bro_cdf = sapply(t, function(x) (ks_cdf(x)))
plot(t, cov_cdf, type = "l")
lines(t, bro_cdf, col = "red")




dens = data.frame(freq = c(cover, brown),
                  method = c(rep("cover", sims),
                             rep("brown", sims)))

ggplot(dens, aes(x = freq, color = method, fill = method)) + 
  geom_density(alpha = 0.2) +
  theme_classic()

