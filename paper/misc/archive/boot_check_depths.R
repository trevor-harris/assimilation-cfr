library(tictoc)
library(reshape2)
library(ggplot2)
library(future)
library(future.apply)

source("../research/assimilation-cfr/code/simulation.R")
source("../research/assimilation-cfr/code/depths.R")

n = 500
x = gp1d(n, l = l, pts = 25)
boot = seq(20, 20000, length.out = 100)
x.boot = matrix(0, n, length(boot))
x.biased = matrix(0, n, length(boot))

for(b in 1:length(boot)) {

  tic("Boot")  
  x.b = future_sapply(1:b, function(t) {
    x.depth = rep(0, n)
    
    group_ind = sample(1:n, n/2)
    x1 = x[, group_ind]
    x2 = x[, -group_ind]
    
    x.depth[group_ind] = xdepth(x1, x2)
    x.depth[-group_ind] = xdepth(x2, x1)
    
    x.depth
  })
  
  x.biased[,b] = rowMeans(x.b)
  x.boot[,b] = rowMeans(x.b) + 1/n
  
  toc()
}

x.depth = xdepth(x, x)
x.diff1 = x.biased[,10] - x.depth

mean(x.diff1)
1/500

plot(x.diff1)
abline(mean(x.diff1), 0)
abline(0, 0, col = "red")




x.diff = apply(x.boot, 2, function(y) mean((y - x.depth)^2))
plot(boot, x.diff, type = "l", ylim = c(0, max(x.diff)))
abline(0, 0, col = "blue")

d = 1
x.diff = sapply(x.boot[d,], function(y) (y - x.depth[d]))
plot(boot, x.diff, type = "l", ylim = c(min(x.diff), max(x.diff)))
abline(mean(x.diff), 0, col = "red")
abline(0, 0, col = "blue")

#seems to be a bias term
d = 90
d.bias = sapply(1:n, function(d) mean(sapply(x.boot[d,], function(y) (y - x.depth[d]))))
plot(1:n, d.bias)
abline(mean(d.bias), 0)
abline(0, 0, col = "red")
mean(d.bias)
1/n
