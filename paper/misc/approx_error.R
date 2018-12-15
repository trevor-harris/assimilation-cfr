library(tictoc)
library(reshape2)
library(ggplot2)

source("../assimilation-cfr/code/simulation.R")
source("../assimilation-cfr/code/depths.R")

# my way
reg_way = function(x, y) {
  dxx = xdepth(x, x)
  dyx = xdepth(y, x)

  cdf.xx = sapply(dxx, function(y) mean(dxx <= y))
  cdf.yx = sapply(dxx, function(y) mean(dyx <= y))

  max(cdf.xx - cdf.yx)
  sapply(dxx, function(y) mean(dxx <= y))
}

# split way
split_way = function(x, y, z) {
  dzz = xdepth(z, z)
  dxz = xdepth(x, z)
  # dyz = xdepth(y, z)
  
  # cdf.xz = sapply(dzz, function(y) mean(dxz <= y))
  # cdf.yz = sapply(dzz, function(y) mean(dyz <= y))
  
  # max(cdf.xz - cdf.yz)
  
  sapply(dzz, function(y) mean(dxz <= y))
}


# just one N
n = 1000
l = 40

sims = 50
diff = rep(0, sims)
for(i in 1:sims) {
  tic(paste0("Sim: ", i))
  
  x = gp1d(n, l = l)
  y = gp1d(n, l = l)
  z = gp1d(500, l = l)
  
  diff[i] = reg_way(x, y) - split_way(x, y, z)
  
  toc()
}

plot(diff)


# several  N
set.seed(1023)

N = c(100, 250, 500, 750, 1000)
l = 40

ns = length(N)
sims = 100
diff = matrix(0, sims, ns)

for(n in 1:ns) {
  for(i in 1:sims) {
    tic(paste0("Sim: ", i))
    
    x = gp1d(N[n], l = l)
    y = gp1d(N[n], l = l)
    z = gp1d(2000, l = l)
    
    diff[i,n] = reg_way(x, y) - split_way(x, y, z)
    
    toc()
  }
}

diffs = melt(diff)
diffs[["N"]] = as.factor(diffs[["Var2"]])
levels(diffs[["N"]]) = N

ggplot(diffs, aes(x = N, y = value, group = N)) + 
  geom_boxplot() +
  geom_hline(yintercept = 0) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5)) +
  ggtitle("Difference in K when splitting v.s. not-splitting")

save_dir = "../research/assimilation-cfr/paper/misc/"
ggsave(paste0(save_dir,"approx_error.png"), width = 5, height = 3.2)






# several  N
set.seed(1023)

N = c(100, 500, 1000, 1500, 2000, 3000, 5000)
l = 40

ns = length(N)
sims = 100
diff = matrix(0, sims, ns)

for(n in 1:ns) {
  for(i in 1:sims) {
    tic(paste0("Sim: ", i))
    
    x = gp1d(N[n], l = l)
    
    g1 = sample(1:N[n], N[n]/2)
    
    x1 = x[, g1]
    x2 = x[, -g1]
    
    xd.split = c(xdepth(x1, x2), xdepth(x2, x1))
    
    diff[i,n] = mean(xdepth(x, x) <= 0.5) - mean(xd.split <= 0.5)
    
    toc()
  }
}


diff2 = sapply(1:length(N), function(x) diff[,x] * sqrt(N[x]))

diffs = melt(diff)
diffs[["N"]] = as.factor(diffs[["Var2"]])
levels(diffs[["N"]]) = N

ggplot(diffs, aes(x = N, y = value, group = N)) + 
  geom_boxplot() +
  geom_hline(yintercept = 0) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5)) +
  ggtitle("mean(xdepth(x, x) <= 0.5) - mean(xd.split <= 0.5)")



boot = 20
x = gp1d(100, l = l)
x.boot = matrix(0, boot, 100)

for(b in 1:boot) {
  g1 = sample(1:100, 50)
  x1 = x[, g1]
  x2 = x[, -g1]
  
  x.boot[b, g1] = xdepth(x1, x2)
  x.boot[b, -g1] = xdepth(x2, x1)
}

x.boot = colMeans(x.boot)
x.depth = xdepth(x, x)
x.split = c(xdepth(x1, x2), xdepth(x2, x1))

xd.split = c(xdepth(x1, x2), xdepth(x2, x1))




## use the bootstrapped depths
set.seed(1023)

N = c(100, 500, 1000, 1500)
l = 40

ns = length(N)
sims = 100
diff = matrix(0, sims, ns)

for(n in 1:ns) {
  for(i in 1:sims) {
    tic(paste0("Sim: ", i))
    
    x = gp1d(N[n], l = l)
    
    boot = 50
    x.boot = matrix(0, boot, N[n])
    
    for(b in 1:boot) {
      g1 = sample(1:N[n], N[n]/2)
      x1 = x[, g1]
      x2 = x[, -g1]
      
      x.boot[b, g1] = xdepth(x1, x2)
      x.boot[b, -g1] = xdepth(x2, x1)
    }
    x.boot = colMeans(x.boot)
    
    diff[i,n] = mean(xdepth(x, x) <= 0.5) - mean(x.boot <= 0.5)
    
    toc()
  }
}

diff2 = sapply(1:length(N), function(x) diff[,x] * sqrt(N[x]))

diffs = melt(diff2)
diffs[["N"]] = as.factor(diffs[["Var2"]])
levels(diffs[["N"]]) = N

ggplot(diffs, aes(x = N, y = value, group = N)) + 
  geom_boxplot() +
  geom_hline(yintercept = 0) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5)) +
  ggtitle("mean(xdepth(x, x) <= 0.5) - mean(xd.split <= 0.5)")





## use the bootstrapped depths
set.seed(0211)

N = c(100, 400, 700, 1000)
l = 40

ns = length(N)
sims = 100
diff = matrix(0, sims, ns)

for(n in 1:ns) {
  for(i in 1:sims) {
    tic(paste0("Sim: ", i))
    
    x = gp1d(N[n], l = l, pts = 25)
    
    boot = 100
    x.boot = matrix(0, boot, N[n])
    
    for(b in 1:boot) {
      g1 = sample(1:N[n], N[n]/2)
      x1 = x[, g1]
      x2 = x[, -g1]
      
      x.boot[b, g1] = xdepth(x1, x2)
      x.boot[b, -g1] = xdepth(x2, x1)
    }
    x.boot = colMeans(x.boot)
    
    diff[i,n] = mean(xdepth(x, x) <= 0.5) - mean(x.boot <= 0.5)
    
    toc()
  }
}

diff2 = sapply(1:length(N), function(x) diff[,x] * sqrt(N[x]))

diffs = melt(diff2)
diffs[["N"]] = as.factor(diffs[["Var2"]])
levels(diffs[["N"]]) = N

ggplot(diffs, aes(x = N, y = value, group = N)) + 
  geom_boxplot() +
  geom_hline(yintercept = 0) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5)) +
  ggtitle("mean(xdepth(x, x) <= 0.5) - mean(xd.split <= 0.5)")



## splitting method
N = 100
x = gp1d(N, l = l, pts = 25)

boot = 50
x.boot = matrix(0, boot, N)

for(b in 1:boot) {
  g1 = sample(1:N, N/2)
  x1 = x[, g1]
  x2 = x[, -g1]
  
  x.boot[b, g1] = xdepth(x1, x2)
  x.boot[b, -g1] = xdepth(x2, x1)
}

x.boot = colMeans(x.boot)
x.depth = xdepth(x, x)

t = seq(0, 1, length.out = 500)
x.boot.cdf = sapply(t, function(y) mean(x.boot <= y))
x.depth.cdf = sapply(t, function(y) mean(x.depth <= y))

mean((x.boot.cdf - x.depth.cdf))

plot(x.boot.cdf, type = "l")
lines(x.depth.cdf, col = "red")
