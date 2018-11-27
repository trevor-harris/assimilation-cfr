library(tictoc)
library(reshape2)
library(ggplot2)

source("../research/assimilation-cfr/code/simulation.R")
source("../research/assimilation-cfr/code/depths.R")


# my way
reg_way = function(x, y) {
  dxx = xdepth(x, x)
  dyx = xdepth(y, x)
  
  cdf.xx = sapply(dxx, function(y) mean(dxx <= y))
  cdf.yx = sapply(dxx, function(y) mean(dyx <= y))
  
  max(cdf.xx - cdf.yx)
}

# split way
split_way = function(x, y, z) {
  dzz = xdepth(z, z)
  dxz = xdepth(x, z)
  dyz = xdepth(y, z)
  
  cdf.xz = sapply(dzz, function(y) mean(dxz <= y))
  cdf.yz = sapply(dzz, function(y) mean(dyz <= y))
  
  max(cdf.xz - cdf.yz)
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
  z = gp1d(n, l = l)
  
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
    z = gp1d(N[n], l = l)
    
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

