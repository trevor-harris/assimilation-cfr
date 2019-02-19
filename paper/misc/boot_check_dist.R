library(tictoc)
library(reshape2)
library(ggplot2)
library(future)
library(future.apply)

source("../research/assimilation-cfr/code/simulation.R")
source("../research/assimilation-cfr/code/depths.R")


plan(multiprocess)
# several  N
N = c(100, 400, 700, 1000, 1300)
l = 40

ns = length(N)
sims = 100
diff = matrix(0, sims, ns)

tic("Total")
for(n in 1:ns) {
  cat(paste0("Number of function - ", N[n], "\n"))
  
  for(i in 1:sims) {
    tic(paste0("Sim - ", i))
    
    x = gp1d(N[n], l = l, pts = 25)
    
    boot = 512
    x.boot = future_sapply(1:boot, function(b) {
      x.depth = rep(0, N[n])
      
      group_ind = sample(1:N[n], N[n]/2)
      x1 = x[, group_ind]
      x2 = x[, -group_ind]
      
      x.depth[group_ind] = xdepth(x1, x2)
      x.depth[-group_ind] = xdepth(x2, x1)
      
      x.depth
    })
    x.boot = rowMeans(x.boot) + 1/N[n]
    x.depth = xdepth(x, x)
    
    diff[i,n] = mean(xdepth(x, x) <= 0.5) - mean(x.boot <= 0.5)
    
    toc()
  }
}
toc()

diff2_0.5 = abs(sapply(1:ns, function(x) diff[,x] * sqrt(N[x])))

diffs = melt(diff2_0.5)
diffs[["N"]] = as.factor(diffs[["Var2"]])
levels(diffs[["N"]]) = N

ggplot(diffs, aes(x = N, y = value, group = N)) + 
  geom_boxplot() +
  geom_hline(yintercept = 0) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5)) +
  ggtitle("sqrt(N) * (mean(xdepth(x, x) <= 0.5) - mean(x.boot <= 0.5))")
