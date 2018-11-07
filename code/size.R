rm(list = ls()); gc()
library(tictoc)

# get args
args = commandArgs(TRUE)
n = as.double(args[1])
d = as.integer(args[1])
l = as.integer(args[2])
i = as.integer(args[3])
pts = as.integer(args[4])
seed = as.integer(args[5]) + i + 1000
sims = as.integer(args[6])

n = 500
d = 30
l = 30
pts = 10
sims = 10
seed = 1

source("../research/assimilation-cfr/code/depth_tests.R")
source("../research/assimilation-cfr/code/depths.R")
source("../research/assimilation-cfr/code/simulation.R")

#### SIZE
set.seed(seed)

k.xd = matrix(0, sims, 2)
k.pd = matrix(0, sims, 2)

q.xd = matrix(0, sims, 2)
q.pd = matrix(0, sims, 2)

for(s in 1:sims) {
  tic("Total")
  cat("Simulation ", s, "\n")
  
  f = gp1d(fields = n, pts = pts, l = l)
  g = gp1d(fields = n, pts = pts, l = l)
  
  f = apply(f, 2, function(x) spline(1:pts, x, n = d*pts)$y)
  g = apply(g, 2, function(x) spline(1:pts, x, n = d*pts)$y)
  
  k.xd[s,] = kolm.xd(f, g)
  k.pd[s,] = kolm.pd(f, g)
  
  q.xd[s,] = quality.xd(f, g)
  q.pd[s,] = quality.pd(f, g)
  
  toc()
}

size = data.frame(stat = c(k.xd[,1], k.pd[,1], q.xd[,1], q.pd[,1]),
                  pval = c(k.xd[,2], k.pd[,2], q.xd[,2], q.pd[,2]),
                  sims = sims,
                  pts = pts,
                  functions = n,
                  range = l,
                  seed = seed)

saveRDS(size, file = paste0("sim",i,".RDS"))

