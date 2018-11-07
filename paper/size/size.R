rm(list = ls()); gc()
library(tictoc)

# get args
args = commandArgs(TRUE)
n = as.double(args[1])
d = as.integer(args[2])
l = as.integer(args[3])
i = as.integer(args[4])
seed = as.integer(args[5]) + i
pts = as.integer(args[6])
sims = as.integer(args[7])

source("/home/trevorh2/assimilation-cfr/code/depth_tests.R")
source("/home/trevorh2/assimilation-cfr/code/depths.R")
source("/home/trevorh2/assimilation-cfr/code/simulation.R")

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
  g = gp1d(fields = n+1, pts = pts, l = l)
  
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
		  method = rep(c("K_XD", "K_PD", "Q_XD", "Q_PD"), each=500),
                  sims = sims,
                  pts = pts,
                  functions = n,
                  range = l,
                  seed = seed)

saveRDS(size, file = paste0("/home/trevorh2/size/output/sim",i,".RDS"))

