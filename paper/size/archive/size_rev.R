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

k.wa = matrix(0, sims, 2)
q.xd = matrix(0, sims, 2)

for(s in 1:sims) {
  tic("Total")
  cat("Simulation ", s, "\n")
  
  f = gp1d(fields = n, pts = pts, l = l)
  g = gp1d(fields = n+1, pts = pts, l = l)

  f = apply(f, 2, function(x) spline(1:pts, x, n = d*pts)$y)
  g = apply(g, 2, function(x) spline(1:pts, x, n = d*pts)$y)
  
  # f = apply(f, 2, function(x) smooth.spline(x, nknots = 20)$y)
  # g = apply(g, 2, function(x) smooth.spline(x, nknots = 20)$y)
 
  k.wa[s,] = kolm.wa(f, g)
  q.xd[s,] = quality.xd(f, g)
  
  toc()
}

size = data.frame(stat = c(k.wa[,1], q.xd[,1]),
                  pval = c(k.wa[,2], q.xd[,2]),
		  method = rep(c("K_WA", "Q_XD"), each=sims),
                  sims = sims,
                  pts = pts,
                  functions = n,
                  range = l,
                  seed = seed)

saveRDS(size, file = paste0("/home/trevorh2/size/output_rev/sim",i,".RDS"))

