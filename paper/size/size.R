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

kvals = matrix(0, sims, 2)
qvals = matrix(0, sims, 2)

for(s in 1:sims) {
  tic("Total")
  cat("Simulation ", s, "\n")
  
  f = gp2d(fields = n, pts = pts, l = l)
  g = gp2d(fields = n+1, pts = pts, l = l)

  f = flatten(f)
  g = flatten(g)

  kvals[s,] = kolm(f, g)
  qvals[s,] = quality(f, g)
  
  toc()
}

size = data.frame(stat = c(kvals[,1], qvals[,1]),
                  pval = c(kvals[,2], qvals[,2]),
		  method = rep(c("K", "Q"), each=sims),
                  sims = sims,
                  pts = pts,
                  functions = n,
                  range = l,
                  seed = seed)

saveRDS(size, file = paste0("/home/trevorh2/size/independent/sim",i,".RDS"))

