rm(list = ls()); gc()
library(tictoc)

# get args
args = commandArgs(TRUE)
n = as.double(args[1])
pts = as.integer(args[2])
mu1 = as.double(args[3])
mu2 = as.double(args[4])
sd1 = as.double(args[5])
sd2 = as.double(args[6])
r1 = as.integer(args[7])
r2 = as.integer(args[8])
i = as.integer(args[9])
seed = as.integer(args[10]) + i
feature = as.character(args[11])
sims = as.integer(args[12])

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
  
  f = gp1d(fields = n, mu = mu1, sd = sd1, l = r1, pts = pts)
  g = gp1d(fields = n+1, mu = mu2, sd = sd2, l = r2, pts = pts)

  kvals[s,] = kolm(f, g)
  qvals[s,] = quality(f, g)
  
  toc()
}

power = data.frame(stat = c(kvals[,1], qvals[,1]),
                  pval = c(kvals[,2], qvals[,2]),
		  method = rep(c("K", "Q"), each=sims),
                  mu1 = mu1,
		  mu2 = mu2,
                  sd1 = sd1,
                  sd2 = sd2,
                  r1 = r1,
                  r2 = r2,
                  feature = feature,
                  sims = sims,
                  pts = pts,
                  functions = n,
                  seed = seed)

saveRDS(power, file = paste0("/home/trevorh2/power/independent1d/",feature,"sim",i,".RDS"))
