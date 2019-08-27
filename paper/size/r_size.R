library(reshape2)
library(tictoc)

source("/home/trevorh2/assimilation-cfr/code/depth_tests.R")
source("/home/trevorh2/assimilation-cfr/code/depths.R")
source("/home/trevorh2/assimilation-cfr/code/simulation.R")

# fixed settings
iter = 50
sims = 1000

# get args
args = commandArgs(TRUE)
n1 = as.double(args[1])
n2 = as.double(args[2])
rng = as.double(args[3])
nu = as.double(args[4])
i = as.integer(args[5])
seed = as.integer(args[6]) + i

set.seed(seed)

# test settings
# iter = 2
# sims = 5
# n1 = 100
# n2 = 100
# rng = 0.3
# nu = 0.8

# initiate output files
kvals = matrix(NA, sims*iter, 3)
qvals = matrix(NA, sims*iter, 3)

j = 1
for(k in 1:iter) {
  tic(paste0("Iter ", k))
  
  for(s in 1:sims) {
    f = gp1d.mat(fields = n1, range = rng, nu = nu)
    g = gp1d.mat(fields = n2, range = rng, nu = nu)
    
    kvals[j,] = c(k, s, kolm(f, g)[2])
    qvals[j,] = c(k, s, quality(f, g)[2])
    
    j = j + 1
  }
  
  toc()
}

ksize = data.frame(iter = kvals[,1],
                   sim = kvals[,2],
                   n1 = n1,
                   n2 = n2,
                   range = rng,
                   nu = nu,
                   stat = "K",
                   size = kvals[,3])

qsize = data.frame(iter = qvals[,1],
                   sim = qvals[,2],
                   n1 = n1,
                   n2 = n2,
                   range = rng,
                   nu = nu,
                   stat = "Q",
                   size = qvals[,3])

size = rbind(ksize, qsize)
saveRDS(size, file = paste0("/home/trevorh2/kstat/size/out/size",i,".RDS"))

