
rm(list = ls()); gc()

library(reshape2)
library(tictoc)
library(future)
library(future.apply)

source("../research/assimilation-cfr/code/simulation.R")
source("../research/assimilation-cfr/code/depths.R")
source("../research/assimilation-cfr/code/depth_tests.R")


args = commandArgs(TRUE)
std = 1
l = 2
n = 50
seed = 123
i = 1
sims = 10

#### SIZE
set.seed(seed)

kvals = matrix(0, sims, 2)
qvals = matrix(0, sims, 2)

for(s in 1:sims) {
  tic("Total")
  cat("Simulation ", s, "\n")
  
  f = gp1d(fields = n, sd  = std, l = l)
  g = gp1d(fields = n+1, sd  = std, l = l)
  
  
  kvals[s,] = kolm(f, g)
  qvals[s,] = quality(f, g)
  
  toc()
}


n = 50
iter = 10
sims = 5
ranges = c(10, 20, 30, 40)
stdevs = c(0.1, 0.5, 1, 2, 3)

ksize = data.frame(size = numeric(iter*length(ranges)*length(stdevs)),
                   range = numeric(iter*length(ranges)*length(stdevs)),
                   std = numeric(iter*length(ranges)*length(stdevs)))

qsize = data.frame(size = numeric(0),
                   range = numeric(0),
                   std = numeric(0))

k = 1
for(r in ranges) {
  for (s in stdevs) {
    tic(paste0("r = ", r, " s = ", s))
    
    klevel = rep(0, iter)
    qlevel = rep(0, iter)
    for(j in 1:iter) {
      
      kpval = rep(0, sims)
      qpval = rep(0, sims)
      for (i in 1:sims) {
        f = gp1d(n, sd = s, l = r)
        g = gp1d(n+1, sd = s, l = r)
        
        kpval[i] = kolm(f, g)[2]
        qpval[i] = quality(f, g)[2]
      }
      
      ksize[k,] = c(mean(kpval <= 0.05), r, s)
      qsize[k,] = c(mean(qpval <= 0.05), r, s)
      
      k = k + 1
    }
    
    toc()
  }
}


##### multicore version

n = 50
iter = 10
sims = 500
ranges = c(10, 20, 30, 40)
stdevs = c(0.1, 0.5, 1, 2, 3)
nfuns = c(50, 100, 200, 300, 400, 500)

ksize = data.frame(size = numeric(0),
                    range = numeric(0),
                    std = numeric(0))

qsize = data.frame(size = numeric(0),
                   range = numeric(0),
                   std = numeric(0))

plan(multiprocess)

k = 1

for(n in nfuns) {
  for(r in ranges) {
    for (s in stdevs) {
      tic(paste0("n = ", n, " r = ", r, " s = ", s))
      
      klevel = rep(0, iter)
      qlevel = rep(0, iter)
      for(j in 1:iter) {
        
        kpval = rep(0, sims)
        qpval = rep(0, sims)
        
        val = future_sapply(1:sims, function(x) {
          f = gp1d(n, sd = s, l = r)
          g = gp1d(n+1, sd = s, l = r)
          
          c(kolm(f, g)[2], quality(f, g)[2], 3)
        })
        
        ksize[k,] = c(mean(val[1,] <= 0.05), r, s)
        qsize[k,] = c(mean(val[2,] <= 0.05), r, s)
        
        k = k + 1
      }
      
      toc()
    }
  }
}


