
rm(list = ls()); gc()

library(reshape2)
library(tictoc)
library(future)
library(future.apply)

source("/home/trevorh2/assimilation-cfr/code/depth_tests.R")
source("/home/trevorh2/assimilation-cfr/code/depths.R")
source("/home/trevorh2/assimilation-cfr/code/simulation.R")

#### SIZE
set.seed(042696)


iter = 10
sims = 500
ranges = c(10, 20, 30, 40)
stdevs = c(0.1, 0.5, 1, 2, 3)
nfuns = c(50, 100, 200, 300, 400, 500)

ksize = data.frame(size = numeric(0),
                    nfun = numeric(0),
                    range = numeric(0),
                    std = numeric(0))

qsize = data.frame(size = numeric(0),
                   nfun = numeric(0),
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
          
          c(kolm(f, g)[2], quality(f, g)[2])
        })
        
        ksize[k,] = c(mean(val[1,] <= 0.05), n, r, s)
        qsize[k,] = c(mean(val[2,] <= 0.05), n, r, s)
        
        k = k + 1
      }
      
      toc()
    }
  }
}


saveRDS(ksize, file = "/home/trevorh2/assimilation-cfr/sim/size2/out/ksize.RDS")
saveRDS(qsize, file = "/home/trevorh2/assimilation-cfr/sim/size2/out/qsize.RDS")
