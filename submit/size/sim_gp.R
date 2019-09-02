rm(list = ls()); gc()

library(tictoc)
library(future.apply)

# FAD test
library(refund)

# Band test
library(roahd)

# set to the top level folder
# set to the top level folder
setwd("/Users/trevorh2/research/assimilation-cfr/submit/")

source("method/depth_tests.R")
source("method/depths.R")
source("method/simulation.R")

#### SIZE
plan(multicore)

i=1
sims = 2000

# do not actually run this entire loop. The inner part was run in parallel on a cluster.
for (n1 in c(50,100,200,300)) {
  for (n2 in c(50,100,200,300)) {
    for (rng in c(0.2,0.3,0.4,0.5)) {
      for (nu in c(0.5,1.0,1.5)) {
        
        set.seed(072393 + i)
        
        vals = future_sapply(1:sims, function(x) {
          f = gp2d.mat(fields = n1, range = rng, nu = nu)
          g = gp2d.mat(fields = n2, range = rng, nu = nu)
          
          f = flatten(f)
          g = flatten(g)
          
          c(kolm(f, g)[2], quality(f, g)[2], fadtest(f, g), bandtest(f, g)[2])
          
        })
        
        size = data.frame(pval = c(vals[1,], vals[2,], vals[3,], vals[4,]),
                          method = rep(c("Kolm", "Qual", "FAD", "MBD"), each=sims),
                          rng = rng,
                          nu = nu,
                          n1 = n1,
                          n2 = n2)
        
        saveRDS(size, file = paste0("out_gp/","size_sim_",i,".RDS"))
        
        i = i + 1
      }
    }
  }
}

