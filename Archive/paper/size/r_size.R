rm(list = ls()); gc()
library(tictoc)
library(future.apply)

# get args
args = commandArgs(TRUE)
n1 = as.double(args[1])
n2 = as.double(args[2])
rng = as.double(args[3])
nu = as.double(args[4])
i = as.integer(args[5])
seed = as.integer(args[6]) + i
sims = as.integer(args[7])

source("/home/trevorh2/assimilation-cfr/code/depth_tests.R")
source("/home/trevorh2/assimilation-cfr/code/depths.R")
source("/home/trevorh2/assimilation-cfr/code/simulation.R")

#### SIZE
set.seed(seed)

plan(multicore)

vals = future_sapply(1:sims, function(x) {
  f = gp2d.mat(fields = n1, range = rng, nu = nu)
  g = gp2d.mat(fields = n2, range = rng, nu = nu)
  
  f = flatten(f)
  g = flatten(g)
  
  c(kolm(f, g)[2], quality(f, g)[2])
  
})

size = data.frame(pval = c(vals[1,], vals[2,]),
                  method = rep(c("K", "Q"), each=sims),
                  rng = rng,
                  nu = nu,
                  n1 = n1,
                  n2 = n2)

saveRDS(size, file = paste0("/home/trevorh2/matern/size/out/","size_sim_",i,".RDS"))

