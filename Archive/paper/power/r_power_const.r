rm(list = ls()); gc()
library(tictoc)
library(future.apply)

# get args
args = commandArgs(TRUE)
n1 = as.double(args[1])
n2 = as.double(args[2])
pts = as.integer(args[3])
mu1 = as.double(args[4])
mu2 = as.double(args[5])
sd1 = as.double(args[6])
sd2 = as.double(args[7])
r1 = as.double(args[8])
r2 = as.double(args[9])
nu1 = as.double(args[10])
nu2 = as.double(args[11])
i = as.integer(args[12])
seed = as.integer(args[13]) + i
feature = as.character(args[14])
sims = as.integer(args[15])

source("/home/trevorh2/assimilation-cfr/code/depth_tests.R")
source("/home/trevorh2/assimilation-cfr/code/depths.R")
source("/home/trevorh2/assimilation-cfr/code/simulation.R")

#### SIZE
set.seed(seed)

plan(multicore)

vals = future_sapply(1:sims, function(x) {
  f = gp2d.mat(fields = n1, mu = mu1, sd = sd1, range = r1, nu = nu1, pts = pts)
  g = gp2d.mat(fields = n2, mu = mu2, sd = sd2, range = r2, nu = nu2, pts = pts)
  
  f = flatten(f)
  g = flatten(g)
  
  c(kolm(f, g)[2], quality(f, g)[2])
})

power = data.frame(pval = c(vals[1,], vals[2,]),
                   method = rep(c("K", "Q"), each=sims),
                   mu1 = mu1,
                   mu2 = mu2,
                   sd1 = sd1,
                   sd2 = sd2,
                   r1 = r1,
                   r2 = r2,
                   nu1 = nu1,
                   nu2 = nu2,
                   n1 = n1,
                   n2 = n2,
                   feature = feature)

saveRDS(power, file = paste0("/home/trevorh2/matern/power/const/",feature,"sim_",i,"_const",".RDS"))

