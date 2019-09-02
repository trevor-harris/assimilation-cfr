rm(list = ls()); gc()

library(tictoc)
library(future.apply)

# FAD test
library(refund)

# Band test
library(roahd)

# set to the top level folder
setwd("/Users/trevorh2/research/assimilation-cfr/submit/")

source("method/depth_tests.R")
source("method/depths.R")
source("method/simulation.R")

# reproducibility
seed=042696

# shared parameters
n1=100
n2=50
pts=20

# X's parameters
mu1=0
sd1=1
r1=0.4
nu1=1.0

i=1
sims = 2000

# for constructing the mean and sd functions
t = seq(0, 1, length.out = pts)

# Mean Power
feature = "mean"
sd2 = sd1
r2 = r1
nu2 = nu1
for (amp in seq(0, 1, by = 0.05)) {
  set.seed(042696 + i)
  
  mu2 = 0.5 * amp * sin(t * pi * 4 - pi/2) + 1
  mu2 = as.vector(outer(mu2, mu2)) - 1
  
  vals = future_sapply(1:sims, function(x) {
    f = gp2d.mat(fields = n1, mu = mu1, sd = sd1, range = r1, nu = nu1, pts = pts)
    g = gp2d.mat(fields = n2, mu = mu2, sd = sd2, range = r2, nu = nu2, pts = pts)
    
    f = flatten(f)
    g = flatten(g)
    
    c(kolm(f, g)[2], quality(f, g)[2], fadtest(f, g), bandtest(f, g)[2])
    
  })
  
  power = data.frame(pval = c(vals[1,], vals[2,], vals[3,], vals[4,]),
                     method = rep(c("Kolm", "Qual", "FAD", "MBD"), each=sims),
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
  
  saveRDS(power, file = paste0("power/out_het_gp/",feature,"_sim_",i,"_const",".RDS"))
  
  i = i + 1
}

# Standard Deviation Power
feature = "sd"
mu2 = mu1
r2 = r1
nu2 = nu1
for (amp in seq(0, 1, by = 0.05)) {
  set.seed(042696 + i)
  
  sd2 = 0.5 * amp * sin(t * pi * 4 - pi/2) + 1
  sd2 = as.vector(outer(sd2, sd2))
  
  vals = future_sapply(1:sims, function(x) {
    f = gp2d.mat(fields = n1, mu = mu1, sd = sd1, range = r1, nu = nu1, pts = pts)
    g = gp2d.mat(fields = n2, mu = mu2, sd = sd2, range = r2, nu = nu2, pts = pts)
    
    f = flatten(f)
    g = flatten(g)
    
    c(kolm(f, g)[2], quality(f, g)[2], fadtest(f, g), bandtest(f, g)[2])
    
  })
  
  power = data.frame(pval = c(vals[1,], vals[2,], vals[3,], vals[4,]),
                     method = rep(c("Kolm", "Qual", "FAD", "MBD"), each=sims),
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
  
  saveRDS(power, file = paste0("power/out_het_gp/",feature,"_sim_",i,"_const",".RDS"))
  
  i = i + 1
}

