rm(list = ls()); gc()
library(tictoc)
library(future.apply)
library(fdasrvf)

# FAD
library(refund)

# BAND
library(roahd)

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
amp = as.double(args[12])
i = as.integer(args[13])
seed = as.integer(args[14]) + i
feature = as.character(args[15])
sims = as.integer(args[16])

source("/home/trevorh2/assimilation-cfr/code/depth_tests.R")
source("/home/trevorh2/assimilation-cfr/code/depths.R")
source("/home/trevorh2/assimilation-cfr/code/simulation.R")


#### SIZE
set.seed(seed)

plan(multicore)

t = seq(0, 1, length.out = pts)
if(feature == "mean") {
  mu2 = 0.5 * amp * sin(t * pi * 4 - pi/2) + 1
  mu2 = as.vector(outer(mu2, mu2)) - 1
}

if(feature == "sd") {
  sd2 = 0.5 * amp * sin(t * pi * 4 - pi/2) + 1
  sd2 = as.vector(outer(sd2, sd2))
}

vals = future_sapply(1:sims, function(x) {
 # f = abs(gp2d.mat(fields = n1, mu = mu1, sd = sd1, range = r1, nu = nu1, pts = pts))
 # g = abs(gp2d.mat(fields = n2, mu = mu2, sd = sd2, range = r2, nu = nu2, pts = pts))
  
 # f = flatten(f)
 # g = flatten(g)

  f = tp2d(fields = n1, mu = mu1, df = 3, phi = sd1, range = r1, nu = nu1, pts = pts)
  g = tp2d(fields = n2, mu = mu2, df = 3, phi = sd2, range = r2, nu = nu2, pts = pts)

  f = flatten(f)
  g = flatten(g)

  FAD = tryCatch(fadtest(f, g), error = function(x) NA)

  c(kolm(f, g)[2], quality(f, g)[2], FAD, bandtest(f, g)[2])
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
                   amp = amp,
                   feature = feature)

saveRDS(power, file = paste0("/home/trevorh2/matern/power/hetngp/",feature,"_sim_",i,"_het",".RDS"))

