rm(list = ls()); gc()
library(tictoc)
library(future.apply)

# get args
n = 400
pts = 50
mu1 = 0
mu2 = seq(-0.5, 0.5, length.out = pts)
sd1 = 1
sd2 = 1
r1 = 20
r2 = 20
i = 1
seed = 102391 + i
feature = "mean"
sims = 250

# source("/home/trevorh2/assimilation-cfr/code/depth_tests.R")
# source("/home/trevorh2/assimilation-cfr/code/depths.R")
# source("/home/trevorh2/assimilation-cfr/code/simulation.R")

source("../../../code/depth_tests.R")
source("../../../code/depths.R")
source("../../../code/simulation.R")

#### POWER MEAN
plan(multiprocess)
set.seed(seed)


nmu = 20
out3 = matrix(0, nmu, 3)

for(m in 1:nmu) {
  mu2 = gp1d(1, mu = 0, sd = m/nmu, l = r1, pts = pts)
  mu2 = mu2 - mean(mu2)
  
  tic("Total")
  cat("Simulation ", m, "\n")
  
  pvals = future_sapply(1:sims, function(x) {
    
    f = gp1d(fields = n, mu = mu1, sd = sd1, l = r1, pts = pts)
    g = gp1d(fields = n+1, mu = mu2, sd = sd2, l = r2, pts = pts)
    
    c(kolm(f, g)[2], quality(f, g)[2])
  })
  
  out3[m,1] = mean(pvals[1,] < 0.05)
  out3[m,2] = mean(pvals[2,] < 0.05)
  out3[m,3] = mean(mu2^2)
  
  toc()
  
}


out4 = as.data.frame(rbind(out3[,c(1, 3)], out3[,c(2, 3)]))
out4["Stat"] = c(rep("K", nmu), rep("Q", nmu))

ggplot(out4, aes(x = V2, y = V1, color = Stat)) +
  geom_line() +
  theme_classic() +
  ggtitle("Stochastic mean change")



### POWER SD
sims = 250

nmu = 40
out3 = matrix(0, nmu, 4)

for(m in 1:nmu) {
  sd2 = gp1d(1, mu = 1, sd = nmu/nmu, l = r1, pts = pts)
  sd2 = abs(sd2 - mean(sd2) + 1)
  
  plot(sd2, type = "l")
  
  tic("Total")
  cat("Simulation ", m, "\n")
  
  pvals = future_sapply(1:sims, function(x) {
    
    f = gp1d(fields = n, mu = mu1, sd = sd1, l = r1, pts = pts)
    g = gp1d(fields = n+1, mu = mu1, sd = sd2, l = r2, pts = pts)
    
    c(kolm(f, g)[2], quality(f, g)[2])
  })
  
  out3[m,1] = mean(pvals[1,] < 0.05)
  out3[m,2] = mean(pvals[2,] < 0.05)
  out3[m,3] = mean((sd2 - 1)^2)
  out3[m,4] = max(abs(sd2 - 1))
  
  toc()
  
}

out4 = as.data.frame(rbind(out3[,c(1, 3)], out3[,c(2, 3)]))
out4["Stat"] = c(rep("K", nmu), rep("Q", nmu))

ggplot(out4, aes(x = V2, y = V1, color = Stat)) +
  geom_line() +
  theme_classic()  +
  ggtitle("Stochastic sd change")




### POWER SD
sims = 250

nmu = 40
out3 = matrix(0, nmu, 4)

for(m in 1:nmu) {
  cr2 = gp1d(1, mu = 0, sd = 20*m/nmu, l = r1, pts = pts)
  cr2 = abs(cr2 - mean(cr2) + 30)
  
  # plot(cr2, type = "l")
  
  tic("Total")
  cat("Simulation ", m, "\n")
  
  pvals = future_sapply(1:sims, function(x) {
    
    f = gp1d(fields = n, mu = mu1, sd = sd1, l = r1, pts = pts)
    g = gp1d(fields = n+1, mu = mu1, sd = sd1, l = cr2, pts = pts)
    
    c(kolm(f, g)[2], quality(f, g)[2])
  })
  
  out3[m,1] = mean(pvals[1,] < 0.05)
  out3[m,2] = mean(pvals[2,] < 0.05)
  out3[m,3] = mean((sd2 - 1)^2)
  out3[m,4] = max(abs(sd2 - 1))
  
  toc()
  
}

out4 = as.data.frame(rbind(out3[,c(1, 3)], out3[,c(2, 3)]))
out4["Stat"] = c(rep("K", nmu), rep("Q", nmu))

ggplot(out4, aes(x = V2, y = V1, color = Stat)) +
  geom_line() +
  theme_classic()  +
  ggtitle("Stochastic sd change")
 
