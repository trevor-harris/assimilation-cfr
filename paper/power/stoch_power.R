rm(list = ls()); gc()
library(plotly)
library(fdasrvf)
library(tictoc)
library(future.apply)

# get args
n = 200
pts = 20
mu1 = 0
mu2 = 0
sd1 = 1
sd2 = 1
r1 = 10
r2 = 10
seed = 102391
sims = 100

# source("/home/trevorh2/assimilation-cfr/code/depth_tests.R")
# source("/home/trevorh2/assimilation-cfr/code/depths.R")
# source("/home/trevorh2/assimilation-cfr/code/simulation.R")

source("../../../code/depth_tests.R")
source("../../../code/depths.R")
source("../../../code/simulation.R")

plan(multiprocess)
set.seed(seed)

params = 20
out = matrix(0, params, 3)


#### MEAN changes
for(m in 1:params) {
  mu2 = flatten(gp2d(1, mu = 0, sd = m/params, l = r1, pts = pts))
  mu2 = mu2 - mean(mu2)
  
  tic(paste0("Param ", m))
  
  pvals = future_sapply(1:sims, function(x) {
    
    f = gp2d(fields = n, mu = mu1, sd = sd1, l = r1, pts = pts)
    g = gp2d(fields = n+1, mu = mu2, sd = sd1, l = r1, pts = pts)
    
    f = flatten(f)
    g = flatten(g) 
    
    c(kolm(f, g)[2], quality(f, g)[2])
  })
  
  out[m,1] = mean(pvals[1,] < 0.05)
  out[m,2] = mean(pvals[2,] < 0.05)
  out[m,3] = mean(mu2^2)
  
  toc()
}

out2 = as.data.frame(rbind(out[,c(1, 3)], out[,c(2, 3)]))
out2["Stat"] = c(rep("K", params), rep("QI", params))
out2["param"] = "mean"

ggplot(out2, aes(x = V2, y = V1, color = Stat)) +
  geom_line() +
  theme_classic() +
  ggtitle("Stochastic mean change")


#### SD changes
sigmoid = function(r) {
  x = seq(-r, r, length.out = 20)
  1 / (1 + exp(-x))
}
relu = function(x) {
  sapply(x, function(t) max(t, 0.01))
}
for(m in 1:params) {

  sd2 = gp2d(1, mu = 0, sd = 5*m/params, l = 100, pts = 20)
  sd2 = sd2 - mean(sd2) + 0.5
  sd2 = matrix(relu(sd2), pts, pts)
  
  plot_ly() %>% add_surface(z = ~sd2)
  sd2 = smooth.data(t(smooth.data(sd2, 3)), 3)
  sd2 = as.vector(sd2)
  
  tic(paste0("Param ", m))
  
  pvals = future_sapply(1:sims, function(x) {
    
    f = gp2d(fields = n, mu = mu1, sd = 1, l = r1, pts = pts)
    g = gp2d(fields = n, mu = mu1, sd = sd2, l = r1, pts = pts)
    
    f = flatten(f)
    g = flatten(g) 
    
    c(kolm(f, g)[2], quality(f, g)[2])
  })
  
  out[m,1] = mean(pvals[1,] < 0.05)
  out[m,2] = mean(pvals[2,] < 0.05)
  out[m,3] = mean((sd2 - 1)^2)
  
  toc()
}

out3 = as.data.frame(rbind(out[,c(1, 3)], out[,c(2, 3)]))
out3["Stat"] = c(rep("K", params), rep("Q", params))
out3["param"] = "Standard Deviation"

ggplot(out3, aes(x = V2, y = V1, color = Stat)) +
  geom_line() +
  theme_classic() +
  ggtitle("Stochastic sd change")


power2d = rbind(out2, out3)
ggplot(power2d, aes(x=V2, y=V1, color=Stat)) +
  geom_line() +
  geom_hline(yintercept = 0.05) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5)) +
  xlab("Parameter Change") +
  ylab("Power") +
  facet_wrap(. ~ param, nrow = 1, scales = "free_x")

ggsave("../stoch_multi2d.png", width = 9, heigh = 3)

saveRDS(power2d, "")



m = 20
sd2 = gp2d(1, mu = 0, sd = 1 + m/params, l = 200, pts = 20)
sd2 = sd2 - mean(sd2) + 1
sd2 = matrix(relu(sd2), pts, pts)
mean((sd2 - 1)^2)
plot_ly() %>% add_surface(z = ~sd2)

sd2 = smooth.data(t(smooth.data(sd2, 2)), 2)

plot_ly() %>% add_surface(z = ~sd2)
