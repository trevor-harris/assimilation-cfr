rm(list = ls())
gc()

library(ggplot2)
library(reshape2)

source("../research/assimilation-cfr/code/depth_tests.R")
source("../research/assimilation-cfr/code/depths.R")
source("../research/assimilation-cfr/code/simulation.R")

hall_cvm = function(f, g) {
  # f = gp1d()
  # g = gp1d()
  
  n = ncol(f)
  m = ncol(g)
  
  z = seq(min(f, g), max(f, g), length.out = n+m)
  f_hat = sapply(z, function(y) mean(apply(f, 2, function(x) sum(x > y) == 0)))
  g_hat = sapply(z, function(y) mean(apply(g, 2, function(x) sum(x > y) == 0)))
  
# plot(f_hat, type = "l")
# lines(g_hat, col = "red")
  
  rate = n*m / (n+m)^2
  tstat = rate*mean((f_hat - g_hat)^2)
  tstat_mu = 1/6 + 1/(6*(n + m))
  tstat_var = 1/45 * ((m + n -1) / (m + n)^2 ) * (4*m*n*(m+n) - 3*(m^2 + n^2) - 2*m*n)/(4*m*n)
  
  tnorm = (tstat - tstat_mu) / sqrt(tstat_var)
  1-pnorm(abs(tnorm))
  
  c(tstat, 1-pnorm(abs(tnorm)))
}

hall_ks = function(f, g) {
  n = ncol(f)
  m = ncol(g)
  
  z = seq(min(f, g), max(f, g), length.out = n+m)
  f_hat = sapply(z, function(y) mean(apply(f, 2, function(x) sum(x > y) == 0)))
  g_hat = sapply(z, function(y) mean(apply(g, 2, function(x) sum(x > y) == 0)))
  
  plot(f_hat, type = "l")
  lines(g_hat, col = "red")
  
  # f_hat = sapply(1:ncol(z), function(t1) {
  #   mean(sapply(1:ncol(f), function(t2) sum(f[,t2] <= z[,t1]) == nrow(f)))
  # })
  # g_hat = sapply(1:ncol(z), function(t1) {
  #   mean(sapply(1:ncol(g), function(t2) sum(g[,t2] <= z[,t1]) == nrow(g)))
  # })
  
  tstat = max(abs(f_hat - g_hat))
  rate = sqrt(n*m / (n+m))
  
  c(tstat, 1-ks_cdf(rate*tstat))
}

kvals = rep(0, 500)
pvals = rep(0, 500)
for (i in 1:500){
  f = gp1d(300)
  g = gp1d(301)
  
  k = hall_cvm(f, g)
  kvals[i] = k[1]
  pvals[i] = k[2]
}
hist(kvals, breaks = 40)

plot(rate*(kvals))

z = seq(0, 2, length.out = 500)
plot(sapply(z, function(x) mean(rate*kvals < x)))
plot(sapply(z, function(x) mean(rate*kvals < x)))

plot(sapply(z, function(x) ks_cdf(x)^0.5))
plot(sapply(z, function(x) ks_cdf(x)^2))


t = seq(0, 1, length.out = 50)
mu = 2*sin(t * pi)

plot(mu)

f = gp1d()
g = gp1d()
hall_cvm(f, g)

plt_funs(f, g)


z = cbind(f, g)

sapply(1:200, function(x) sum(f[,1] < z[,x]) == 50)

f_hat = sapply(1:ncol(z), function(t1) {
  mean(sapply(1:ncol(f), function(t2) sum(f[,t2] <= z[,t1]) == nrow(f)))
})

plot(f[,1], type = "l")
lines(z[,5], col = "red")
sum(f[,1] <= z[,5])

g_hat = sapply(1:ncol(z), function(t1) {
  mean(sapply(1:ncol(g), function(t2) sum(g[,t2] <= z[,t1]) == nrow(g)))
})

plot(sort(f_hat), type = "l")
lines(sort(g_hat), col = "red")

f_hat = sapply(z, function(y) mean(apply(f, 2, function(x) sum(x > y) == 0)))
g_hat = sapply(z, function(y) mean(apply(g, 2, function(x) sum(x > y) == 0)))


cdf = cbind(f_hat, g_hat)
cdf = melt(cdf)

ggplot() + 
  geom_path(data = cdf, aes(x = Var1, y = value, color = Var2),
            size = 1) +
  theme_classic() +
  ylab("") +
  xlab("Function") +
  scale_color_manual(labels = c("F", "G"), values = c("#F8766D", "#00BFC4")) +
  guides(color=guide_legend("Distribution"))

rate = ncol(f)*ncol(g) / (ncol(f) + ncol(g))^2
rate*mean((f_hat - g_hat)^2)

n = ncol(f)
m = ncol(g)

tstat = mean((f_hat - g_hat)^2)
tstat_mu = 1/6 + 1/(6*(n + m))
tstat_var = 1/45 * ((m + n -1) / (m + n)^2 ) * (4*m*n*(m+n) - 3*(m^2 + n^2) - 2*m*n)/(4*m*n)

tnorm = (tstat - tstat_mu) / sqrt(tstat_var)
pnorm(abs(tnorm))
