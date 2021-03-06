rm(list = ls()); gc()
library(extdepth)
library(tictoc)

# sims
gp1d = function(fields = 100, mu = 0, l = 50, pts = 50) {
  
  grid = 1:pts
  distmat = as.matrix(dist(grid))
  
  # calc sigma with cov kernel
  sigma = exp(-distmat / l)
  
  sigma.eig = eigen(sigma)
  sigma.half = sigma.eig$vectors %*% diag(sqrt(sigma.eig$values)) %*% t(sigma.eig$vectors)
  
  gps = matrix(0, pts, fields)
  for(f in 1:fields) {
    gps[,f] = (sigma.half %*% rnorm(pts)) + mu
  }
  return(gps)
}
plt_funs = function(f, g) {
  f = as.matrix(f)
  g = as.matrix(g)
  plot(f[,1], type = "l", ylim = c(min(cbind(f, g)), max(cbind(f, g))), col = "red")
  for(i in 2:ncol(f)) {
    lines(f[,i], col = "red")
  }
  for(i in 1:ncol(g)) {
    lines(g[,i], col = "blue")
  }
}

# depth CDF
depth = function(g, fmat) {
  
  # Computes the depth values of a function with respect to a set of functions (fmat)
  fn = ncol(fmat)
  depth = rep(0, length(g))
  
  for (row in 1:nrow(fmat)) {
    diff = abs(sum(sign(g[row] - fmat[row,])))
    depth[row] = 1 - (diff / fn)
  }
  
  return(depth)
}
meandepth = function(gmat, fmat) {
  apply(gmat, 2, function(x) mean(depth(x, fmat)))
}
ks.mean = function(f, g) {
  fed = meandepth(f, f)
  ged = meandepth(g, g)
  
  f.surv = rev(c(0, sort(fed)))
  g.cdf = sapply(1:length(f.surv), function(x) mean(ged > f.surv[x]))
  f.cdf = sapply(1:length(f.surv), function(x) mean(fed > f.surv[x]))
  
  # plot(g.cdf, type = "l", col = "blue")
  # lines(f.cdf, col = "red")
  
  ks = max(abs(f.cdf - g.cdf))
  rate = sqrt((ncol(g)*ncol(f)) / (ncol(g) + ncol(f)))
  ks_pval(rate*ks)
}
ks_pval = function(t, n = 20) {
  2*(sum(sapply(1:n, function(x) (-1)^(x-1) * exp(-2*(x^2)*t^2))))
}


ks1.mean = function(f, g) {
  fed = meandepth(f, f)
  ged = meandepth(g, f)
  
  f.surv = rev(c(0, sort(fed)))
  g.cdf = sapply(1:length(f.surv), function(x) mean(ged >= f.surv[x]))
  f.cdf = seq(0, 1, length.out = length(g.cdf))
  
  # plot(g.cdf, type = "l", col = "blue")
  # lines(f.cdf, col = "red")
  
  ks = max(abs(f.cdf - g.cdf))
  rate = sqrt(ncol(g))
  ks_pval(rate*ks)
}


#### SIZE
set.seed(1)

sims = 1000
kdist = rep(0, sims)
for(s in 1:sims) {
  tic("Total")
  cat("Simulation ", s, "\n")
  
  gp1 = gp1d(500, mu = 0, pts = 50, l = 10)
  gp2 = gp1d(500, mu = 0, pts = 50, l = 10)
  
  kdist[s] = ks.mean(gp1, gp2) 
  
  toc()
  cat("\n")
}
mean(kdist < 0.05)
plot(kdist)

plt_funs(gp1, gp2)



#### SIN AND COSINE MEANS
mu1 = cos(seq(-2*pi, 2*pi, length.out = 50))
mu2 = sin(seq(-pi, pi, length.out = 50))
plot(mu1, type = "l")
lines(mu2)

set.seed(1)

sims = 1000
kdist = rep(0, sims)
for(s in 1:sims) {
  tic("Total")
  cat("Simulation ", s, "\n")
  
  gp1 = gp1d(500, mu = mu1, pts = 50, l = 10)
  gp2 = gp1d(500, mu = mu2, pts = 50, l = 10)
  
  kdist[s] = ks.mean(gp1, gp2) 
  
  toc()
  cat("\n")
}
wave = mean(kdist < 0.05)
plot(kdist)

plt_funs(gp1, gp2)



#### TREND MEANS
mu1 = seq(0, 1, length.out = 50)
mu2 = seq(0, 0.25, length.out = 50)

plot(mu1, type = "l")
lines(mu2)

set.seed(1)

sims = 1000
kdist = rep(0, sims)
for(s in 1:sims) {
  tic("Total")
  cat("Simulation ", s, "\n")
  
  gp1 = gp1d(500, mu = mu1, pts = 50, l = 10)
  gp2 = gp1d(500, mu = mu2, pts = 50, l = 10)
  
  kdist[s] = ks.mean(gp1, gp2) 
  
  toc()
  cat("\n")
}
trend = mean(kdist < 0.05)
plot(kdist)

plt_funs(gp1, gp2)



#### POWER MEANS
mu1 = (seq(0, 1, length.out = 50))
mu2 = (seq(0, 1, length.out = 50))^4

plot(mu1, type = "l")
lines(mu2)

set.seed(1)

sims = 1000
kdist = rep(0, sims)
for(s in 1:sims) {
  tic("Total")
  cat("Simulation ", s, "\n")
  
  gp1 = gp1d(500, mu = mu1, pts = 50, l = 10)
  gp2 = gp1d(500, mu = mu2, pts = 50, l = 10)
  
  kdist[s] = ks.mean(gp1, gp2) 
  
  toc()
  cat("\n")
}
mean(kdist < 0.05)
plot(kdist)

plt_funs(gp1, gp2)




#### GP MEANS
pts = 1000
set.seed(5)
mu1 = rep(0, pts)
mu2 = gp1d(fields = 1, pts = pts, l = 500)
mu2 = rep(0, pts)

plot(mu1, type = "l")
lines(mu2)

set.seed(1)

sims = 100
kdist = rep(0, sims)
for(s in 1:sims) {
  tic("Total")
  cat("Simulation ", s, "\n")
  
  gp1 = gp1d(100, mu = mu1, pts = pts, l = 10)
  gp2 = gp1d(100, mu = mu2, pts = pts, l = 10)
  
  kdist[s] = ks.mean(gp1, gp2) 
  
  toc()
  cat("\n")
}
mean(kdist < 0.05)
plot(kdist)

plt_funs(gp1, gp2)