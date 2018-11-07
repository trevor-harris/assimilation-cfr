# T TEST STUFF

rm(list = ls()); gc()
library(extdepth)
library(e1071)

# sims
gp1d = function(fields = 100, mu = 0, sd = 1, l = 50, pts = 50) {
  
  grid = 1:pts
  distmat = as.matrix(dist(grid))
  
  # calc sigma with cov kernel
  sigma = exp(-distmat / l)
  
  sigma.eig = eigen(sigma)
  sigma.half = sigma.eig$vectors %*% diag(sqrt(sigma.eig$values)) %*% t(sigma.eig$vectors)
  
  gps = matrix(0, pts, fields)
  for(f in 1:fields) {
    gps[,f] = (sigma.half %*% rnorm(pts, sd = sd)) + mu
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
xdepth = function(gmat, fmat) {
  apply(gmat, 2, function(x) mean(depth(x, fmat)))
}

# OG KS
ks_pval = function(t, n = 20) {
  2*(sum(sapply(1:n, function(x) (-1)^(x-1) * exp(-2*(x^2)*t^2))))
}
coverage = function(f, g) {
  ffxd = xdepth(f, f)
  gfxd = xdepth(g, f)
  
  ffr = sapply(sort(ffxd), function(y) mean(ffxd <= y))
  gfr = sapply(sort(ffxd), function(y) mean(gfxd <= y))
  
  rate = sqrt((ncol(g)*ncol(f)) / (ncol(g) + ncol(f)))
  ks = max(abs(ffr - gfr))
  (rate * ks)
}

set.seed(1023)
sims = 25000
boot = rep(0, sims)
trevor = rep(0, sims)
for(s in 1:sims) {
  tic("Total")
  cat("Simulation ", s, "\n")
  
  gp1 = gp1d(l = 10)
  gp2 = gp1d(l = 10)
  
  boot[s] = max(abs(rbridge(frequency = 100)))
  trevor[s] = coverage(gp1, gp2)
  
  toc()
  cat("\n")
}

dens = data.frame(freq = c(boot, trevor),
                  method = c(rep("boot", sims),
                             rep("trevor", sims)))

ggplot(dens, aes(x = freq, color = method, fill = method)) + 
  geom_density(alpha = 0.2)