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
ks_pval = function(t, n = 20) {
  2*(sum(sapply(1:n, function(x) (-1)^(x-1) * exp(-2*(x^2)*t^2))))
}

joint.ks = function(f, g) {
  ff = meandepth(f, f)
  fg = meandepth(f, g)
  gf = meandepth(g, f)
  gg = meandepth(g, g)
  
  ffr = rank(ff) / length(ff)
  fgr = rank(fg) / length(fg)
  gfr = rank(gf) / length(gf)
  ggr = rank(gg) / length(gg)
  
  
  alpha = seq(1, 0, length.out = 1000)
  jcdf = sapply(alpha, function(t) mean(c(((ffr > t) & (fgr > t)), ((ggr > t) & (gfr > t)))))
  
  ks = max(abs(1-alpha - jcdf))
  rate = sqrt(sqrt(2)*length(jcdf))
  ks_pval(rate*ks)
}


f = gp1d(500, mu = 0, pts = 50, l = 10)
g = gp1d(500, mu = 0, pts = 50, l = 10)

joint.ks(f, g)



#### SIZE
set.seed(1)

sims = 2000
kdist = rep(0, sims)
for(s in 1:sims) {
  tic("Total")
  cat("Simulation ", s, "\n")
  
  gp1 = gp1d(500, mu = 0, pts = 50, l = 10)
  gp2 = gp1d(500, mu = 0, pts = 50, l = 10)
  
  kdist[s] = joint.ks(gp1, gp2) 
  
  toc()
  cat("\n")
}
mean(kdist < 0.05)
plot(kdist)

plt_funs(gp1, gp2)