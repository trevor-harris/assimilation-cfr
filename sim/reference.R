# T TEST STUFF

# rm(list = ls()); gc()
# library(extdepth)
# library(tictoc)

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
gp2d = function(fields = 100, mu = 0, sd = 1, l = 30, pts = 30) {
  grid = 1:pts
  grid = expand.grid(grid, grid)
  distmat = as.matrix(dist(grid))
  
  # calc sigma with cov kernel
  sigma = exp(-distmat / l)
  
  sigma.eig = eigen(sigma)
  sigma.half = sigma.eig$vectors %*% diag(sqrt(sigma.eig$values)) %*% t(sigma.eig$vectors)
  
  gps = array(0, dim=c(pts, pts, fields))
  for(f in 1:fields) {
    gps[,,f] = (sigma.half %*% rnorm(pts^2, sd = sd)) + mu
  }
  return(gps)
}
flatten = function(mat) {
  matrix(mat, prod(dim(mat)[1:2]), dim(mat)[3])
}
plt_funs = function(f, g, main = "Functions") {
  
  if(missing(g)) {
    f = as.matrix(f)
    plot(f[,1], type = "l", ylim = c(min(f), max(f)), col = "red", main = main)
    for(i in 2:ncol(f)) {
      lines(f[,i], col = "red")
    }
  }
  else {
    f = as.matrix(f)
    g = as.matrix(g)
    plot(f[,1], type = "l", ylim = c(min(cbind(f, g)), max(cbind(f, g))), col = "red", main = main)
    for(i in 2:ncol(f)) {
      lines(f[,i], col = "red")
    }
    for(i in 1:ncol(g)) {
      lines(g[,i], col = "blue")
    }
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

# SYM KS
ks_cdf = function(x, n = 100) {
  if(x < 0.05) return(0)
  1 - 2*(sum(sapply(1:n, function(k) ((-1)^(k-1)) * exp(-2*(k^2)*(x^2)))))
}

sk.test = function(f, g) {
  ff.xd = xdepth(f, f)
  fg.xd = xdepth(f, g)
  
  gg.xd = xdepth(g, g)
  gf.xd = xdepth(g, f)
  
  tf = seq(0, 1, length.out = max(1000, 3*length(ff.xd)))  
  ff.cdf = sapply(tf, function(y) mean(ff.xd <= y))
  gf.cdf = sapply(tf, function(y) mean(gf.xd <= y))
  
  tg = seq(0, 1, length.out = max(1000, 3*length(gg.xd)))
  fg.cdf = sapply(tg, function(y) mean(fg.xd <= y))
  gg.cdf = sapply(tg, function(y) mean(gg.xd <= y))
  
  rate = sqrt((ncol(g)*ncol(f)) / (ncol(g) + ncol(f)))
  
  ks1 = rate*max(abs(ff.cdf - fg.cdf))
  ks2 = rate*max(abs(gf.cdf - gg.cdf))
  ks3 = rate*max(abs(ff.cdf - gg.cdf))
  ks4 = rate*max(abs(fg.cdf - gf.cdf))
  
  1 - ks_cdf(max(ks1, ks2, ks3, ks4))
}

quality = function(f, g) {
  fxd = xdepth(f, f)
  gxd = xdepth(g, f)
  
  r = sapply(gxd, function(y) mean(y >= fxd))
  1 - pnorm((mean(r) - 0.5) / sqrt((1/ncol(f) + 1/ncol(g))/12))
}

