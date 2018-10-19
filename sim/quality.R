# T TEST STUFF

rm(list = ls()); gc()
library(extdepth)
library(tictoc)

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


f = gp1d(10)
g = gp1d(10)

fxd = xdepth(f, f)
gxd = xdepth(g, f)

r = sapply(gxd, function(y) mean(y >= fxd))
1-pnorm((mean(r) - 0.5) / sqrt((1/ncol(f) + 1/ncol(g))/12))

plot(r)
plot(sort(r))

alpha = seq(0, 1, length.out = ncol(f))
t = sapply(alpha, function(a) mean(gxd <= a))
tr = sapply(alpha, function(a) mean(fxd <= a))
ks.test(t, tr)$p.value

t.test(fxd, gxd)$p.value

plot(alpha, type = "l")
lines(sort(r), col = "red")

plot(tr, type = "l")
lines(t, col = "blue")

quality = function(f, g) {
  fxd = xdepth(f, f)
  gxd = xdepth(g, f)
  
  r = sapply(gxd, function(y) mean(y >= fxd))
  1-pnorm((mean(r) - 0.5) / sqrt((1/ncol(f) + 1/ncol(g))/12))
}

# T
coverage = function(f, g) {
  fxd = xdepth(f, f)
  gxd = xdepth(g, f)
  
  t.test(fxd, gxd)$p.value
}

# OG KS
ks_pval = function(t, n = 20) {
  2*(sum(sapply(1:n, function(x) (-1)^(x-1) * exp(-2*(x^2)*t^2))))
}
coverage = function(f, g) {
  fxd = xdepth(f, f)
  gxd = xdepth(g, f)
  
  ffr = sort(sapply(fxd, function(y) mean(y >= fxd)))
  gfr = sort(sapply(gxd, function(y) mean(y >= fxd)))
  
  
  rate = sqrt((ncol(g)*ncol(f)) / (ncol(g) + ncol(f)))
  ks = max(abs(ffr - gfr))
  ks_pval(rate * ks)
}

# SYM KS
ks_cdf = function(x, n = 10) {
  1 - 2*(sum(sapply(1:n, function(k) ((-1)^(k-1)) * exp(-2*(k^2)*(x^2)))))
}

coverage = function(f, g) {
  ffxd = xdepth(f, f)
  gfxd = xdepth(g, f)
  
  ffr = sapply(sort(ffxd), function(y) mean(ffxd <= y))
  gfr = sapply(sort(ffxd), function(y) mean(gfxd <= y))
  
  fgxd = xdepth(f, g)
  ggxd = xdepth(g, g)
  
  fgr = sapply(sort(ggxd), function(y) mean(fgxd <= y))
  ggr = sapply(sort(ggxd), function(y) mean(ggxd <= y))
  
  rate = sqrt((ncol(g)*ncol(f)) / (ncol(g) + ncol(f)))
  
  ksf = rate*max(abs(ffr - gfr))
  ksg = rate*max(abs(fgr - ggr))
  
  1 - ks_cdf(max(ksf, ksg))^2
}


coverage2 = function(f, g) {
  ffxd = xdepth(f, f)
  gfxd = xdepth(g, f)
  fgxd = xdepth(f, g)
  ggxd = xdepth(g, g)

  tf = seq(0, 1, length.out = max(1000, 3*length(ffxd)))  
  ffr = sapply(tf, function(y) mean(ffxd <= y))
  gfr = sapply(tf, function(y) mean(gfxd <= y))
  
  tg = seq(0, 1, length.out = max(1000, 3*length(ggxd)))
  fgr = sapply(tf, function(y) mean(fgxd <= y))
  ggr = sapply(tf, function(y) mean(ggxd <= y))
  
  rate = sqrt((ncol(g)*ncol(f)) / (ncol(g) + ncol(f)))
  
  ksf = rate*max(abs(ffr - gfr))
  ksg = rate*max(abs(fgr - ggr))
  
  1 - ks_cdf(max(ksf, ksg))^2
}

set.seed(0723)
sims = 1000
regina = rep(0, sims)
trevor = rep(0, sims)
for(s in 1:sims) {
  tic("Total")
  cat("Simulation ", s, "\n")
  
  gp1 = gp1d(50, l = 30)
  gp2 = gp1d(50, l = 30)
  
  regina[s] = coverage(gp1, gp2)
  trevor[s] = coverage2(gp1, gp2)
  
  toc()
  cat("Regina: ", mean(regina[1:s] < 0.05), "\n")
  cat("Trevor: ", mean(trevor[1:s] < 0.05), "\n")
  cat("\n")
}

