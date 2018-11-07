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

# 1way t test version
tcover1 = function(f, g) {
  ff.xd = xdepth(f, f)
  gf.xd = xdepth(g, f)
  
  t.test(ff.xd, gf.xd)$p.value
}

tcover2 = function(f, g) {
  ff.xd = xdepth(f, f)
  gf.xd = xdepth(g, f)
  
  gg.xd = xdepth(g, g)
  fg.xd = xdepth(f, g)
  
  stat = max(t.test(ff.xd, gf.xd)$statistic, t.test(gg.xd, fg.xd)$statistic)
  1 - pnorm(stat)^2
}

nscore = function(xd1, xd2) {
  xd = c(xd1, xd2)
  nxd = qnorm(rank(xd) / (length(xd) + 1))
  return(list(xd1 = nxd[1:length(xd1)],
              xd2 = nxd[(length(xd1) + 1):length(xd)]))
}

tscover = function(f, g) {
  ff.xd = xdepth(f, f)
  gf.xd = xdepth(g, f)
  
  gg.xd = xdepth(g, g)
  fg.xd = xdepth(f, g)
  
  fscores = nscore(ff.xd, gf.xd)
  gscores = nscore(gg.xd, fg.xd)
  
  stat = max(t.test(fscores$xd1, fscores$xd2)$statistic,
             t.test(gscores$xd1, gscores$xd2)$statistic)
  1 - pnorm(stat)^2
}


gp1 = gp1d()
gp2 = gp1d()
gp2 = cbind(gp1d(fields = 50, mu = 2, l = 10),
            gp1d(fields = 50, mu = -2, l = 10))

plt_funs(gp1, gp2)

set.seed(1022)
sims = 500
cov1 = rep(0, sims)
cov2 = rep(0, sims)
cov3 = rep(0, sims)
for(s in 1:sims) {
  tic("Total")
  cat("Simulation ", s, "\n")
  
  gp1 = gp1d(fields = 500, l = 10)
  gp2 = gp1d(fields = 500, l = 10)
  
  # modify to test either location, scale, amplitude, phase
  # gp2 = cbind(gp1d(fields = 50, mu = 2, l = 10),
  #             gp1d(fields = 50, mu = -2, l = 10))
  
  cov1[s] = tcover1(gp1, gp2) 
  cov2[s] = tcover2(gp1, gp2)
  cov3[s] = tscover(gp1, gp2) 
  
  toc()
  cat(mean(cov1[1:s] < 0.05), "\n")
  cat(mean(cov2[1:s] < 0.05), "\n")
  cat(mean(cov3[1:s] < 0.05), "\n")
  cat("\n")
}
mean(cov1 < 0.05)
mean(cov2 < 0.05)
mean(cov3 < 0.05)

plot(cov1)
plot(cov2)
plot(cov3)

