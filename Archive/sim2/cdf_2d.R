rm(list = ls()); gc()
library(tictoc)

# get args
args = commandArgs(TRUE)
f1 = as.double(args[1])
f2 = as.integer(args[2])
d = as.integer(args[3])
l = as.integer(args[4])
i = as.integer(args[5])
seed = as.integer(args[6]) + i
sims = as.integer(args[7])

#### SIMULATION ####
gp2d = function(fields = 100, mu = 0, l = 30, pts = 30) {
  
  grid = 1:pts
  grid = expand.grid(grid, grid)
  distmat = as.matrix(dist(grid))
  
  # calc sigma with cov kernel
  sigma = exp(-distmat / l)
  
  sigma.eig = eigen(sigma)
  sigma.half = sigma.eig$vectors %*% diag(sqrt(sigma.eig$values)) %*% t(sigma.eig$vectors)
  
  gps = array(0, dim=c(pts, pts, fields))
  for(f in 1:fields) {
    gps[,,f] = (sigma.half %*% rnorm(pts^2)) + mu
  }
  return(gps)
}
flatten = function(mat) {
  matrix(mat, prod(dim(mat)[1:2]), dim(mat)[3])
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
ks_cdf = function(x, n = 100) {
  if(x < 0.3) return(0)
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
  
  ksf = rate*max(abs(ff.cdf - fg.cdf))
  ksg = rate*max(abs(gf.cdf - gg.cdf))
  
  1 - ks_cdf(max(ksf, ksg))
}

#### SIZE
set.seed(seed)

kmain = rep(0, sims)
for(s in 1:sims) {
  tic("Total")
  cat("Simulation ", s, "\n")

  gp1 = gp2d(f1, mu = 0, pts = d, l = l)
  gp2 = gp2d(f2, mu = 0, pts = d, l = l)

  gp1 = flatten(gp1)
  gp2 = flatten(gp2)

  kmain[s] = sk.test(gp1, gp2)

  cat("Size: ", mean(kmain[1:s] < 0.05), "\n")
  toc()
  cat("\n")
}
size = mean(kmain < 0.05)

out = c(seed, sims, f1, f2, d, l, size)
save(out, file = paste0("out/sim",i,".RData"))

