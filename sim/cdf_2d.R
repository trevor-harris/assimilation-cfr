rm(list = ls()); gc()
library(extdepth)
library(tictoc)

# get args
args = commandArgs(TRUE)
f1 = as.double(args[1])
f2 = as.integer(args[2])
d = as.integer(args[3])
l = as.integer(args[4])
i = as.integer(args[5])
seed = as.integer(args[6])
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
meandepth = function(gmat, fmat) {
  apply(gmat, 2, function(x) mean(depth(x, fmat)))
}
ks.mean = function(f, g) {
  fed = meandepth(f, f)
  ged = meandepth(g, f)
  
  f.surv = rev(c(0, sort(fed)))
  g.cdf = sapply(1:length(f.surv), function(x) mean(ged > f.surv[x]))
  f.cdf = sapply(1:length(f.surv), function(x) mean(fed > f.surv[x]))
  
  ks = max(abs(f.cdf - g.cdf))
  rate = sqrt((ncol(g)*ncol(f)) / (ncol(g) + ncol(f)))
  ks_pval(rate*ks)
}
ks_pval = function(t, n = 5) {
  2*(sum(sapply(1:n, function(x) (-1)^(x-1) * exp(-2*(x^2)*t^2))))
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

  kmain[s] = ks.mean(gp1, gp2)

  cat("Size: ", mean(kmain[1:s] < 0.05), "\n")
  toc()
  cat("\n")
}
size = mean(kmain < 0.05)

out = c(seed, sims, f1, f2, d, l, size)
save(out, file = paste0("out/sim",i,".RData"))

