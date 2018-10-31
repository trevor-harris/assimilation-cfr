rm(list = ls()); gc()
library(tictoc)

# get args
args = commandArgs(TRUE)
i = as.integer(args[1])
seed = as.integer(args[2])
sims = as.integer(args[3])

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
ks_cdf = function(x, n = 10) {
  1 - 2*(sum(sapply(1:n, function(k) ((-1)^(k-1)) * exp(-2*(k^2)*(x^2)))))
}
perm.pval = function(est, table) {
  mean(est > table)
}

ksd = function(f, g) {
  ffxd = xdepth(f, f)
  gfxd = xdepth(g, f)
  fgxd = xdepth(f, g)
  ggxd = xdepth(g, g)
  
  tf = seq(0, 1, length.out = max(1000, 3*length(ffxd)))  
  ffr = sapply(tf, function(y) mean(ffxd <= y))
  gfr = sapply(tf, function(y) mean(gfxd <= y))
  
  tg = seq(0, 1, length.out = max(1000, 3*length(ggxd)))
  fgr = sapply(tg, function(y) mean(fgxd <= y))
  ggr = sapply(tg, function(y) mean(ggxd <= y))
  
  rate = sqrt((ncol(g)*ncol(f)) / (ncol(g) + ncol(f)))
  
  ksf = rate*max(abs(ffr - gfr))
  ksg = rate*max(abs(fgr - ggr))
  
  max(ksf, ksg)
}
ksd.perm = function(f, g, perms=500) {
  h = cbind(f, g)
  
  hn = ncol(h)
  fn = ncol(f)
  
  ksd.dist = rep(0, perms)
  
  for(i in 1:perms) {
    hstar = h[,sample(1:hn, hn)]
    ksd.dist[i] = ksd(hstar[,1:ncol(f)], hstar[,-(1:ncol(f))])
  }
  ksd.dist
}

#### SIZE
set.seed(seed + i)

nfun = 50
mu = 0
std = 1
pts = 15
rnge = 30
nperm = 500

empiri = rep(0, sims)
asympt = rep(0, sims)
for(s in 1:sims) {
  tic("Total")
  cat("Simulation ", s, "\n")

  gp1 = gp2d(nfun, mu = mu, pts = pts, l = rnge)
  gp2 = gp2d(nfun, mu = mu, pts = pts, l = rnge)

  gp1 = flatten(gp1)
  gp2 = flatten(gp2)
  
  ks = ksd(gp1, gp2)

  empiri[s] = perm.pval(ks, ksd.perm(gp1, gp2, nperm))
  asympt[s] = 1 - ks_cdf(ks)

  cat("Emprically: ", mean(empiri[1:s] < 0.05), "\n")
  cat("Asymptotic: ", mean(asympt[1:s] < 0.05), "\n")
  toc()
  cat("\n")
}
empiri = mean(empiri < 0.05)
asympt = mean(asympt < 0.05)

out = data.frame(seed = seed+i, i = i, nfun = nfun, mu = mu, std = std, rnge = rnge,
                 nperm = nperm, empiri = empiri, asympt = asympt)
save(out, file = paste0("out/size",i,".RData"))

