rm(list = ls()); gc()
library(tictoc)

# get args
args = commandArgs(TRUE)
i = as.integer(args[1])
seed = as.integer(args[2])
sims = as.integer(args[3])

#### SIMULATION ####
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

#### COMPARE CDFS


nfun = 2000
mu = 0
std = 1
pts = 15
rnge = 30
nperm = 200

# this seed ensures the functions are the same across runs
set.seed(seed)
gp1 = gp1d(nfun, mu = mu, pts = pts, l = rnge)
gp2 = gp1d(nfun, mu = mu, pts = pts, l = rnge)

# this seed increment ensures the permutations are different across runs
set.seed(seed + i)
perms = ksd.perm(gp1, gp2, nperm)
out = data.frame(i = i, nfun = nfun, mu = mu, std = std, pts = pts, rnge = rnge, nperm = nperm, 
                 perms = perms)

save(out, file = paste0("out/size",i,".RData"))

