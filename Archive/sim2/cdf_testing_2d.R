rm(list = ls()); gc()
library(extdepth)
library(tictoc)


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
# ks.mean = function(f, g) {
#   fed = meandepth(f, f)
#   ged = meandepth(g, f)
#   
#   f.surv = rev(c(0, sort(fed)))
#   g.cdf = sapply(1:length(f.surv), function(x) mean(ged > f.surv[x]))
#   f.cdf = sapply(1:length(f.surv), function(x) mean(fed > f.surv[x]))
#   
#   ks = max(abs(f.cdf - g.cdf))
#   rate = sqrt((ncol(g)*ncol(f)) / (ncol(g) + ncol(f)))
#   ks_pval(rate*ks)
# }
ks.mean = function(f, g) {
  fed = meandepth(f, f)
  ged = meandepth(g, g)
  
  f.surv = rev(c(0, sort(fed)))
  gf.cdf = sapply(1:length(f.surv), function(x) mean(ged > f.surv[x]))
  
  g.surv = rev(c(0, sort(ged)))
  fg.cdf = sapply(1:length(f.surv), function(x) mean(fed > g.surv[x]))
  
  uni = seq(0, 1, length.out = length(gf.cdf))
  
  rate = sqrt((ncol(g)*ncol(f)) / (ncol(g) + ncol(f)))
  
  ksf = max(abs(uni - gf.cdf))
  ksg = max(abs(uni - fg.cdf))
  ks_pval(rate*max(ksf, ksg))
}
ks_pval = function(t, n = 5) {
  2*(sum(sapply(1:n, function(x) (-1)^(x-1) * exp(-2*(x^2)*t^2))))
}

#### SIZE
set.seed(1)

sims = 500
kmain = rep(0, sims)
for(s in 1:sims) {
  tic("Total")
  cat("Simulation ", s, "\n")
  
  gp1 = gp2d(500, mu = 0, pts = 30, l = 10)
  gp2 = gp2d(500, mu = 0, pts = 30, l = 10)
  
  gp1 = flatten(gp1)
  gp2 = flatten(gp2)
  
  kmain[s] = ks.mean(gp1, gp2)
  
  cat("Size: ", mean(kmain[1:s] < 0.05), "\n")
  toc()
  cat("\n")
}
mean(kmain < 0.05)
plot(kmain)

# plt_funs(gp1, gp2)
