rm(list = ls()); gc()

library(fdasrvf)

plt_funs = function(fmat) {
  plot(fmat[,1], type = "l", ylim = c(min(fmat), max(fmat)))
  for(i in 2:ncol(fmat)) {
    lines(fmat[,i])
  }
}
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

crossings = function(f1, f2) {
  f1 = as.vector(f1)
  f2 = as.vector(f2)
  
  h1 = tail(f1, -1)
  h2 = tail(f2, -1)
  
  sum(((f1 > f2)[-length(f1)] & (h1 <= h2)) | 
        ((f1 < f2)[-length(f1)] & (h1 >= h2)))
}
cross = function(fmat) {
  nf = ncol(fmat)
  crossmat = matrix(0, nf, nf)
  for(i in 1:nf) {
    for(j in i:nf) {
      crossmat[i, j] = crossmat[j, i] = crossings(fmat[,i], fmat[,j])
    }
  }
  sum(crossmat)
}


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
vardepth = function(gmat, fmat) {
  apply(gmat, 2, function(x) var(depth(x, fmat)))
}

# cross(simu_data$f)
# exp(ncol(simu_data$f))

gp01 = flatten(gp2d(500, l = 1))
gp05 = flatten(gp2d(500, l = 5))
gp10 = flatten(gp2d(500, l = 10))
gp100 = flatten(gp2d(500, l = 100))
gp1000 = flatten(gp2d(500, l = 1000))

# c01 = cross(gp01)
# c05 = cross(gp05)
# c10 = cross(gp10)
# c100 = cross(gp100)
# c1000 = cross(gp1000)

v01 = meandepth(gp01, gp01)
v05 = meandepth(gp05, gp05)
v10 = meandepth(gp10, gp10)
v100 = meandepth(gp100, gp100)
v1000 = meandepth(gp1000, gp1000)

v = c(v01, v05, v10, v100, v1000)
ylims = c(min(v), max(v))

plot(v01, ylim = ylims)
points(v05, col = "blue")
points(v10, col = "red")
points(v100, col = "green")
points(v1000, col = "purple")

sd(v01)
sd(v05)
sd(v10)
sd(v100)
sd(v1000)

# sd of a uniform(0, 1)
sqrt(1/12)

# closer your data is to the theoretical the better test will work?



rescale = function(x) {
  z = x
  
  # z = (x - mean(x)) / sd(x)
  # plot(z)
  
  d = (z - min(z)) / (max(z) - min(z))
  # plot(d)
  return(d)
  
}
n01 = rescale(v01)
n05 = rescale(v05)
n10 = rescale(v10)

plot(rescale(v01), ylim = c(0, 1))
points(rescale(v05), col = "blue")
points(rescale(v10), col = "red")

hist(rescale(v10))


# does rescaling help?
ks_pval = function(t, n = 5) {
  2*(sum(sapply(1:n, function(x) (-1)^(x-1) * exp(-2*(x^2)*t^2))))
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
ks.mean2 = function(f, g) {
  fed = rescale(meandepth(f, f) / sqrt(vardepth(f, f)))
  ged = rescale(meandepth(g, f) / sqrt(vardepth(g, f)))
  
  f.surv = rev(c(0, sort(fed)))
  g.cdf = sapply(1:length(f.surv), function(x) mean(ged > f.surv[x]))
  f.cdf = sapply(1:length(f.surv), function(x) mean(fed > f.surv[x]))
  
  ks = max(abs(f.cdf - g.cdf))
  rate = sqrt((ncol(g)*ncol(f)) / (ncol(g) + ncol(f)))
  ks_pval(rate*ks)
}

k1 = rep(0, 100)
k2 = rep(0, 100)
for(s in 1:100) {
  tic("Total")
  cat("Simulation ", s, "\n")

  gp1 = flatten(gp2d(l = 1, pts = 10))
  gp2 = flatten(gp2d(l = 1, pts = 10))
  
  k1[s] = ks.mean(gp1, gp2)
  k2[s] = ks.mean2(gp1, gp2)
  
  cat("Size1: ", mean(k1[1:s] < 0.05), "\n")
  cat("Size2: ", mean(k2[1:s] < 0.05), "\n")
  toc()
  cat("\n")
}

plot(k1)
points(k2, col = "red")



