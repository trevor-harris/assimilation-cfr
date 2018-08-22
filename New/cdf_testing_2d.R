rm(list = ls()); gc()
library(extdepth)
library(tictoc)


#### SIMULATION ####
create_basis = function(field, basis.pts = 10, r = 1.5) {
  
  x.dim = dim(field)[1]
  y.dim = dim(field)[2]
  
  x.loc = as.integer(seq(1, x.dim, length.out = basis.pts))
  y.loc = as.integer(seq(1, y.dim, length.out = basis.pts))
  
  # find the radius (1.5 times the minimal distance)
  # will either be the distance b/w the first point and the next (horizontally or vertically)
  rad = r * min(abs(x.loc[1] - x.loc[2]), abs(y.loc[1] - y.loc[2]))
  
  n.loc = basis.pts * basis.pts
  grid = expand.grid(x.loc, y.loc)
  basis = matrix(0, x.dim*y.dim, n.loc)
  
  for (i in 1:n.loc) {
    
    # dist = (lat - center.lat)^2 + (lon - center.lon)^2
    dist = c(outer((1:x.dim - grid[i,][[1]])^2, (1:y.dim - grid[i,][[2]])^2, "+"))
    
    # use the bisquare (radial) basis
    basis[,i] = (1 - dist/rad^2)^2 * (sqrt(dist) < rad)
  }
  
  return(basis)
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
ks.mean = function(f, g) {
  fed = meandepth(f, f)
  ged = meandepth(g, f)
  
  f.surv = rev(c(0, sort(fed)))
  g.cdf = sapply(1:length(f.surv), function(x) mean(ged > f.surv[x]))
  f.cdf = sapply(1:length(f.surv), function(x) mean(fed > f.surv[x]))
  
  # plot(g.cdf, type = "l", col = "blue")
  # lines(f.cdf, col = "red")
  
  ks = max(abs(f.cdf - g.cdf))
  rate = sqrt((ncol(g)*ncol(f)) / (ncol(g) + ncol(f)))
  ks_pval(rate*ks)
}
ks_pval = function(t, n = 20) {
  2*(sum(sapply(1:n, function(x) (-1)^(x-1) * exp(-2*(x^2)*t^2))))
}


ks1.mean = function(f, g) {
  fed = meandepth(f, f)
  ged = meandepth(g, f)
  
  f.surv = rev(c(0, sort(fed)))
  g.cdf = sapply(1:length(f.surv), function(x) mean(ged >= f.surv[x]))
  f.cdf = seq(0, 1, length.out = length(g.cdf))
  
  # plot(g.cdf, type = "l", col = "blue")
  # lines(f.cdf, col = "red")
  
  ks = max(abs(f.cdf - g.cdf))
  rate = sqrt(ncol(g))
  ks_pval(rate*ks)
}

# convert 2d fields into 1D coefficient functions
basis = create_basis(matrix(0, 30, 30), basis.pts = 10)
proj = solve(t(basis) %*% basis) %*% t(basis)

#### SIZE
set.seed(1)

sims = 100
kdist = rep(0, sims)
for(s in 1:sims) {
  tic("Total")
  cat("Simulation ", s, "\n")
  
  nfield = 500
  
  gp1 = gp2d(nfield, mu = 0, pts = 30, l = 10)
  gp2 = gp2d(nfield, mu = 0, pts = 30, l = 10)
  
  # gp1 = sapply(1:nfield, function(x) proj %*% as.vector(gp1[,,x]))
  # gp2 = sapply(1:nfield, function(x) proj %*% as.vector(gp2[,,x]))
  
  gp1 = flatten(gp1)
  gp2 = flatten(gp2)
  
  kdist[s] = ks.mean(gp1, gp2) 
  
  toc()
  cat("\n")
}
mean(kdist < 0.05)
plot(kdist)

plt_funs(gp1, gp2)
