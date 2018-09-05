rm(list = ls()); gc()

library(extdepth)
library(tictoc)
library(fda)

# reformat
flatten = function(mat) {
  matrix(mat, prod(dim(mat)[1:2]), dim(mat)[3])
}

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

#### SIMULATION ####
sim_gp = function(fields = 100, mu = 0, l = 30, pts = 30) {
  
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
num.int = function(y, x = seq_len(length(y)) / length(y)) {
  sum(diff(x) * (head(y, -1) + tail(y, -1))) / 2
}
l2 = function(f, g) {
  sqrt(num.int((f-g)^2))
}

depth <- function(g, fmat) {
  
  # Computes the depth values of a function with respect to a set of functions (fmat)
  fn = ncol(fmat)
  depth = rep(0, length(g))
  
  for (row in 1:nrow(fmat)) {
    diff = abs(sum(sign(g[row] - fmat[row,])))
    depth[row] = 1 - (diff / fn)
  }
  
  return(depth)
}
edepth_multi = function(gmat, fmat_sorted) {
  fmat = fmat_sorted
  fdepths = depth_set(fmat)
  gdepths = apply(gmat, 2, function(x) depth(x, fmat))
  
  # get the allowed r values (for calculating the dCDF)
  rvals = sort(unique(c(gdepths, fdepths)))
  
  ged = rep(1, ncol(gmat))
  for(g in 1:ncol(gmat)) {
    for(f in 1:ncol(fmat)) {
      for(r in rvals) {
        dg = sum(gdepths[,g] <= r)
        df = sum(fdepths[,f] <= r)
        if(dg != df) {
          break;
        }
      }
      if (dg > df) {
        ged[g] = (f-1) / ncol(gmat)
        break;
      }
    }
  }
  ged
}
ks_pval = function(t, n = 10) {
  2*(sum(sapply(1:n, function(x) (-1)^(x-1) * exp(-2*(x^2)*t))))
}
dks = function(f, g) {
  
  f = edepth_sort(f)
  fed = (1:ncol(f)) / ncol(f)
  ged = edepth_multi(g, f)
  
  f.surv = rev(c(0, sort(fed)))
  g.cdf = sapply(1:length(f.surv), function(x) mean(ged >= f.surv[x]))
  f.cdf = rev(f.surv)
  
  ks = max(abs(f.cdf - g.cdf))
  ks_pval(sqrt(ncol(g))*ks)
}

edepth_multi(post_flat, prior_flat)


# get args
batch = 10
sims = 50
set.seed(batch + 1023)

ks = rep(0, sims)

for(s in 1:sims) {
  tic("Total")
  cat("Simulation ", s, "\n")
  
  ### Generate Fields
  tic("Generating fields")
  
  ## generate
  nfield = 200
  fields = sim_gp(nfield, mu = 0, l = 45)
  
  ## smoothe
  basis = create_basis(matrix(0, 30, 30), basis.pts = 10)
  proj = solve(t(basis) %*% basis) %*% t(basis)
  
  fields = vapply(1:nfield, function(x) basis %*% (proj %*% as.vector(fields[,,x])),
                  FUN.VALUE = matrix(0, 30, 30))
  
  ## demean
  # fields = vapply(1:nfield, function(x) fields[,,x] - mean(fields[,,x]), FUN.VALUE = matrix(0, 30, 30))
  
  ## split
  ind = sample(1:nfield, nfield/2, replace = F)
  prior = fields[,,ind]
  posterior = fields[,,-ind]
  toc()
  
  ### Pointwise Comparisons
  tic("ED")
  prior_flat = flatten(prior)
  post_flat = flatten(posterior)
  
  ks[s] = dks(prior_flat, post_flat)

  toc()
  cat("\n")
}
plot(ks)
mean(ks < 0.05)
