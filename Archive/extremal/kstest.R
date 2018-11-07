rm(list = ls())

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

permute_fields = function(prior, posterior) {
  # permutes 2 samples of 2D regionalized functions
  # prior.split = 4D array. dim1 = lat, dim2 = lon, dim3 = region no, dim4 = sample number
  # post.split = same as prior.split
  # returns a list = [permuted.prior, permuted.post]
  
  nlat = dim(prior)[1]
  nlon = dim(prior)[2]
  ens = dim(prior)[3]
  
  prior.new = prior
  post.new = posterior
  
  prior.ind = sample(c(0, 1), size = ens, replace = TRUE)
  
  for (e in 1:ens) {
    if (prior.ind[e] == 1) {
      prior.new[,,e] = posterior[,,e]
      post.new[,,e] = prior[,,e]
    } else {
      prior.new[,,e] = prior[,,e]
      post.new[,,e] = posterior[,,e]
    }
  }
  
  return(list(prior.new, post.new))
}

num.int = function(y, x = seq_len(length(y)) / length(y)) {
  sum(diff(x) * (head(y, -1) + tail(y, -1))) / 2
}

l2 = function(f, g) {
  sqrt(num.int((f-g)^2))
}

# reformat
flatten = function(mat) {
  matrix(mat, prod(dim(mat)[1:2]), dim(mat)[3])
}

# get args
args = commandArgs(TRUE)
batch = as.double(args[1])
batch = 10

sims = 500

beta_obs = array(0, c(30, 30, sims))
ed_lower = array(0, c(30, 30, sims))
ed_upper = array(0, c(30, 30, sims))

ed = rep(0, sims)
fdr = array(0, c(30, 30, sims))

set.seed(batch + 1023)

for(s in 1:sims) {
  # s = 1
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
  fields = vapply(1:nfield, function(x) fields[,,x] - mean(fields[,,x]), FUN.VALUE = matrix(0, 30, 30))
  
  ## split
  ind = sample(1:nfield, nfield/2, replace = F)
  prior = fields[,,ind]
  posterior = fields[,,-ind]
  toc()
  
  ### Pointwise Comparisons
  tic("ED")
  prior_flat = flatten(prior)
  post_flat = flatten(posterior)
  
  ### MEAN
  prior_eds = edepth_set(prior_flat)
  prior_med = prior_flat[,which.max(prior_eds)]
  
  prior_dist = rep(0, 100)
  post_dist = rep(0, 100)
  for(f in 1:100) {
    prior_dist[f] = l2(prior_flat[,f], prior_med)
    post_dist[f] = l2(post_flat[,f], prior_med)
  }
  
  toc()
  
  prior_dist = prior_dist[prior_dist > 0]
  ed[s] = ks.test(prior_dist, post_dist)$p.value
  cat("\n")
}

plot(ed)
mean(ed < 0.05)



