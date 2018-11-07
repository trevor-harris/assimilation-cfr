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


# reformat
flatten = function(mat) {
  matrix(mat, prod(dim(mat)[1:2]), dim(mat)[3])
}

# get args
args = commandArgs(TRUE)
batch = as.double(args[1])
batch = 1

sims = 10
perms = 1000

beta_obs = array(0, c(30, 30, sims))
ed_lower = array(0, c(30, 30, sims))
ed_upper = array(0, c(30, 30, sims))

ed = rep(0, sims)
fdr = array(0, c(30, 30, sims))

set.seed(batch + 1023)

for(s in 1:1) {
  # s = 3
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
  posterior = fields[,,-ind]+1.5
  toc()

  ### Pointwise Comparisons
  tic("ED")
  prior_flat = flatten(prior)
  post_flat = flatten(posterior)
  # post_eds = sapply(1:(nfield/2), function(x) edepth(post_flat[,x], prior_flat))
  # post_eds = sapply(1:(nfield/2), function(x) edepth_set(cbind(post_flat[,x], prior_flat))[1])
  # fdr_eds = p.adjust(post_eds, "BH")

  ### MEAN
  post_mean = rowMeans(post_flat) / (apply(post_flat, 1, sd))
  post_eds = edepth(post_mean, prior_flat)

  cr = central_region(prior_flat, edepth_set(prior_flat))
  upper = cr$upper
  lower = cr$lower
  post_eds

  toc()

  plot(post_mean, type = "l", ylim = c(-3, 3), col = "red")
  lines(upper)
  lines(lower)

  cat("\n")
}

######################################################


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


# reformat
flatten = function(mat) {
  matrix(mat, prod(dim(mat)[1:2]), dim(mat)[3])
}

# get args
args = commandArgs(TRUE)
batch = as.double(args[1])
batch = 1

sims = 10
perms = 1000

beta_obs = array(0, c(30, 30, sims))
ed_lower = array(0, c(30, 30, sims))
ed_upper = array(0, c(30, 30, sims))

ed = rep(0, sims)
fdr = array(0, c(30, 30, sims))

set.seed(batch + 1023)

r = 0.9
for(s in 1:1) {
  # s = 3
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
  posterior = fields[,,-ind] + r
  toc()

  prior_flat = flatten(prior)
  post_flat = flatten(posterior)

  # prior_time = vapply(1:10, function(x) prior_flat + matrix(rnorm(900*100), 900, 100), FUN.VALUE = array(0, c(900, 100)))
  # post_time = vapply(1:10, function(x) post_flat + matrix(rnorm(900*100), 900, 100), FUN.VALUE = array(0, c(900, 100)))

  prior_time = vapply(1:10, function(x) prior_flat + (x/10)^3, FUN.VALUE = array(0, c(900, 100)))
  post_time = vapply(1:10, function(x) post_flat + (x/10)^3, FUN.VALUE = array(0, c(900, 100)))

  prior_time = aperm(prior_time, c(1, 3, 2))
  post_time = aperm(post_time, c(1, 3, 2))
  ### Pointwise Comparisons
  tic("ED")

  prior_flat = flatten(prior_time)
  post_flat = flatten(post_time)

  ### MEAN
  post_mean = rowMeans(post_flat)
  post_eds = edepth(post_mean, prior_flat)

  cr = central_region(prior_flat, edepth_set(prior_flat))
  upper = cr$upper
  lower = cr$lower
  post_eds

  toc()

  plot(post_mean, type = "l", ylim = c(-4, 4), col = "red")
  lines(upper)
  lines(lower)

  cat("\n")
}



######################################################



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


# reformat
flatten = function(mat) {
  matrix(mat, prod(dim(mat)[1:2]), dim(mat)[3])
}

# get args
args = commandArgs(TRUE)
batch = as.double(args[1])
batch = 1

sims = 10
perms = 1000

beta_obs = array(0, c(30, 30, sims))
ed_lower = array(0, c(30, 30, sims))
ed_upper = array(0, c(30, 30, sims))

ed = rep(0, sims)
fdr = array(0, c(30, 30, sims))

set.seed(batch + 1023)

r = 0.9
for(s in 1:1) {
  # s = 3
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
  posterior = fields[,,-ind] + r
  toc()

  prior_flat = flatten(prior)
  post_flat = flatten(posterior)

  prior_time = vapply(1:10, function(x) prior_flat + matrix(rnorm(900*100, x/10), 900, 100), FUN.VALUE = array(0, c(900, 100)))
  post_time = vapply(1:10, function(x) post_flat + matrix(rnorm(900*100, x/10), 900, 100), FUN.VALUE = array(0, c(900, 100)))

  # prior_time = vapply(1:10, function(x) prior_flat + (x/10)^3, FUN.VALUE = array(0, c(900, 100)))
  # post_time = vapply(1:10, function(x) post_flat + (x/10)^3, FUN.VALUE = array(0, c(900, 100)))

  prior_time = aperm(prior_time, c(1, 3, 2))
  post_time = aperm(post_time, c(1, 3, 2))
  ### Pointwise Comparisons
  tic("ED")

  prior_flat = flatten(prior_time)
  post_flat = flatten(post_time)

  ### MEAN
  post_mean = rowMeans(post_flat)
  post_eds = edepth(post_mean, prior_flat)

  cr = central_region(prior_flat, edepth_set(prior_flat))
  upper = cr$upper
  lower = cr$lower
  post_eds

  toc()

  plot(post_mean, type = "l", ylim = c(-4, 4), col = "red")
  lines(upper)
  lines(lower)

  cat("\n")
}

