rm(list = ls())

library(extdepth)
library(tictoc)

#### FUNCTIONS ####
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
ed_mwu = function(fmat, gmat) {
  fcol = ncol(fmat)
  gcol = ncol(gmat)
  
  fed = 1-edepth_set(fmat) 
  ged = 1-sapply(1:gcol, function(g) edepth(gmat[,g], cbind(gmat[,g], fmat)))
  combined_ranks = rank(c(fed, ged))
  
  U_f = sum(combined_ranks[1:fcol]) - fcol*(fcol + 1)/2
  U_g = sum(combined_ranks[-(1:fcol)]) - gcol*(gcol + 1)/2
  
  mwu = min(U_f, U_g)
  mwu.mean = (fcol * gcol) / 2
  mwu.sd = sqrt((fcol * gcol) * (fcol + gcol + 1) / 12)
  
  mwu.z = abs((mwu - mwu.mean) / mwu.sd)
  return(1-pnorm(mwu.z))
}

# reformat
flatten = function(mat) {
  matrix(mat, prod(dim(mat)[1:2]), dim(mat)[3])
}


#### SIMULATION ####
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

prior_flat = flatten(prior)
post_flat = flatten(posterior)


test = post_flat[1:100,1:20]

test.ed = edepth_set(test)

test.rnk = sapply(1:20, function(f) edepth(test[,f], test, depth_function = "rank"))
edepth(test[,1], test[,c(1,2)], depth_function = "rank")

small = test[,1:2]
plot(small[,1], type = "l", col = "red")
lines(small[,2])

small.d = depth_set(small)
plot(small.d[,1], type = "l", col = "red")
lines(small.d[,2])
