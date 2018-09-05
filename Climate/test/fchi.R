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
ks.dist = function(x, y) {
  nx = ncol(x)
  ny = ncol(y)
  
  x_dist = rep(0, nx)
  y_dist = rep(0, ny)
  
  x_ed = edepth_set(x)
  x_med = x[,which.max(x_ed)]
  
  for(f in 1:nx) {
    x_dist[f] = l2(x[,f], x_med)
  }
  for(f in 1:ny) {
    y_dist[f] = l2(y[,f], x_med)
  }
  
  x_dist = x_dist[x_dist > 0]
  
  ks.test(x_dist, y_dist)$p.value
  # wilcox.test(x_dist, y_dist)$p.value
  # t.test(x_dist, y_dist)$p.value
}

plt_funs = function(fmat) {
  plot(fmat[,1], type = "l", ylim = c(min(fmat), max(fmat)))
  for(i in 2:ncol(fmat)) {
    lines(fmat[,i])
  }
}
plt_2funs = function(f, g) {
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



#### TESTS ####
f = post_flat[1:100,1:10]
g = -sign(f) * abs(f)

plt_2funs(f, g)

ks.dist(f, g)

t = seq(0, 1, length.out = 100)
fq = f_to_srvf(f, t)
gq = f_to_srvf(g, t)

ks.dist(fq, gq)

# differnce of bounds
cr_set = function(f) {
  ed = edepth_set(f)
  
  ed.seq = seq(1, 0, length.out = ncol(f))
  upper = matrix(0, nrow(f), length(ed.seq))
  for(i in 2:ncol(upper)) {
    upper[,i] = central_region(f, ed, ed.seq[i])$upper
  }
  upper[,-1]
}

fcr = cr_set(f)
gcr = cr_set(g)
plt_2funs(fcr, gcr)

upper.dist = function(f, g) {
  fcr = cr_set(f)
  gcr = cr_set(g)
  
  up.dist = rep(0, ncol(fcr))
  for(i in 1:ncol(fcr)) {
    up.dist[i] = l2(fcr[,i], gcr[,i])
  }
  max(up.dist)
}

upper.dist(f, g)


# two sample chi square using the depths to define functional bins
is.contained = function(x, lower, upper) {
  min((x <= upper) & (x >= lower))
}


fchi = function(f, g) {
  ed = edepth_set(f)
  
  ed.seq = rev(sort(ed))[-1]
  
  f.count = rep(0, length(ed.seq))
  g.count = rep(0, length(ed.seq))
  
  for(i in 1:length(ed.seq)) {
    fset = f[, ed >= ed.seq[i]]
    fset = as.matrix(fset)
    # lower
    lower = sapply(1:nrow(fset), function(x) min(fset[x,]))
    
    # upper
    upper = sapply(1:nrow(fset), function(x) max(fset[x,]))
    
    f.count[i] = sum(sapply(1:ncol(f), function(x) is.contained(f[,x], lower, upper))) / ncol(f)
    g.count[i] = sum(sapply(1:ncol(g), function(x) is.contained(g[,x], lower, upper))) / ncol(g)
  }
  
  chi.values = (f.count - g.count)^2 / (f.count + g.count)
  chi = sum(chi.values)
  return(1-pchisq(chi, df = (length(chi.values) - 1)))
}

fchi(f, g)



