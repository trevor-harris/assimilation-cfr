rm(list = ls())
gc()

library(extdepth)
library(ncdf4)
library(dplyr)
library(reshape2)

# prepare data
prep_prior = function(nc.prior) {
  
  n.lon = nc.prior$dim$lon$len
  n.lat = nc.prior$dim$lat$len
  n.ens = nc.prior$dim$time2$len
  
  # extract data from the ncdf4 objects
  prior = ncvar_get(nc.prior, attributes(nc.prior$var)$names[1], start = c(1, 1, 1), count = c(-1, -1, -1))
  
  # transpose for intuitive (to me) layout
  prior = aperm(prior, c(2, 1, 3))
  
  # remove lat means
  # prior = vapply(1:n.ens, function(x) prior[,,x] - rowMeans(prior[,,x]), FUN.VALUE = matrix(0, nrow = n.lat, ncol = n.lon))
  
  # normalize
  lats = as.vector(nc.prior$dim$lat$vals)
  latmat = matrix(rep(lats, n.lon), n.lat, n.lon)
  latmat = sqrt(abs(cos(latmat*pi/180)))
  
  prior = vapply(1:n.ens, function(x) prior[,,x]*latmat, FUN.VALUE = matrix(0, nrow = n.lat, ncol = n.lon))
  
  return(prior)
}
prep_post = function(nc.post, t) {
  
  n.lon = nc.post$dim$lon$len
  n.lat = nc.post$dim$lat$len
  n.ens = nc.post$dim$sub_ens$len
  
  # extract data from the ncdf4 objects
  ens = ncvar_get(nc.post, attributes(nc.post$var)$names[1], start = c(1, 1, t, 1), count = c(-1, -1, 1, -1))
  
  # transpose for intuitive (to me) layout
  ens = aperm(ens, c(2, 1, 3))
  
  # remove lat means
  # ens = vapply(1:n.ens, function(x) ens[,,x] - rowMeans(ens[,,x]), FUN.VALUE = matrix(0, nrow = n.lat, ncol = n.lon))
  
  # normalize
  lats = as.vector(nc.post$dim$lat$vals)
  latmat = matrix(rep(lats, n.lon), n.lat, n.lon)
  latmat = sqrt(abs(cos(latmat*pi/180)))
  
  ens = vapply(1:n.ens, function(x) ens[,,x]*latmat, FUN.VALUE = matrix(0, nrow = n.lat, ncol = n.lon))
  
  return(ens)
}
create_basis = function(nlat.centers, nlon.centers, nc) {
  
  n.center = nlat.centers * nlon.centers
  
  lats = as.vector(nc$dim$lat$vals)
  lons = as.vector(nc$dim$lon$vals)
  n.lon = nc$dim$lon$len
  n.lat = nc$dim$lat$len
  
  # lats.cen = seq(min(lats), max(lats), length.out = nlat.centers)
  lons.cen = seq(min(lons), max(lons), length.out = nlon.centers)
  
  # cosine normalize to add more basis towards tropics and away from poles
  lats.pi = cos((lats / max(lats)) * (pi/2))
  lats.cen = seq(min(-lats.pi), max(lats.pi), length.out = nlat.centers)
  lats.cen = (acos(lats.cen) / (pi/2) * max(lats)) - 90
  lats.cen = -lats.cen
  
  
  # find the radius (1.5 times the minimal distance)
  # will either be the distance b/w the first point and the next (horizontally or vertically)
  rad = 1.5 * min(abs(lats.cen[1] - lats.cen[2]), abs(lons.cen[1] - lons.cen[2]))
  
  latlon = as.matrix(expand.grid(1:nlat.centers, 1:nlon.centers))
  basis = matrix(0, n.lat*n.lon, n.center)
  
  for (i in 1:n.center) {
    
    # dist = (lat - center.lat)^2 + (lon - center.lon)^2
    dist = c(outer((lats - lats.cen[latlon[i,][1]])^2, (lons - lons.cen[latlon[i,][2]])^2, "+"))
    
    # use the bisquare (radial) basis
    basis[,i] = (1 - dist/rad^2)^2 * (sqrt(dist) < rad)
  }
  
  return(basis)
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

# read ensembles and prior ncdf4 objects
nc.post = nc_open('/Users/trevh/research/assimilation-cfr/data/tas_ens_da_hydro_r.1000-2000_d.16-Feb-2018.nc')
nc.prior = nc_open('/Users/trevh/research/assimilation-cfr/data/tas_prior_da_hydro_r.1000-2000_d.16-Feb-2018.nc')

prior = prep_prior(nc.prior)
post = prep_post(nc.post, 100)

# convert matrix representation into vector rep
# prior = flatten(prior)
# posterior = flatten(posterior)

# convert field representation into basis coefficient representation
prior_basis = create_basis(40, 40, nc.prior)
post_basis = create_basis(40, 40, nc.post)

prior_proj = solve(t(prior_basis) %*% prior_basis) %*% t(prior_basis)
post_proj = solve(t(post_basis) %*% post_basis) %*% t(post_basis)

prior = sapply(1:dim(prior)[3], function(x) prior_proj %*% as.vector(prior[,,x]))
post = sapply(1:dim(post)[3], function(x) post_proj %*% as.vector(post[,,x]))

# test
ks.mean(prior, post)
