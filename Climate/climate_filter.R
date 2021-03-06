rm(list = ls())
gc()

library(extdepth)
library(ncdf4)
library(dplyr)
library(reshape2)
library(OpenImageR)

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
ks.da = function(f, g, ged) {
  fed = meandepth(f, f)
  ged = ged
  
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
ks_pval = function(t, n = 20) {
  2*(sum(sapply(1:n, function(x) (-1)^(x-1) * exp(-2*(x^2)*t^2))))
}

# read ensembles and prior ncdf4 objects
# nc.post = nc_open('/Users/trevh/research/assimilation-cfr/data/tas_ens_da_hydro_r.1000-2000_d.16-Feb-2018.nc')
# nc.prior = nc_open('/Users/trevh/research/assimilation-cfr/data/tas_prior_da_hydro_r.1000-2000_d.16-Feb-2018.nc')


nc.post = nc_open('../research/climate_data/tas_ens_da_hydro_r.1000-2000_d.16-Feb-2018.nc')
nc.prior = nc_open('../research/climate_data/tas_prior_da_hydro_r.1000-2000_d.16-Feb-2018.nc')

prior = prep_prior(nc.prior)

# times = as.integer(seq(1, 998, length.out = 20))
times = 1:998
ksp = rep(0, length(times))

prior.ed = meandepth(prior, prior)

for(t in 1:length(times)) {
  post = prep_post(nc.post, times[t])
  
  # apply low pass filter to remove "noise" (really just to speed it up)
  prior.t = sapply(1:dim(prior)[3], function(x) down_sample_image(prior[,,x], 4))
  post.t = sapply(1:dim(post)[3], function(x) down_sample_image(post[,,x], 4))
  
  # test
  ksp[t] = ks.da(post.t, prior.t, prior.ed)
  cat(times[t], ":", ksp[t], "\n")
}

ksp.adj = p.adjust(ksp, n = 998)

plot(ksp.adj)

summary(ksp.adj)

write.csv(ksp, "../research/assimilation-cfr/cfr/pvals")
write.csv(ksp.adj, "../research/assimilation-cfr/cfr/adjusted_pvals")
