rm(list = ls())
gc()

library(ncdf4)
library(dplyr)
library(plotly)
library(mgcv)
library(irlba)


#### PRIOR PREPROCESS
nc = nc_open('/Users/trevh/Research/assimilation-cfr/data/tas_prior_da_hydro_r.1000-2000_d.30-Nov-2017.nc')

prior = ncvar_get(nc, attributes(nc$var)$names[1], count = c(-1, -1, -1))
prior = lapply(1:dim(prior)[3], function(x) t(as.matrix(prior[ , , x])))

lats = as.vector(nc$dim$lat$vals)
nlon = nc$dim$lon$len
nlat = nc$dim$lat$len

latmat = matrix(rep(lats, nlon), nlat, nlon)
latmat = sqrt(abs(cos(latmat*pi/180)))

for (c in 1:length(prior)) {
  prior[[c]] = data.frame(x = rep(seq_len(nlon), each = nlat),
                          y = rep(seq_len(nlat), times = nlon),
                          z = c(prior[[c]]))
  # gam smoother
  mod = gam(z ~ te(x, y, k = c(10,12)), data = prior[[c]])$fitted
  
  # return to matrix format
  prior[[c]] = matrix(mod, nlat, nlon)
  
  # remove lat mean
  prior.t = prior[[c]] - rowMeans(prior[[c]])
  
  # cosine normalize
  prior[[c]] = prior.t*latmat
  
}

# cleanup
rm(c, lats, nlat, nlon, latmat, prior.t, mod)

# check em
plot_ly(showscale = FALSE) %>% 
  add_surface(z = ~prior[[1]])


#### FIND SVD BASIS

# convert to big matrix
prior.mat = matrix(0, length(prior), length(prior[[1]]))
for (t in 1:length(prior)) {
  prior.mat[t, ] = as.vector(unlist(prior[[t]]))
}

# compute full covariance of the functions
prior.cov = cov(prior.mat)

# partial SVD to find the firt nv eigenvalues and eigenfunctions
prior.eig = irlba(prior.cov, nv=50)

