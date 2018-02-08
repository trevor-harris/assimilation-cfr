rm(list = ls())
gc()

library(ncdf4)
library(dplyr)
library(plotly)

# eventually it'll be a loooop
t = 1

# read ensembles and prior ncdf4 objects
nc.ens = nc_open('/Users/trevh/Research/assimilation-cfr/data/tas_ens_da_hydro_r.1000-2000_d.30-Nov-2017.nc')
nc.prior = nc_open('/Users/trevh/Research/assimilation-cfr/data/tas_prior_da_hydro_r.1000-2000_d.30-Nov-2017.nc')

# get data dim
n.lon = nc.ens$dim$lon$len
n.lat = nc.ens$dim$lat$len
n.ens = nc.ens$dim$sub_ens$len

# extract data from the ncdf4 objects
ens = ncvar_get(nc.ens, attributes(nc.ens$var)$names[1], count = c(-1, -1, t, -1))
prior = ncvar_get(nc.prior, attributes(nc.prior$var)$names[1], count = c(-1, -1, t))

# transpose for intuitive layout
ens = aperm(ens, c(2, 1, 3))
prior = t(prior)

# remove lat means
ens = vapply(1:n.ens, function(x) ens[,,x] - rowMeans(ens[,,x]), FUN.VALUE = matrix(0, nrow = n.lat, ncol = n.lon))
prior = prior - rowMeans(prior)

# normalize
lats = as.vector(nc.ens$dim$lat$vals)
latmat = matrix(rep(lats, n.lon), n.lat, n.lon)
latmat = sqrt(abs(cos(latmat*pi/180)))

ens = vapply(1:n.ens, function(x) ens[,,x]*latmat, FUN.VALUE = matrix(0, nrow = n.lat, ncol = n.lon))
prior = prior * latmat

# represent in the bisquare basis
l2.norm = function(x1, x2, y1, y2) {
  return (sqrt((x1 + x2)^2 + (y1 + y2)^2))
}

pts = as.integer(seq(from = 1, to = nrow(prior) * ncol(prior), length.out = 200))
pts

plot_ly(showscale = FALSE) %>% 
  add_surface(z = ~prior) %>%
  add_surface(z = ~ens[,,1])
