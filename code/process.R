rm(list = ls())
gc()

library(extdepth)
library(ncdf4)
library(dplyr)
library(plotly)
library(RcppEigen)
library(RColorBrewer)
library(fields)

##### CONNECT TO DATA #####
# read ensembles and prior ncdf4 objects
nc.ens = nc_open('/Users/trevh/Research/assimilation-cfr/data/tas_ens_da_hydro_r.1000-2000_d.30-Nov-2017.nc')
nc.prior = nc_open('/Users/trevh/Research/assimilation-cfr/data/tas_prior_da_hydro_r.1000-2000_d.30-Nov-2017.nc')

# get data dim
n.lon = nc.ens$dim$lon$len
n.lat = nc.ens$dim$lat$len
n.ens = nc.ens$dim$sub_ens$len



##### CONSTRUCT THE BASIS FUNCTIONS #####
# incase this gets changed do 10 and 15
nlat.centers = 10
nlon.centers = 15
n.center = nlat.centers * nlon.centers

lats = as.vector(nc.ens$dim$lat$vals)
lons = as.vector(nc.ens$dim$lon$vals)
lats.cen = seq(min(lats), max(lats), length.out = nlat.centers)
lons.cen = seq(min(lons), max(lons), length.out = nlon.centers)

# find the radius (1.5 times the minimal distance)
# will either be the distance b/w the first point and the next (horizontally or vertically)
rad = 1.5 * min(abs(lats.cen[1] - lats.cen[2]), abs(lons.cen[1] - lons.cen[2]))

# CHECK THAT THIS LINES UP WITH X
latlon = as.matrix(expand.grid(1:nlat.centers, 1:nlon.centers))
basis = matrix(0, n.lat*n.lon, n.center)

for (i in 1:n.center) {
  
  # dist = (lat - center.lat)^2 + (lon - center.lon)^2
  dist = c(outer((lats - lats.cen[latlon[i,][1]])^2, (lons - lons.cen[latlon[i,][2]])^2, "+"))
  
  # use the bisquare (radial) basis
  basis[,i] = (1 - dist/rad^2)^2 * (sqrt(dist) < rad)
}


t = 700

# extract data from the ncdf4 objects
ens = ncvar_get(nc.ens, attributes(nc.ens$var)$names[1], start = c(1, 1, t, 1), count = c(-1, -1, 1, -1))
prior = ncvar_get(nc.prior, attributes(nc.prior$var)$names[1], start = c(1, 1, t), count = c(-1, -1, 1))

# transpose for intuitive (to me) layout
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



##### FIT BASIS AND FIND ED #####
prior.alpha = coef(fastLmPure(basis, as.vector(prior)))
ens.alpha = sapply(1:n.ens, function(x) coef(fastLmPure(basis, as.vector(ens[,,x]))))

# calculate extremal depth
ed = sapply(1:ncol(ens.alpha), function(x) edepth(ens.alpha[,x], ens.alpha))
central = central_region(ens.alpha, ed, 0.05)


##### DIAGNOSTICS #####

# lines(central[[1]])
# lines(central[[2]])


# # what do the ens coeff look like
# plot(ens.alpha[,1], type = "l")
# for (i in 2:100) {
#   lines(ens.alpha[,i])
# }
# # plus prior
# lines(prior.alpha, col = "red")
# 
# # check the smoothed approximation
# prior.smooth = matrix(fastLm(basis, as.vector(prior))$fitted.values, 96, 144)
# 
# plot_ly(showscale = FALSE) %>% 
#   add_surface(z = ~prior)
# 
# plot_ly(showscale = FALSE) %>% 
#   add_surface(z = ~prior.smooth)



##### PREPROCESS THE FIELDS #####
# eventually it'll be a loooop
# years = 20
# ed = rep(0, years)
# 
# for (t in 1:years) {
#   # extract data from the ncdf4 objects
#   ens = ncvar_get(nc.ens, attributes(nc.ens$var)$names[1], start = c(1, 1, t, 1), count = c(-1, -1, 1, -1))
#   prior = ncvar_get(nc.prior, attributes(nc.prior$var)$names[1], start = c(1, 1, t), count = c(-1, -1, 1))
#   
#   # transpose for intuitive (to me) layout
#   ens = aperm(ens, c(2, 1, 3))
#   prior = t(prior)
#   
#   # remove lat means
#   ens = vapply(1:n.ens, function(x) ens[,,x] - rowMeans(ens[,,x]), FUN.VALUE = matrix(0, nrow = n.lat, ncol = n.lon))
#   prior = prior - rowMeans(prior)
#   
#   # normalize
#   lats = as.vector(nc.ens$dim$lat$vals)
#   latmat = matrix(rep(lats, n.lon), n.lat, n.lon)
#   latmat = sqrt(abs(cos(latmat*pi/180)))
#   
#   ens = vapply(1:n.ens, function(x) ens[,,x]*latmat, FUN.VALUE = matrix(0, nrow = n.lat, ncol = n.lon))
#   prior = prior * latmat
#   
#   
#   
#   ##### FIT BASIS AND FIND ED #####
#   prior.alpha = coef(fastLmPure(basis, as.vector(prior)))
#   ens.alpha = sapply(1:n.ens, function(x) coef(fastLmPure(basis, as.vector(ens[,,x]))))
#   
#   # # calculate extremal depth
#   # ed = sapply(1:ncol(ens.alpha), function(x) edepth(ens.alpha[,x], ens.alpha))
#   # central = central_region(ens.alpha, ed)
#   
#   ed[t] = edepth(prior.alpha, ens.alpha)
# }
# plot(ed)





