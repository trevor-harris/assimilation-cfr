##### FIND OPTIMAL NUMBER OF BASIS FUNCTIONS #####
rm(list = ls())
gc()

library(ncdf4)
library(dplyr)
library(plotly)

##### CONNECT TO DATA #####
# read ensembles and prior ncdf4 objects
nc.prior = nc_open('/Users/trevh/Research/assimilation-cfr/data/tas_prior_da_hydro_r.1000-2000_d.30-Nov-2017.nc')

# get data dim
n.lon = nc.prior$dim$lon$len
n.lat = nc.prior$dim$lat$len
n.time = nc.prior$dim$time2$len

##### AVERAGE OVER ALL PRIOR FIELDS #####
prior.avg = matrix(0, n.lat, n.lon)
for (t in 1:n.time) {
  # extract data from the ncdf4 objects
  prior = ncvar_get(nc.prior, attributes(nc.prior$var)$names[1], start = c(1, 1, t), count = c(-1, -1, 1))
  
  # transpose for intuitive (to me) layout
  prior = t(prior)
  
  # remove lat means
  prior = prior - rowMeans(prior)
  
  # normalize
  lats = as.vector(nc.prior$dim$lat$vals)
  latmat = matrix(rep(lats, n.lon), n.lat, n.lon)
  latmat = sqrt(abs(cos(latmat*pi/180)))
  
  prior = prior * latmat
  prior.avg = prior.avg + prior
}
prior.avg = prior.avg / n.time



##### FIND OPTIMAL NUMBER OF BASIS FUNCTIONS #####
lon.range = 10:30
lat.range = 10:30

MSE = matrix(0, length(lon.range), length(lat.range))
runtime = matrix(0, length(lon.range), length(lat.range))
for (nlo in lon.range) {
  for (nla in lat.range) {
    start.time = Sys.time()
    nlat.centers = nlo
    nlon.centers = nla
    n.center = nlat.centers * nlon.centers
    
    lats = as.vector(nc.prior$dim$lat$vals)
    lons = as.vector(nc.prior$dim$lon$vals)
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
    
    # check the smoothed approximation
    prior.smooth = matrix(fastLmPure(basis, as.vector(prior.avg))$fitted.values, 96, 144)
 
    MSE[nlo-min(lon.range)+1, nla-min(lat.range)+1] = sum((prior.avg-prior.smooth)^2)
    runtime[nlo-min(lon.range)+1, nla-min(lat.range)+1] = Sys.time() - start.time
    
  }
}

# MSE over lat conditional on longittude
plot(MSE[1,], type = "l", ylim = c(20000, 100000))
for(i in 2:nrow(MSE)) {
  lines(MSE[i,])
}

# MSE surface
plot_ly(showscale = FALSE) %>% 
  add_surface(z = ~MSE) %>%
  layout(
    scene = list(
      xaxis = list(title = "lon"),
      yaxis = list(title = "lat"),
      zaxis = list(title = "MSE")
    )
  )

# runtime surface
plot_ly(showscale = FALSE) %>% 
  add_surface(z = ~runtime) %>%
  layout(
    scene = list(
      xaxis = list(title = "lon"),
      yaxis = list(title = "lat"),
      zaxis = list(title = "Run Time")
    )
  )

# for the top 20 lowest MSE what are the lat and lon?
for(i in 0:20) {
  minN.mse = sort(MSE, decreasing = TRUE)[length(MSE) - i]
  minN.ind = which(MSE == minN.mse, arr.ind = TRUE)
  cat(" nlon = ", minN.ind[,2] + min(lon.range)-1, " nlat = ", minN.ind[,1] + min(lat.range)-1,  minN.mse, "\n")
}

min.mse = which(MSE == min(MSE), arr.ind = TRUE)
cat(" nlon = ", min.mse[,2] + min(lon.range)-1, "\n", "nlat = ", min.mse[,1] + min(lat.range)-1)
