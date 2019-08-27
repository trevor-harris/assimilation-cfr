
rm(list = ls())
gc()

years = 851:1848
ens = 100

library(extdepth)
library(ncdf4)
library(dplyr)
library(reshape2)
library(ggplot2)
library(OpenImageR)
library(future)
library(future.apply)


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

# plots
remove_cr = function(cr, gmat, downsamp=1) {
  lower = cr$lower
  upper = cr$upper
  out = rowMeans((gmat - lower)*(lower > gmat) + (gmat - upper)*(upper < gmat))
  
  matrix(out, 96/downsamp, 144/downsamp)
}
field_plot2 <- function(field, nc) {
  
  lats = as.vector(nc$dim$lat$vals)
  lons = as.vector(nc$dim$lon$vals)
  dimnames(field) = list(lats, ifelse(lons >= 180, lons - 360, lons))
  
  field.gg = melt(field)
  colnames(field.gg) = c("lat", "lon", "Temp")
  
  world = map_data("world")
  world = world[world$long <= 178, ]
  
  zlim = c(-max(abs(field)), max(abs(field)))
  
  plt = ggplot() +
    geom_raster(data = field.gg, aes(x=lon, y=lat, fill=Temp), interpolate = F) +
    geom_polygon(data = world, aes(x=long, y=lat, group=group), fill = NA, color="black") +
    coord_cartesian() +
    scale_fill_distiller(palette = "RdBu") +
    theme_void()
  # labs(fill = "Anomaly")
  # theme(legend.position="none")
  
  print(plt)
}


##### Actual data
nc.post = nc_open('../research/climate_data/tas_ens_da_hydro_r.1000-2000_d.16-Feb-2018.nc')
nc.prior = nc_open('../research/climate_data/tas_prior_da_hydro_r.1000-2000_d.16-Feb-2018.nc')

prior = prep_prior(nc.prior)
post = prep_post(nc.post, 1)

mu = matrix(0, nrow(prior), ncol(prior))
for(i in 1:nrow(prior)) {
  for(j in 1:ncol(prior)) {
    mu[i, j] = mean(prior[i,j,])
  }
}

# prior.anom = prior[,,1] - mu
# post.anom = post[,,2] - mu

##### ACTUAL CLIMATE ##### 
# field_plot2(prior.anom, nc.post)
# ggsave("../research/assimilation-cfr/paper/misc/prior_anom.png", width = 5, height = 3.2)
# 
# 
# field_plot2(post.anom, nc.post)
# ggsave("../research/assimilation-cfr/paper/misc/post_anom.png", width = 5, height = 3.2)



nc = nc.post
prior.anom = prior[,,2]
post.anom = post[,,2]


lats = as.vector(nc$dim$lat$vals)
lons = as.vector(nc$dim$lon$vals)
dimnames(prior.anom) = list(lats, ifelse(lons >= 180, lons - 360, lons))
dimnames(post.anom) = list(lats, ifelse(lons >= 180, lons - 360, lons))

prior.gg = melt(prior.anom)
prior.gg[["DA"]] = "Background"
colnames(prior.gg) = c("lat", "lon", "Temp", "DA")

post.gg = melt(post.anom)
post.gg[["DA"]] = "Analysis"
colnames(post.gg) = c("lat", "lon", "Temp", "DA")

field.gg = rbind(prior.gg, post.gg)
field.gg[["DA"]] = as.factor(field.gg[["DA"]])
levels(field.gg[["DA"]]) = levels(field.gg[["DA"]])[c(2, 1)]

world = map_data("world")
world = world[world$long <= 178, ]

ggplot() +
  geom_raster(data = field.gg, aes(x=lon, y=lat, fill=Temp), interpolate = F) +
  geom_polygon(data = world, aes(x=long, y=lat, group=group), fill = NA, color="black") +
  coord_cartesian() +
  scale_fill_distiller(palette = "RdBu") +
  theme_void() +
  theme(strip.text = element_text(size = 14)) +
  labs(fill="Temp. (C)") +
  facet_grid(. ~ DA)

ggsave("../research/assimilation-cfr/paper/misc/background_analysis.png", width = 10, height = 3.2)



lats = as.vector(nc$dim$lat$vals)
lons = as.vector(nc$dim$lon$vals)
dimnames(prior.anom) = list(lats, ifelse(lons >= 180, lons - 360, lons))
dimnames(post.anom) = list(lats, ifelse(lons >= 180, lons - 360, lons))

prior.gg = melt(prior.anom - mu)
prior.gg[["DA"]] = "Background"
colnames(prior.gg) = c("lat", "lon", "Temp", "DA")

post.gg = melt(post.anom - mu)
post.gg[["DA"]] = "Analysis"
colnames(post.gg) = c("lat", "lon", "Temp", "DA")

field.gg = rbind(prior.gg, post.gg)
field.gg[["DA"]] = as.factor(field.gg[["DA"]])
levels(field.gg[["DA"]]) = levels(field.gg[["DA"]])[c(2, 1)]

world = map_data("world")
world = world[world$long <= 178, ]

min(field.gg$Temp)

ggplot() +
  geom_raster(data = field.gg, aes(x=lon, y=lat, fill=Temp), interpolate = F) +
  geom_polygon(data = world, aes(x=long, y=lat, group=group), fill = NA, color="black") +
  coord_cartesian() +
  scale_fill_distiller(palette = "RdBu", limits =c(-2.5, 4)) +
  theme_void() +
  theme(strip.text = element_text(size = 14)) +
  labs(fill="Temp.\nAnomaly (C)") +
  facet_grid(. ~ DA) 

ggsave("../research/assimilation-cfr/paper/misc/background_analysis_anomaly.png", width = 10.2, height = 3.2)

