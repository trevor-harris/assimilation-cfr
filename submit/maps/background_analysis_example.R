rm(list = ls())
gc()

library(ncdf4)
library(tictoc)
library(ggplot2)
library(dplyr)
library(reshape2)

# set to the top level folder
setwd("/Users/trevorh2/research/assimilation-cfr/submit/")


source("method/depth_tests.R")
source("method/depths.R")
source("method/simulation.R")

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

#### Regionlized significant differences
# prior
nc.prior = nc_open('data/tas_prior_da_hydro_r.1000-2000_d.16-Feb-2018.nc')
prior_ind = read.csv("data/prior_ens.txt", header = F)$V1

prior.all = prep_prior(nc.prior)
prior.all = prior.all[,,prior_ind]
prior.mean = apply(prior.all, c(1, 2), mean)

# post
nc.post = nc_open('data/tas_ens_da_hydro_r.1000-2000_d.16-Feb-2018.nc')
post.all = prep_post(nc.post, 1)


prior = prior.all[,,70] - prior.mean
post = post.all[,,70] - prior.mean

lats = as.vector(nc.prior$dim$lat$vals)
lons = as.vector(nc.prior$dim$lon$vals)

dimnames(prior) = list(lats, ifelse(lons >= 180, lons - 360, lons))
dimnames(post) = list(lats, ifelse(lons >= 180, lons - 360, lons))

prior = melt(prior)
post = melt(post)

colnames(prior) = c("lat", "lon", "val")
colnames(post) = c("lat", "lon", "val")

prior[["dist"]] = "Background ensemble member 70"
post[["dist"]] = "Analysis ensemble member 70 (850 CE)"

world = map_data("world")
world = world[world$long <= 178, ]

combined = rbind(prior, post)

combined$dist = factor(combined$dist, levels = c(
  "Background ensemble member 70", "Analysis ensemble member 70 (850 CE)"
))

ggplot() +
  geom_raster(data = combined, aes(x=lon, y=lat, fill=val)) +
  geom_polygon(data = world, aes(x=long, y=lat, group=group), fill = NA, color="black") +
  coord_cartesian() +
  theme_void() +
  scale_fill_distiller(palette = "RdBu",
                       limits = c(-2.5, 4),
                       name = "Temp.\nAnomaly (C)") +
  theme(plot.title = element_text(hjust = 0.5),
        strip.text = element_text(size = 16,
                                  margin=margin(b=2))) +
  facet_grid(cols = vars(dist))

# ggsave("maps/background_analysis_anomaly.png", width = 10, height = 4)
