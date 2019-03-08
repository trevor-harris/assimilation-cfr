rm(list = ls())
gc()

library(ncdf4)
library(reshape2)
library(ggplot2)
library(OpenImageR)

source("../research/assimilation-cfr/code/depth_tests.R")
source("../research/assimilation-cfr/code/depths.R")
source("../research/assimilation-cfr/code/simulation.R")

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
read_era = function(dir, file) {
  cbind(readRDS(paste0(dir, file)), era = as.numeric(strsplit(file, "\\D+")[[1]][-1]))
}
region_plot <- function(field, nc, reg_names, cc = cc, main = "", zlim = c(-max(abs(field)), max(abs(field)))) {
  
  lats = as.vector(nc$dim$lat$vals)
  lons = as.vector(nc$dim$lon$vals)
  dimnames(field) = list(lats, ifelse(lons >= 180, lons - 360, lons))
  
  field.gg = melt(field)
  colnames(field.gg) = c("lat", "lon", "Region")
  field.gg$Region = as.factor(field.gg$Region)
  levels(field.gg$Region) = reg_names
  
  world = map_data("world")
  world = world[world$long <= 178, ]
  
  zlim = c(-max(abs(field)), max(abs(field)))
  
  alphas = ifelse(field.gg$Region == 10000, 1, 0.8)
  
  ggplot() +
    geom_raster(data = field.gg, aes(x=lon, y=lat, fill=Region), interpolate = F, alpha = alphas) +
    geom_polygon(data = world, aes(x=long, y=lat, group=group), fill = NA, color="black") +
    coord_cartesian() +
    scale_fill_manual(values = cc) +
    theme_void() +
    ggtitle(main) +
    theme(legend.position="none") +
    theme(plot.title = element_text(hjust = 0.5))
}

region_plot(mask, nc.prior, reg_names, cc)


save_dir = "../research/assimilation-cfr/paper/results/"

#### Regionlized significant differences
nc.post = nc_open('../research/assimilation-cfr/data/tas_ens_da_hydro_r.1000-2000_d.16-Feb-2018.nc')
nc.prior = nc_open('../research/assimilation-cfr/data/tas_prior_da_hydro_r.1000-2000_d.16-Feb-2018.nc')
prior_ind = read.csv("../research/assimilation-cfr/data/prior_ens.txt", header = F)$V1

# read in the mask file
mask = read.csv("../research/assimilation-cfr/data/mask2.csv", stringsAsFactors = F)[,2:145]
mask = as.matrix(mask)

reg_names = c(c("Arctic Ocean", "Indian Ocean", "Pacific Ocean", "Atlantic Ocean", "Southern Ocean"), 
              c("Antarctica", "South America", "North America", "Africa", "Europe", "Asia", "Australia", "Borders"))

cc = c("#F8766D","#E18A00","#BE9C00","#8CAB00","#24B700","#00BE70",
       "#00C1AB","#00BBDA","#00ACFC","#8B93FF","#D575FE","#F962DD",
       "#000000")

region_plot(mask, nc.prior, reg_names, cc)
ggsave(paste0(save_dir, "regions.png"), width = 9, height = 6)
