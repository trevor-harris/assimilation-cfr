rm(list = ls())
gc()

library(ncdf4)
library(dplyr)
library(reshape2)
library(ggplot2)

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


region_plot <- function(field, nc, main = "", zlim = c(-max(abs(field)), max(abs(field)))) {
  
  lats = rev(as.vector(nc$dim$lat$vals))
  lons = as.vector(nc$dim$lon$vals)
  dimnames(field) = list(lats, sort(ifelse(lons >= 180, lons - 360, lons)))
  
  field.gg = melt(field)
  colnames(field.gg) = c("lat", "lon", "Region")
  field.gg[["Region"]] = as.factor(field.gg[["Region"]])
  
  lcuts = field.gg %>% group_by(lon, Region) %>% summarize(lcut = max(lat))
  indy = tail(unique(lcuts$lcut), -1)
  indy = indy - median(indy)
  
  lcuts = field.gg %>% group_by(lat, Region) %>% summarize(lcut = max(lon))
  indx = head(unique(lcuts$lcut), -1) + 1
  
  world = map_data("world")
  world = world[world$long <= 178, ]
  
  zlim = c(-max(abs(field)), max(abs(field)))
  
  ggplot() +
    geom_raster(data = field.gg, aes(x=lon, y=lat, fill=Region), alpha = 0.4) +
    geom_polygon(data = world, aes(x=long, y=lat, group=group), fill = NA, color="black") +
    geom_hline(yintercept = indy, color = "White") +
    geom_vline(xintercept = indx, color = "white") +
    coord_cartesian() +
    theme_void() +
    ggtitle(main) +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(legend.position="none")
}

region_plot_simple <- function(field, nc, main = "", zlim = c(-max(abs(field)), max(abs(field)))) {
  
  lats = rev(as.vector(nc$dim$lat$vals))
  lons = as.vector(nc$dim$lon$vals)
  dimnames(field) = list(lats, sort(ifelse(lons >= 180, lons - 360, lons)))
  
  field.gg = melt(field)
  colnames(field.gg) = c("lat", "lon", "Region")
  
  lcuts = field.gg %>% group_by(lon, Region) %>% summarize(lcut = max(lat))
  indy = tail(unique(lcuts$lcut), -1)
  indy = indy - median(indy)
  
  lcuts = field.gg %>% group_by(lat, Region) %>% summarize(lcut = max(lon))
  indx = head(unique(lcuts$lcut), -1) + 1
  
  world = map_data("world")
  world = world[world$long <= 178, ]
  
  zlim = c(-max(abs(field)), max(abs(field)))
  
  ggplot() +
    geom_polygon(data = world, aes(x=long, y=lat, group=group), fill = NA, color="black") +
    geom_hline(yintercept = indy, color = "#F8766D") +
    geom_vline(xintercept = indx, color = "#F8766D") +
    coord_cartesian() +
    theme_void() +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(legend.position="none")
}

# regionalizations
matcombiner = function(mlist, r) {
  c = length(mlist) / r
  do.call(rbind, lapply(0:(c-1), function(x) do.call(cbind, mlist[(1:r) + x*c])))
}

# equal area regionalization
rescale = function(vec, l, u) {
  r = ((vec - min(vec)) / (max(vec) - min(vec))) * (u - l) + l
  as.integer(r)
}
area_cuts = function(arr, regions = 16, rx = sqrt(regions)) {
  mrow = dim(arr)[1]
  mcol = dim(arr)[2]
  mens = dim(arr)[3]
  
  rx = rx
  ry = regions / rx
  
  indx = 1:mrow - (1+mrow)/2
  indx = -acos(indx / (max(indx)))
  indx = indx[seq(1, mrow, length.out = rx+1)]
  indx = rescale(indx, 1, mrow)
  
  indy = as.integer(seq(1, mcol, length.out = ry+1))
  return(list(indx = indx, indy = indy))
}
area_splitter = function(arr, regions = 16, rx = sqrt(regions)) {
  mrow = dim(arr)[1]
  mcol = dim(arr)[2]
  mens = dim(arr)[3]
  
  rx = rx
  ry = regions / rx
  
  indx = 1:mrow - (1+mrow)/2
  indx = -acos(indx / (max(indx)))
  indx = indx[seq(1, mrow, length.out = rx+1)]
  indx = rescale(indx, 1, mrow)
  
  indy = as.integer(seq(1, mcol, length.out = ry+1))
  
  out = lapply(1:rx, function(x) {
    lapply(1:ry, function(y) arr[(indx[x] + as.integer(x>1)):indx[x+1], 
                                 (indy[y] + as.integer(y>1)):indy[y+1],])
  })
  unlist(out, recursive = FALSE)
}

#### Regaionlized significant differences
nc.prior = nc_open('../research/assimilation-cfr/data/tas_prior_da_hydro_r.1000-2000_d.16-Feb-2018.nc')
prior_ind = read.csv("../research/assimilation-cfr/data/prior_ens.txt", header = F)$V1

# 16 region setting
regions = 16
rs = sqrt(regions)

# import all of prior
prior = prep_prior(nc.prior)
prior = prior[,,prior_ind]
prior.split = area_splitter(prior, regions)

region.map = lapply(1:regions, function(x) matrix(x, nrow(prior.split[[x]][,,1]), ncol(prior.split[[x]][,,1])))
region.map = matcombiner(region.map, rs)
region_plot(region.map, nc.prior, main = "Geographic Regions")
ggsave(paste0("../research/assimilation-cfr/paper/misc/", "fancy_regions.png"), width = 5, height = 3.2)

region_plot_simple(region.map, nc.prior)
ggsave(paste0("../research/assimilation-cfr/paper/misc/", "regions.png"), width = 5, height = 3.2)
