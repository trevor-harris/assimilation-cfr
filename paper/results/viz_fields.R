rm(list = ls())
gc()

library(extdepth)
library(ncdf4)
library(dplyr)
library(reshape2)
library(ggplot2)
library(OpenImageR)

source("research/assimilation-cfr/code/depth_tests.R")
source("research/assimilation-cfr/code/depths.R")
source("research/assimilation-cfr/code/simulation.R")


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
field_plot <- function(field, nc, main = "", downsamp = 1, zlim = c(-max(abs(field)), max(abs(field)))) {
  
  lats = as.vector(nc$dim$lat$vals)[seq(1, 96, by=downsamp)]
  lons = as.vector(nc$dim$lon$vals)[seq(1, 144, by=downsamp)]
  dimnames(field) = list(lats, ifelse(lons >= 180, lons - 360, lons))
  
  field.gg = melt(field)
  colnames(field.gg) = c("lat", "lon", "value")
  
  world = map_data("world")
  world = world[world$long <= 178, ]
  
  zlim = c(-max(abs(field)), max(abs(field)))
  
  ggplot() +
    geom_raster(data = field.gg, aes(x=lon, y=lat, fill=value), interpolate = TRUE) +
    geom_polygon(data = world, aes(x=long, y=lat, group=group), fill = NA, color="black") +
    coord_cartesian() +
    scale_fill_gradient2(midpoint=0, low="blue", mid="white", high="red") +
    # scale_fill_gradient2(midpoint=0, low="#4393c3", mid="white", high="#d6604d") +
    # scale_fill_gradient2(midpoint=0, low="#2166ac", mid="white", high="#b2182b") +
    theme_void() +
    ggtitle(main) +
    theme(plot.title = element_text(hjust = 0.5))
}


###### Read in the temp data results (for K values)
read_era = function(dir, file) {
  cbind(readRDS(paste0(dir, file)), era = as.numeric(strsplit(file, "\\D+")[[1]][-1]))
}

# import raw size data
dir = "../temp/power/independent/"
dir = "research/assimilation-cfr/paper/results/results/"
files = list.files(dir)

temperature = read_era(dir, files[1])
for(f in 3:length(files)) {
  temperature = rbind(temperature, read_era(dir, files[f]))
}
temperature = rbind(temperature, read_era(dir, files[2]))
temperature$stat = temperature$stat / (sqrt((100*100)/(100 + 100)))
temperature[["time"]] = 1:998


##### Actual data
nc.post = nc_open('research/assimilation-cfr/data/tas_ens_da_hydro_r.1000-2000_d.16-Feb-2018.nc')
nc.prior = nc_open('research/assimilation-cfr/data/tas_prior_da_hydro_r.1000-2000_d.16-Feb-2018.nc')

prior_ind = read.csv("research/assimilation-cfr/data/prior_ens.txt", header = F)$V1

prior = prep_prior(nc.prior)
prior = flatten(prior[,,prior_ind])

prior.depths = xdepth(prior, prior)


##### Exceedence plots
save_dir = "research/assimilation-cfr/paper/results/"

# CDF of prior to get central regions
prior.ranks = rank(prior.depths) / length(prior.depths)
cr = central_region(prior, prior.ranks, 0.05)

# times = as.integer(seq(1, 998, length.out = 5))
times = c(2, 300, 600, 998)
for(t in times) {
  post = flatten(prep_post(nc.post, t))
  
  k.t = formatC(temperature$stat[t], digits = 2)
  print(field_plot(remove_cr(cr, post),
                   nc.prior,
                   main = paste0("Year ", t, " (K = ", k.t, ")")))
  ggsave(paste0(save_dir, "year", t, "field.png"), width = 5, height = 3.2)
}


##### ACTUAL CLIMATE ##### 
field_plot2 <- function(field, nc, main = "", downsamp = 1, zlim = c(-max(abs(field)), max(abs(field)))) {
  
  lats = as.vector(nc$dim$lat$vals)[seq(1, 96, by=downsamp)]
  lons = as.vector(nc$dim$lon$vals)[seq(1, 144, by=downsamp)]
  dimnames(field) = list(lats, ifelse(lons >= 180, lons - 360, lons))
  
  field.gg = melt(field)
  colnames(field.gg) = c("lat", "lon", "value")
  
  world = map_data("world")
  world = world[world$long <= 178, ]
  
  zlim = c(-max(abs(field)), max(abs(field)))
  
  ggplot() +
    geom_raster(data = field.gg, aes(x=lon, y=lat, fill=value), interpolate = TRUE) +
    geom_polygon(data = world, aes(x=long, y=lat, group=group), fill = NA, color="black") +
    coord_cartesian() +
    # scale_fill_gradient2(midpoint=0, low="#4393c3", mid="white", high="#d6604d") +
    # scale_fill_gradient2(midpoint=0, low="#2166ac", mid="white", high="#b2182b") +
    scale_fill_distiller(palette = "RdBu") +
    theme_void() +
    ggtitle(main) +
    theme(legend.position="none") + 
    theme(plot.title = element_text(hjust = 0.5))
}


field_plot2(matrix(prior[,1], 96, 144), nc.post, main = "Prior ensemble member", downsamp = 1)
ggsave(paste0(save_dir, "prior.png"), width = 5, height = 3.2)


post = prep_post(nc.post, times[4])
field_plot2(post[,,10], nc.post, main = "Posterior ensemble member (Year 998)", downsamp = 1)
ggsave(paste0(save_dir, "posterior.png"), width = 5, height = 3.2)

