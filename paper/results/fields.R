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


source("../research/assimilation-cfr/code/depth_tests.R")
source("../research/assimilation-cfr/code/depths.R")
source("../research/assimilation-cfr/code/simulation.R")

plan(multiprocess)


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
  
  field = prior[,,1]
  nc = nc.prior
  downsamp = 1
  
  lats = as.vector(nc$dim$lat$vals)[seq(1, 96, by=downsamp)]
  lons = as.vector(nc$dim$lon$vals)[seq(1, 144, by=downsamp)]
  dimnames(field) = list(lats, ifelse(lons >= 180, lons - 360, lons))
  
  field.gg = melt(field)
  colnames(field.gg) = c("lat", "lon", "Temp")
  
  world = map_data("world")
  world = world[world$long <= 178, ]
  
  zlim = c(-max(abs(field)), max(abs(field)))
  
  ggplot() +
    geom_raster(data = field.gg, aes(x=lon, y=lat, fill=Temp), interpolate = TRUE) +
    geom_polygon(data = world, aes(x=long, y=lat, group=group), fill = NA, color="black") +
    coord_cartesian() +
    scale_fill_gradient2(midpoint=0, low="blue", mid="white", high="red") +
    # scale_fill_gradient2(midpoint=0, low="#4393c3", mid="white", high="#d6604d") +
    # scale_fill_gradient2(midpoint=0, low="#2166ac", mid="white", high="#b2182b") +
    theme_void() +
    # ggtitle(main) +
    theme(plot.title = element_text(hjust = 0.5))
}


##### Actual data
nc.post = nc_open('../research/assimilation-cfr/data/tas_ens_da_hydro_r.1000-2000_d.16-Feb-2018.nc')
nc.prior = nc_open('../research/assimilation-cfr/data/tas_prior_da_hydro_r.1000-2000_d.16-Feb-2018.nc')

prior_ind = read.csv("../research/assimilation-cfr/data/prior_ens.txt", header = F)$V1


prior = prep_prior(nc.prior)
prior = flatten(prior[,,prior_ind])
prior.depths = xdepth(prior, prior)

# also pull in K values
dir = "../research/assimilation-cfr/paper/results/results/"
files = list.files(dir)
read_era = function(dir, file) {
  cbind(readRDS(paste0(dir, file)), era = as.numeric(strsplit(file, "\\D+")[[1]][-1]))
}

temperature = read_era(dir, files[1])
for(f in 3:length(files)) {
  temperature = rbind(temperature, read_era(dir, files[f]))
}
temperature = rbind(temperature, read_era(dir, files[2]))
temperature[["time"]] = years


##### Exceedence plots
save_dir = "research/assimilation-cfr/paper/results/"

# CDF of prior to get central regions
prior.ranks = rank(prior.depths) / length(prior.depths)
cr = central_region(prior, prior.ranks, 0.05)

# Import post years of interest and remove the cr
post = data.frame()


times = c(2, 300, 600, 998)
kt = numeric(0)
for(t in times) {
  post.t = melt(remove_cr(cr, flatten(prep_post(nc.post, t))))
  post.t[["time"]] = t
  names(post.t) = c("lat", "lon", "val", "time")
  
  post = rbind(post, post.t)
  
  kt = c(kt, formatC(temperature$stat[t] /  sqrt(100*100 / 200), digits = 2))
}

# format into lats and lons properly
nc = nc.prior
lats = as.vector(nc$dim$lat$vals)
lons = as.vector(nc$dim$lon$vals)

lats = rep(lats, 144*length(times))
lons = rep(rep(lons, each = 96), 4)
lons = ifelse(lons >= 180, lons - 360, lons)

post[["lat"]] = lats
post[["lon"]] = lons

post[["time"]] = factor(post[["time"]], levels = c("2", "300", "600", "998"), 
                        labels = c(paste0(years[2], "CE ","(K = ", kt[1], ")"),
                                   paste0(years[300], "CE ","(K = ", kt[2], ")"),
                                   paste0(years[600], "CE ","(K = ", kt[3], ")"),
                                   paste0(years[998], "CE ","(K = ", kt[4], ")")))

# get world "underlay"
world = map_data("world")
world = world[world$long <= 178, ]

ggplot() +
  geom_raster(data = post, aes(x=lon, y=lat, fill=val), interpolate = TRUE) +
  geom_polygon(data = world, aes(x=long, y=lat, group=group), fill = NA, color="black") +
  coord_cartesian() +
  scale_fill_gradient2(midpoint=0, low="blue", mid="white", high="red") +
  theme_void() +
  facet_wrap(. ~ time, nrow = 2) +
  theme(plot.title = element_text(hjust = 0.5))
ggsave("../research/assimilation-cfr/paper/results/multiyear_fields.png")
