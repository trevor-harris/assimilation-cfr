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
    ggtitle(main) +
    theme(plot.title = element_text(hjust = 0.5))
}
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
    theme(plot.title = element_text(hjust = 0.5)) +
    guide_legend(title="Temperature")
}


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

# import all of post
# post = vapply(1:998, function(t) flatten(prep_post(nc.post, t)), matrix(0, 13824, 100))

#### average warming/cooling over time
post.gain = vapply(1:998, function(t) remove_cr(cr, post[,,t]), matrix(0, 96, 144))

# significant trends
sig.trend = t(flatten(post.gain))

dt.gg = melt(sig.trend)
dt.gg[["Year"]] = dt.gg[["Var1"]] + 850

ggplot(dt.gg, aes(Year, value, group = Var2)) +
  geom_line(alpha = 0.2, size = 0.2) +
  geom_hline(yintercept = 0, color = "red") +
  theme_classic() + 
  xlab("Year") +
  ylab("Temperature in C")
ggsave(paste0(save_dir,"trend.png"), width = 5, height = 3.2)


# Average map
post.gain.avg = apply(post.gain, c(1, 2), mean)

field_plot(post.gain.avg, nc.post, main = "Average exceedence (all years)")
ggsave(paste0(save_dir, "average.png"), width = 5, height = 3.2)




#### Example times
times = c(2, 300, 600, 998)
for(t in times) {
  post.t = post[,,t]
  
  k.t = formatC(temperature$stat[t], digits = 2)
  print(field_plot(remove_cr(cr, post.t),
                   nc.prior,
                   main = paste0("Year ", t, " exceedence ","(K = ", k.t, ")")))
  ggsave(paste0(save_dir, "year", t, "field.png"), width = 5, height = 3.2)
}

# trends
trends = vapply(1:(96*144), function(x) colMeans(post[x,,]), rep(0, 998))
f.trend = vapply(1:(96*144), function(x) smooth.spline(trends[,x], nknots = 25)$y, rep(0, 998))

prior.mean = rowMeans(prior)
d.trend = vapply(1:(96*144), function(x) f.trend[,x] - prior.mean[x], rep(0, 998))

dt.gg = melt(d.trend)
ggplot(dt.gg, aes(Var1, value, group = Var2)) +
  geom_line(alpha = 0.2) +
  theme_classic() + 
  xlab("Time") +
  ylab("Temperature")



##### ACTUAL CLIMATE ##### 
field_plot2(matrix(prior[,1], 96, 144), nc.post, main = "Prior ensemble member", downsamp = 1)
# ggsave(paste0(save_dir, "prior.png"), width = 5, height = 3.2)


post = prep_post(nc.post, times[4])
field_plot2(post[,,10], nc.post, main = "Posterior ensemble member (Year 998)", downsamp = 1)
# ggsave(paste0(save_dir, "posterior.png"), width = 5, height = 3.2)

