rm(list = ls())
gc()

library(extdepth)
library(ncdf4)
library(dplyr)
library(reshape2)
library(ggplot2)
library(OpenImageR)

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

# depth CDF
depth = function(g, fmat) {
  
  # Computes the depth values of a function with respect to a set of functions (fmat)
  fn = ncol(fmat)
  depth = rep(0, length(g))
  
  for (row in 1:nrow(fmat)) {
    diff = abs(sum(sign(g[row] - fmat[row,])))
    depth[row] = 1 - (diff / fn)
  }
  
  return(depth)
}
meandepth = function(gmat, fmat) {
  apply(gmat, 2, function(x) mean(depth(x, fmat)))
}
ks.mean = function(f, g) {
  fed = meandepth(f, f)
  ged = meandepth(g, f)
  
  f.surv = rev(c(0, sort(fed)))
  g.cdf = sapply(1:length(f.surv), function(x) mean(ged > f.surv[x]))
  f.cdf = sapply(1:length(f.surv), function(x) mean(fed > f.surv[x]))
  
  ks = max(abs(f.cdf - g.cdf))
  rate = sqrt((ncol(g)*ncol(f)) / (ncol(g) + ncol(f)))
  ks_pval(rate*ks)
}
ks_pval = function(t, n = 20) {
  2*(sum(sapply(1:n, function(x) (-1)^(x-1) * exp(-2*(x^2)*t^2))))
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
  dimnames(field) = list(lats, lons-180)
  
  field.gg = melt(field)
  colnames(field.gg) = c("lat", "lon", "value")
  
  world = map_data("world")
  world = world[world$long <= 178, ]
  
  zlim = c(-max(abs(field)), max(abs(field)))
  
  ggplot() +
    geom_raster(data = field.gg, aes(x=lon, y=lat, fill=value), interpolate = TRUE) +
    # geom_contour(data = field.gg, aes(x=lon, y=lat, z=value)) +
    geom_polygon(data = world, aes(x=long, y=lat, group=group), fill = NA, color="black") +
    coord_cartesian() +
    # scale_fill_distiller(palette = "Spectral", limits = zlim) +
    scale_fill_gradient2(midpoint=0, low="blue", mid="white", high="red", guide = FALSE) +
    theme_void() +
    ggtitle(main) +
    theme(plot.title = element_text(hjust = 0.5))
}

# function for creating file name with leading zeros
# makes it easier to process them sequentially
rename <- function(x, txt=''){
  if (x < 10) {
    return(name <- paste('000',x, txt, 'plot.png',sep=''))
  }
  if (x < 100 && x >= 10) {
    return(name <- paste('00',x, txt, 'plot.png', sep=''))
  }
  if (x >= 100) {
    return(name <- paste('0', x, txt, 'plot.png', sep=''))
  }
}

# load p vals and data
pvals = read.csv("../Climate/adjusted_pvals")$x

nc.post = nc_open('../data/tas_ens_da_hydro_r.1000-2000_d.16-Feb-2018.nc')
nc.prior = nc_open('../data/tas_prior_da_hydro_r.1000-2000_d.16-Feb-2018.nc')

# no downsampling
downsamp = 4

prior = prep_prior(nc.prior)
prior = sapply(1:dim(prior)[3], function(x) down_sample_image(prior[,,x], downsamp))

prior.depths = meandepth(prior, prior)
prior.ranks = rank(prior.depths) / length(prior.depths)
cr = central_region(prior, prior.ranks, 0.05)

#### AVERAGE SIGNIFICANT GIF
times = (201:998)
# times = as.integer(seq(1, 998, length.out = 100))

for(t in times) {
  
  post = prep_post(nc.post, t)
  post = sapply(1:dim(post)[3], function(x) down_sample_image(post[,,x], downsamp))
  post.depths = meandepth(post, prior)
  
  field_plot(remove_cr(cr, post, downsamp),
             nc.prior,
             downsamp = downsamp,
             main = paste0("Year ", t))
  
  ggsave(paste0("animation/", rename(t)), width = 5, height = 3.2)
  
}

# creates the gif
home = getwd()
setwd("animation/")
my_command <- 'convert *.png -delay 10 -loop 0 climate.gif'
system(my_command)
setwd(home)

##### Separate section
library(fda)
library(extdepth)

cw = CanadianWeather

temps = unname(as.matrix(cw$dailyAv[,,1]))
temps = cbind(temps, temps + 5, temps - 5)

plot_cr(temps, alpha = 0.95, main = "5% CR")
ggsave(paste0("animation2/", rename(1, "_cr_")), width = 5, height = 3.2)

plot_cr(temps, alpha = 0.75, main = "25% CR")
ggsave(paste0("animation2/", rename(2, "_cr_")), width = 5, height = 3.2)

plot_cr(temps, alpha = 0.50, main = "50% CR")
ggsave(paste0("animation2/", rename(3, "_cr_")), width = 5, height = 3.2)

plot_cr(temps, alpha = 0.25, main = "75% CR")
ggsave(paste0("animation2/", rename(4, "_cr_")), width = 5, height = 3.2)

plot_cr(temps, alpha = 0.05, main = "95% CR")
ggsave(paste0("animation2/", rename(5, "_cr_")), width = 5, height = 3.2)

# home = getwd()
# setwd("animation2/")
# my_command <- 'convert *.png -delay 1000 cr.gif'
# system(my_command)
# setwd(home)

