rm(list = ls())
gc()

years = 851:1848
ens = 100

library(extdepth)
library(ncdf4)
library(dplyr)
library(reshape2)
library(ggplot2)
library(plotly)

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
    geom_raster(data = field.gg, aes(x=lon, y=lat, fill=value), interpolate = F) +
    geom_polygon(data = world, aes(x=long, y=lat, group=group), fill = NA, color="black") +
    coord_cartesian() +
    scale_fill_distiller(palette = "RdBu") +
    theme_void() +
    ggtitle(main) +
    theme(legend.position="none") + 
    theme(plot.title = element_text(hjust = 0.5))
}


##### Actual data
nc.post = nc_open('../research/assimilation-cfr/data/tas_ens_da_hydro_r.1000-2000_d.16-Feb-2018.nc')
nc.prior = nc_open('../research/assimilation-cfr/data/tas_prior_da_hydro_r.1000-2000_d.16-Feb-2018.nc')

prior_ind = read.csv("../research/assimilation-cfr/data/prior_ens.txt", header = F)$V1

prior = prep_prior(nc.prior)
post = prep_post(nc.post, t = 900)

plot_ly(showscale = FALSE) %>% 
  add_surface(z = ~prior[,,1]) %>%
  add_surface(z = ~prior[,,10] + 50) %>%
  add_surface(z = ~prior[,,20] + 100)


##### ACTUAL CLIMATE ##### 
field_plot2(prior[,,10], nc.post, main = "Model Output", downsamp = 1)
ggsave("../research/assimilation-cfr/Presentation/model_output.png")


field_plot2(post[,,10], nc.post, main = "Data Assimilation Output", downsamp = 1)
ggsave("../research/assimilation-cfr/Presentation/data_assim.png")

field_plot2(prior[,,10], nc.post, main = "")

plot_ly(showscale = F) %>% add_surface(z = prior[,,1])






prep_post2 = function(nc.post, t) {
  
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
  
  ens = vapply(1:n.ens, function(x) ens[,,x], FUN.VALUE = matrix(0, nrow = n.lat, ncol = n.lon))
  
  return(ens)
}

post2 = prep_post2(nc.post, 900)
xax <- list(title = "Longitude")
yax <- list(title = "Latitutde")
zax <- list(title = "Temperature")

plot_ly(showscale = F) %>% 
  add_surface(z = post2[,,1]) %>%
  layout(scene = list(
      xaxis = list(title = "Longitude"),
      yaxis = list(title = "Latitude"),
      zaxis = list(title = "Temperature")
    ))

field_plot2(post2[,,1], nc.post, main = "", downsamp = 1)
ggsave("../research/assimilation-cfr/Presentation/post.png")
