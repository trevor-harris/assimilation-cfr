rm(list = ls())
gc()




########### READ ME #############

# you must change the working directory to be the submitted_code folder
# none of this will work otherwise
# mine is left here as an example

########## Example
# setwd("/Users/trevh/research/assimilation-cfr/submitted_code/")

#################################





library(ncdf4)
library(tictoc)
library(ggplot2)
library(dplyr)
library(reshape2)

devtools::install_github('trevor-harris/kstat')
library(kstat)

# code for simulating guassian processes, t-processes, and plotting functions
source("util/simulation.R")

# code for importing and processing the ensemble data
source("util/import_data.R")


# plots the regions outlined in the paper
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


#### Regionlized significant differences
# prior (background)
nc.prior = nc_open('data/tas_prior_da_hydro_r.1000-2000_d.16-Feb-2018.nc')
prior_ind = read.csv("data/prior_ens.txt", header = F)$V1

# post (analysis)
nc.post = nc_open('data/tas_ens_da_hydro_r.1000-2000_d.16-Feb-2018.nc')

# read in the mask file. Mask 2 include a "region" which denotes the borders of other regions
mask = read.csv("data/mask2.csv", stringsAsFactors = F)[,2:145]
mask = as.matrix(mask)
mask = apply(mask, 2, rev)

reg_names = c(c("Arctic Ocean", "Indian Ocean", "Pacific Ocean", "Atlantic Ocean", "Southern Ocean"), 
              c("Antarctica", "South America", "North America", "Africa", "Europe", "Asia", "Australia", "Borders"))

cc = c("#4bb1bf", "#235354", "#5c9c9a", "#8be0e0", "#102829",
       "#c4c2b3", "#84553b", "#cea672", "#834210", "#af8c78", "#A98E44", "#3e1d0d",
       "#000000")

region_plot(mask, nc.prior, reg_names, cc)
# ggsave("maps/regions.png", width = 9, height = 6)
