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
library(gridExtra)
library(cowplot)

devtools::install_github('trevor-harris/kstat')
library(kstat)

# code for simulating guassian processes, t-processes, and plotting functions
source("util/simulation.R")

# code for importing and processing the ensemble data
source("util/import_data.R")


# helper function for multiplot
get_legend<-function(myggplot){
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
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

# convert to anomalies with respect to the background
prior = prior.all[,,70] - prior.mean
post = post.all[,,70] - prior.mean

# pull out coordinates for plotting
lats = as.vector(nc.prior$dim$lat$vals)
lons = as.vector(nc.prior$dim$lon$vals)

dimnames(prior) = list(lats, ifelse(lons >= 180, lons - 360, lons))
dimnames(post) = list(lats, ifelse(lons >= 180, lons - 360, lons))

# convert to ggplot2 compliant form
prior = melt(prior)
post = melt(post)

colnames(prior) = c("lat", "lon", "val")
colnames(post) = c("lat", "lon", "val")

prior[["dist"]] = "Background ensemble member 70"
post[["dist"]] = "Analysis ensemble member 70 (850 CE)"

# import world map overlay
world = map_data("world")
world = world[world$long <= 178, ]

### add proxy locations
prox = readMat("maps/proxydata_aprmar_lmr_v0.2.0_pages2k_v2.mat")
smap = prox$lmr2k.data[850,]
smap = 1-sapply(smap, is.nan)
smap = which(smap == 1)

# extract lat and lon and reparameterize lon to match ggplot
lat = prox$p.lat[smap]
lon = prox$p.lon[smap]
lon = ifelse(lon >= 180, lon - 360, lon)
prox_loc = data.frame(lat = lat, lon = lon)

# 1. Create the plots
# background plot (lp)
lp <- ggplot() +
  geom_raster(data = prior, aes(x=lon, y=lat, fill=val)) +
  geom_polygon(data = world, aes(x=long, y=lat, group=group),
               fill = NA, color="grey30") +
  coord_cartesian() +
  theme_void() +
  scale_fill_distiller(palette = "RdBu",
                       limits = c(-2.5, 4),
                       name = "Temp.\nAnomaly (C)") +
  ggtitle("Background ensemble member 70") +
  theme(plot.title = element_text(hjust = 0.5),
        text = element_text(size = 14))


# analysis plot (rp)
rp <- ggplot() +
  geom_raster(data = post, aes(x=lon, y=lat, fill=val)) +
  geom_polygon(data = world, aes(x=long, y=lat, group=group), 
               fill = NA, color="grey30") +
  geom_point(data = prox_loc, aes(x=lon, y=lat), 
             color = "black", fill = 'red', shape = 25, size = 3) +
  coord_cartesian() +
  theme_void() +
  scale_fill_distiller(palette = "RdBu",
                       limits = c(-2.5, 4),
                       name = "Temp.\nAnomaly (C)") +
  ggtitle("Analysis ensemble member 70 (850 CE)") +
  theme(plot.title = element_text(hjust = 0.5),
        text = element_text(size = 14)) +
  theme(legend.position = "none")

# Save the legend
legend <- get_legend(lp)
lp <- lp + theme(legend.position="none")

# Arrange plot
ba_anom = grid.arrange(lp, rp, legend, ncol=3, widths=c(2.4, 2.4, 0.8))

save_plot("maps/background_analysis_anomaly.png", ba_anom, base_aspect_ratio = 2.9)
