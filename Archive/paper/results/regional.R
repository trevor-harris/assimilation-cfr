rm(list = ls())
gc()

years = 851:1848
ens = 100

library(extdepth)
library(ncdf4)
library(dplyr)
library(reshape2)
library(ggplot2)
library(RColorBrewer)
library(OpenImageR)
library(tictoc)
library(future)
library(future.apply)
library(fdasrvf)

setwd("C:/Users/trevorh2/research/")
source("assimilation-cfr/code/depth_tests.R")
source("assimilation-cfr/code/depths.R")
source("assimilation-cfr/code/simulation.R")

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
region_plot <- function(field, nc, reg_names, main = "", downsamp = 1, zlim = c(-max(abs(field)), max(abs(field)))) {
  
  lats = as.vector(nc$dim$lat$vals)[seq(1, 96, by=downsamp)]
  lons = as.vector(nc$dim$lon$vals)[seq(1, 144, by=downsamp)]
  dimnames(field) = list(lats, ifelse(lons >= 180, lons - 360, lons))
  
  field.gg = melt(field)
  colnames(field.gg) = c("lat", "lon", "Region")
  field.gg$Region = as.factor(field.gg$Region)
  levels(field.gg$Region) = reg_names
  
  world = map_data("world")
  world = world[world$long <= 178, ]
  
  zlim = c(-max(abs(field)), max(abs(field)))
  
  ggplot() +
    geom_raster(data = field.gg, aes(x=lon, y=lat, fill=Region), interpolate = TRUE) +
    geom_polygon(data = world, aes(x=long, y=lat, group=group), fill = NA, color="black") +
    coord_cartesian() +
    scale_fill_hue(l=c(seq(55, 80,length.out = 10))) +
    theme_void() +
    ggtitle(main) +
    theme(legend.position="none") +
    theme(plot.title = element_text(hjust = 0.5))
}
smooth.mat = function(m) {
  t(smooth.data(t(smooth.data(m, 3)), 3))
}

save_dir = "../research/assimilation-cfr/paper/results/"

#### Regionlized significant differences
nc.post = nc_open('../research/assimilation-cfr/data/tas_ens_da_hydro_r.1000-2000_d.16-Feb-2018.nc')
nc.prior = nc_open('../research/assimilation-cfr/data/tas_prior_da_hydro_r.1000-2000_d.16-Feb-2018.nc')
prior_ind = read.csv("../research/assimilation-cfr/data/prior_ens.txt", header = F)$V1

# read in the mask file
mask = read.csv("../research/assimilation-cfr/data/mask.csv", stringsAsFactors = F)
mask = as.matrix(mask)
mask = apply(mask, 2, rev)

reg_ind = c(c(10, 20, 30, 40, 50), c(1000, 2000, 3000, 4000, 5000, 6000, 7000))
reg_names = c(c("Arctic Ocean", "Indian Ocean", "Pacific Ocean", "Atlantic Ocean", "Southern Ocean"), 
              c("Antarctica", "South America", "North America", "Africa", "Europe", "Asia", "Australia"))


times = 998
lats = 96
lons = 144
ens = 100
reg = length(reg_ind)
kfield = matrix(0, reg, times)
pfield = matrix(0, reg, times)


# import all of prior
prior = prep_prior(nc.prior)
prior = prior[,,prior_ind]

# smooth out prior
prior = vapply(1:ens, function(t) smooth.mat(prior[,,t]), matrix(0, lats, lons))

prior_reg = vector("list", length(reg_ind)) 
for(r in 1:length(reg_ind)) {
  prior_reg[[r]] =  matrix(prior[mask == reg_ind[r]], ncol = 100)
}

# import posterior at time t
post = vapply(1:times, function(t) prep_post(nc.post, t), array(0, dim=c(lats, lons, ens)))

# smooth out posterior
for(t in 1:times) {
  for(e in 1:ens) {
    post[,,e,t] = smooth.mat(post[,,e,t])
  }
}

post_reg = vector("list", length(reg_ind))
for(r in 1:length(reg_ind)) {
  reg.r = post[mask == reg_ind[r]]
  post_reg[[r]] =  array(post[mask == reg_ind[r]], dim = c(length(reg.r)/(ens*times), ens, times))
}


# perform tests
plan(multiprocess)
options(future.globals.maxSize = 4000*1024^2)

for(r in 1:length(reg_ind)) {
  post_reg.r = post_reg[r][[1]]
  prior_reg.r = prior_reg[r][[1]]
  
  ktest = future_sapply(1:times, function(t) kolm(prior_reg.r, post_reg.r[,,t]))
  kfield[r,] = ktest[1,]
  pfield[r,] = ktest[2,]
}


# reassemble into series
kseries = melt(t(kfield)) %>% 
  mutate(Region = as.factor(Var2),
         Year = rep(years, 12),
         value = value)
levels(kseries$Region) = reg_names

pseries = melt(t(pfield)) %>%
  mutate(Region = as.factor(Var2),
         Year = rep(years, 12),
         value = p.adjust(value, method = "BY"))
levels(pseries$Region) = reg_names

# plot K
ggplot(kseries, aes(Year, value, color = Region)) +
  geom_point(size = 0.2) +
  geom_smooth(color = "black") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.position="none") +
  scale_x_continuous(breaks = c(1000, 1400, 1800)) +
  scale_y_continuous(breaks = c(0.1, 0.5, 1)) +
  # scale_color_manual(values = colorRampPalette(brewer.pal(11,"Paired"))(12)) +
  xlab("Year") +
  ylab("K") +
  facet_wrap(vars(Region), 3, 4)
ggsave(paste0(save_dir, "k_region.png"), width = 6, height = 3.5)

# plot Pval of K
ggplot(pseries, aes(Year, value, color = Region)) +
  geom_point(size = 0.2) +
  theme_classic() + 
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.position="none") +
  scale_x_continuous(breaks = c(1000, 1400, 1800)) +
  scale_y_continuous(breaks = c(0, 0.5, 1)) +
  xlab("Year") +
  ylab("p-value") +
  facet_wrap(vars(Region), 3, 4)
ggsave(paste0(save_dir, "pval_region.png"), width = 6, height = 3.5)


# reassemble into maps
# kmaps = vapply(1:times, function(t) {

#   matcombiner(rmap(kfield[,t], lats/4, lons/4), 4)
# }, matrix(0, lats, lons))
# 
# pmaps = vapply(1:times, function(t) {
#   matcombiner(rmap(pfield[,t], lats/4, lons/4), 4)
# }, matrix(0, lats, lons))
# 
# field_plot(kmaps[,,3], nc.post)
# field_plot(pmaps[,,3], nc.post)

