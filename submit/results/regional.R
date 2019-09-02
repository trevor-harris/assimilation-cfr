rm(list = ls())
gc()

# library(extdepth)
library(ncdf4)
library(dplyr)
library(reshape2)
library(ggplot2)
library(RColorBrewer)
library(tictoc)
library(future)
library(future.apply)
library(fdasrvf)
library(tictoc)


# set to the top level folder
setwd("/Users/trevorh2/research/assimilation-cfr/submit/")


source("method/depth_tests.R")
source("method/depths.R")
source("method/simulation.R")


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

# reg_ind = numbers used to identify regions in the mask file
reg_ind = c(c(10, 20, 30, 40, 50), c(1000, 2000, 3000, 4000, 5000, 6000, 7000))
reg_names = c(c("Arctic Ocean", "Indian Ocean", "Pacific Ocean", "Atlantic Ocean", "Southern Ocean"), 
              c("Antarctica", "South America", "North America", "Africa", "Europe", "Asia", "Australia"))

# data dimensions
years = 851:1848 
times = 998
lats = 96
lons = 144
ens = 100
reg = length(reg_ind)
kfield = matrix(0, reg, times)
pfield = matrix(0, reg, times)


###### TESTING


# open connection to NCDF4 data files
nc.post = nc_open('data/tas_ens_da_hydro_r.1000-2000_d.16-Feb-2018.nc')
nc.prior = nc_open('data/tas_prior_da_hydro_r.1000-2000_d.16-Feb-2018.nc')

# this file contains a list of ensemble members used to subset the background (prior)
prior_ind = read.csv("data/prior_ens.txt", header = F)$V1

# Mask file identifies regions
mask = read.csv("data/mask.csv", stringsAsFactors = F)
mask = as.matrix(mask)
mask = apply(mask, 2, rev)

# import all of prior
prior = prep_prior(nc.prior)
prior = prior[,,prior_ind]

# convert to list of matricies. Each list element is the entire reconstruction for a given region.
prior_reg = vector("list", length(reg_ind))
for(r in 1:length(reg_ind)) {
  prior_reg[[r]] =  matrix(prior[mask == reg_ind[r]], ncol = 100)
}

# import entire posterior (takes a long time)
post = vapply(1:times, function(t) prep_post(nc.post, t), array(0, dim=c(lats, lons, ens)))

# Convert to region list
post_reg = vector("list", length(reg_ind))
for(r in 1:length(reg_ind)) {
  reg.r = post[mask == reg_ind[r]]
  post_reg[[r]] =  array(post[mask == reg_ind[r]], dim = c(length(reg.r)/(ens*times), ens, times))
}

# perform tests
plan(multiprocess)
options(future.globals.maxSize = 4000*1024^2)

for(r in 1:length(reg_ind)) {
  tic(reg_names[r])
  
  post_reg.r = post_reg[r][[1]]
  prior_reg.r = prior_reg[r][[1]]
  
  ktest = future_sapply(1:times, function(t) kolm(prior_reg.r, post_reg.r[,,t]))
  kfield[r,] = ktest[1,]
  pfield[r,] = ktest[2,]
  
  toc()
}

###### PLOTTING

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


# plot KD values
ggplot(kseries, aes(Year, value)) +
  geom_point(size = 0.2) +
  geom_smooth(color = "red") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.position="none") +
  scale_x_continuous(breaks = c(1000, 1400, 1800)) +
  scale_y_continuous(breaks = c(0.1, 0.5, 1)) +
  scale_color_identity(guide = 'legend') +
  xlab("Year") +
  ylab("KD") +
  facet_wrap(vars(Region), 3, 4, scales = "free_x")
# ggsave("results/kd_region.png", width = 6, height = 3.8)


# plot p-values
nonsigpseries = pseries[pseries$value > 0.05,]

ggplot() +
  geom_point(data = pseries, aes(Year, value), size = 0.2, color = "black") +
  geom_point(data = nonsigpseries, aes(Year, value), size = 0.2, color = "pink") +
  theme_classic() + 
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.position="none") +
  scale_x_continuous(breaks = c(1000, 1400, 1800)) +
  scale_y_continuous(breaks = c(0, 0.5, 1)) +
  scale_color_identity(guide = 'legend') +
  xlab("Year") +
  ylab("p-value") +
  facet_wrap(vars(Region), 3, 4, scales = "free_x")
# ggsave("results/pval_region.png", width = 6, height = 3.8)
