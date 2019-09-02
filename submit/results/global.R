rm(list = ls())
gc()

library(ncdf4)
library(tictoc)
library(ggplot2)
library(dplyr)
library(reshape2)
library(latex2exp)

setwd("../research/assimilation-cfr/paper/submit")

source("code/depth_tests.R")
source("code/depths.R")
source("code/simulation.R")


###### TESTING


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

nc.post = nc_open('data/tas_ens_da_hydro_r.1000-2000_d.16-Feb-2018.nc')
nc.prior = nc_open('data/tas_prior_da_hydro_r.1000-2000_d.16-Feb-2018.nc')

prior_ind = read.csv("data/prior_ens.txt", header = F)$V1

# import all of prior
prior = prep_prior(nc.prior)
prior = flatten(prior[,,prior_ind])

times = 1:998
years = 851:1848
k.t = matrix(0, length(times), 3)

# run test for each year in reconstruction
for(t in 1:length(times)) {
  post.t = flatten(prep_post(nc.post, times[t]))
  k.t[t,] = kolm(prior, post.t)
}


###### PLOTTING

kstats = data.frame(stat = k.t[,1], pval = k.t[,2], Year = years)

ggplot(kstats, aes(x=Year, y=stat)) +
  geom_point(size = 2) +
  geom_smooth(color = "red", size = 2, se = F) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5), 
        text = element_text(size=24)) +
  ylab("KD")
ggsave("results/kd_global.png", width = 8, heigh = 6)


ggplot(kstats, aes(x=Year, y=pval)) +
  geom_point(size = 2) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5), 
        text = element_text(size=24)) +
  ylab("p-values")
ggsave("results/pval_global.png", width = 8, heigh = 6)