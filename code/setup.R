# rm(list = ls())
# gc()

library(extdepth)
library(ncdf4)
library(dplyr)
library(reshape2)
library(plotly)

source('code/prep_functions.R')

##### CONNECT TO DATA #####
# read ensembles and prior ncdf4 objects
nc.post = nc_open('/Users/trevh/Research/assimilation-cfr/data/tas_ens_da_hydro_r.1000-2000_d.16-Feb-2018.nc')
nc.prior = nc_open('/Users/trevh/Research/assimilation-cfr/data/tas_prior_da_hydro_r.1000-2000_d.16-Feb-2018.nc')

# read ensembles and prior ncdf4 objects
# nc.post = nc_open('/Users/Trevor/research/assimilation-cfr/data/tas_ens_da_hydro_r.1000-2000_d.16-Feb-2018.nc')
# nc.prior = nc_open('/Users/Trevor/research/assimilation-cfr/data/tas_prior_da_hydro_r.1000-2000_d.16-Feb-2018.nc')


lat.basis = 40
lon.basis = 26

# basis = create_basis(lat.basis, lon.basis, nc.prior)
# proj = solve(t(basis) %*% basis) %*% t(basis)

# prep prior
prior.sub = read.csv("data/prior_ens.txt", header = F)
prior.sub = as.vector(prior.sub[,1])

prior.ens = prep_prior(nc.prior)
# prior.ens = prior.ens[,,prior.sub]