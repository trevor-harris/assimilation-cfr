rm(list = ls())
gc()

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

# prep prior
prior.sub = read.csv("data/prior_ens.txt", header = F)
prior.sub = as.vector(prior.sub[,1])

prior.ens = prep_prior(nc.prior)
prior.ens = prior.ens[,,prior.sub]

prior = prior.ens
posterior = prep_post_ens(nc.post, 20)


lat = 96
lon = 144
betas = matrix(0, lat, lon)
for(i in 1:lat) {
  for(j in 1:lon) {
    betas[i,j] = mean(posterior[i,j,]-prior[i,j,])
    # betas[i,j] = t.test(prior[i,j,], posterior[i,j,], paired = T)$statistic
  }
}

perms = 1000
permuted_betas = array(0, c(lat, lon, perms))

for(p in 1:perms) {
  permuted = permute_fields(prior, posterior)
  prior.p = permuted[[1]]
  posterior.p = permuted[[2]]
  
  betas.p = matrix(0, lat, lon)
  for(i in 1:lat) {
    for(j in 1:lon) {
      betas.p[i,j] = mean(posterior.p[i,j,]-prior.p[i,j,])
      # betas[i,j] = t.test(prior.p[i,j,], posterior.p[i,j,], paired = T)$statistic
    }
  }
  
  permuted_betas[,,p] = betas.p
}

# observed betas
beta_obs[,,s] = betas

# ED based upper and lower
betas_flat = flatten(permuted_betas)
betas_ed = edepth_set(betas_flat)
betas_cr = central_region(betas_flat, betas_ed)

ed_lower[,,s] = matrix(betas_cr[[1]], 30, 30)
ed_upper[,,s] = matrix(betas_cr[[2]], 30, 30)