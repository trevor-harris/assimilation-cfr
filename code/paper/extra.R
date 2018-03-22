source('code/setup.R')

##### CONNECT TO DATA #####
# read ensembles and prior ncdf4 objects
nc.post = nc_open('/Users/trevh/Research/assimilation-cfr/data/tas_ens_da_hydro_r.1000-2000_d.16-Feb-2018.nc')
nc.prior = nc_open('/Users/trevh/Research/assimilation-cfr/data/tas_prior_da_hydro_r.1000-2000_d.16-Feb-2018.nc')

# read ensembles and prior ncdf4 objects
# nc.post = nc_open('/Users/Trevor/research/assimilation-cfr/data/tas_ens_da_hydro_r.1000-2000_d.16-Feb-2018.nc')
# nc.prior = nc_open('/Users/Trevor/research/assimilation-cfr/data/tas_prior_da_hydro_r.1000-2000_d.16-Feb-2018.nc')


# commonly use 17x30 and 10x15
lat.basis = 32
lon.basis = 28

basis = create_basis(lat.basis, lon.basis, nc.prior)
proj = solve(t(basis) %*% basis) %*% t(basis)

# prep prior
prior.sub = read.csv("data/prior_ens.txt", header = F)
prior.sub = as.vector(prior.sub[,1])

prior.ens = prep_prior(nc.prior)
prior.ens = prior.ens[,,prior.sub]


#  example of approximation
library(plotly)
prior = basis %*% (proj %*% as.vector(prior.ens[,,1]))
prior = matrix(prior, 96, 144)

plot_ly(showscale = FALSE) %>% 
   add_surface(z = ~prior.ens[,,2])

plot_ly(showscale = FALSE) %>% 
  add_surface(z = ~prior)


post = prep_post_ens(nc.post, 100)
post.coef = sapply(1:dim(post)[3], function(x) proj %*% as.vector(post[,,x]))

coef.rng = 550:575
plot(post.coef[coef.rng, 1], type = "l", xlab = "Coefficient", ylab = "Value")
for (i in 2:100) {
  lines(post.coef[coef.rng,i])
}


post = prep_post_time(nc.post, 10)
prior.coef = proj %*% as.vector(prior.ens[,,10])
post.coef = sapply(1:dim(post)[3], function(x) proj %*% as.vector(post[,,x]))

coef.rng = 550:575
plot(post.coef[coef.rng, 1], type = "l", xlab = "Coefficient", ylab = "Value")
for (i in 2:100) {
  lines(post.coef[coef.rng,i])
}
lines(prior.coef[coef.rng],col="red")

# What if we used the SRFS??
# its actually very different looking. Might fit ED assumptions better
prior.srsf = sign(diff(prior.coef, 1)) * sqrt(abs(diff(prior.coef, 1)))
post.srsf = sign(diff(post.coef, 1)) * sqrt(abs(diff(post.coef, 1)))

prior.srsf1 = sqrt(abs(diff(prior.coef, 1)))
post.srsf1 = sqrt(abs(diff(post.coef, 1)))

 # prior.srsf = sqrt(abs(diff(prior.coef, 1)))
# post.srsf = sqrt(abs(diff(post.coef, 1)))

plot(post.srsf[coef.rng, 1], type = "l", xlab = "Coefficient", ylab = "Value")
for (i in 2:100) {
  lines(post.srsf[coef.rng,i])
}
lines(prior.srsf[coef.rng],col="red")

plot(post.srsf1[coef.rng, 1], type = "l", xlab = "Coefficient", ylab = "Value")
for (i in 2:100) {
  lines(post.srsf1[coef.rng,i])
}
lines(prior.srsf1[coef.rng],col="red")

