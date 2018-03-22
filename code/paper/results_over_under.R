# source('code/setup.R')


##### CONNECT TO DATA #####
# read ensembles and prior ncdf4 objects
nc.post = nc_open('/Users/Trevor/research/assimilation-cfr/data/tas_ens_da_hydro_r.1000-2000_d.16-Feb-2018.nc')
nc.prior = nc_open('/Users/Trevor/research/assimilation-cfr/data/tas_prior_da_hydro_r.1000-2000_d.16-Feb-2018.nc')


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

# compute the central regions and find the ED of each prior ensemble member in the list
# also save the CR difference fields

lower = array(0, dim = c(32*28, 100))
upper = array(0, dim = c(32*28, 100))
priors = array(0, dim = c(32*28, 100))
for(e in 1:100) {
  # slice prior to member e
  prior = prior.ens[,,e]
  
  # prep posterior at member e
  post = prep_post_time(nc.post, e)
  
  ##### FIT BASIS AND FIND ED #####
  prior.coef =  proj %*% as.vector(prior)
  post.coef = sapply(1:dim(post)[3], function(x) proj %*% as.vector(post[,,x]))
  
  # calculate extremal depth
  ed = edepth_set(post.coef)
  central = central_region(post.coef, ed, 0.5)
  
  priors[,e] = prior.coef
  lower[,e] = central[[1]]
  upper[,e] = central[[2]]
}

under = sapply(1:896, function(x) sum(priors[x,] < lower[x,]))
over = sapply(1:896, function(x) sum(upper[x,] < priors[x,]))

under = matrix(basis %*% -under, nrow(prior), ncol(prior))
over = matrix(basis %*% over, nrow(prior), ncol(prior))
outer = over - under

field_plot(under, nc.prior, main = "Under the Central Region")
ggsave("paper/figures/under_cr.png", width = 5, height = 3.2)

field_plot(over, nc.prior, main = "Over the Central Region")
ggsave("paper/figures/over_cr.png", width = 5, height = 3.2)

field_plot(outer, nc.prior, main = "Outside the Central Region")
ggsave("paper/figures/outer_cr.png", width = 5, height = 3.2)

