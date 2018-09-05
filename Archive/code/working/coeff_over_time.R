source('code/setup.R')

##### CONNECT TO DATA #####
# read ensembles and prior ncdf4 objects
nc.post = nc_open('/Users/trevh/Research/assimilation-cfr/data/tas_ens_da_hydro_r.1000-2000_d.16-Feb-2018.nc')
nc.prior = nc_open('/Users/trevh/Research/assimilation-cfr/data/tas_prior_da_hydro_r.1000-2000_d.16-Feb-2018.nc')


##### CONNECT TO DATA #####
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

# compute the central regions and find the ED of each prior ensemble member in the list
# also save the CR difference fields
ens_num = seq(0, 100, by=10)
ens_num[1] = 1

ens_num = c(10, 33, 44, 76)
# eds_sm = ens_num
cr_diff_fields = array(0, dim = c(dim(prior.ens)[1], dim(prior.ens)[2], length(ens_num)))
for(e in ens_num) {

  # slice prior to member e
  prior = prior.ens[,,e]
  
  # prep posterior at member e
  post = prep_post_time(nc.post, e)
  
  ##### FIT BASIS AND FIND ED #####
  prior.coef = proj %*% as.vector(prior)
  post.coef = sapply(1:dim(post)[3], function(x) proj %*% as.vector(post[,,x]))
  
  # calculate extremal depth
  ed = edepth_set(post.coef)
  central = central_region(post.coef, ed, 0.5)
  
  # Use the x% central regions to measure outlyingness
  alpha.dev = rep(0, length(prior.coef))
  for (a in 1:length(alpha.dev)) {
    if (central[[1]][a] > prior.coef[a]) {
      alpha.dev[a] = prior.coef[a] - central[[1]][a]
    } 
    if (prior.coef[a] > central[[2]][a]) {
      alpha.dev[a] = prior.coef[a] - central[[2]][a]
    }
  }
  prior.vec = basis %*% alpha.dev
  prior.diff = matrix(prior.vec , nrow(prior), ncol(prior))
  
  cr_diff_fields[,,which(ens_num == e)] = prior.diff
  # eds_sm[which(ens_num == e)] = edepth(prior.coef, post.coef)
}

# plot the difference fields
for(e in 1:length(ens_num)) {
  print(field_plot(cr_diff_fields[,,e], nc.prior, main = paste0("Ensemble ", ens_num[e]), zlim = c(-1, 1)))
  ggsave(paste0("paper/figures/ens_diff_", ens_num[e], ".png"), width = 5, height = 3.2)
}
