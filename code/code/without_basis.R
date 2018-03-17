source('code/setup.R')

e = 100

# prep prior
prior.sub = read.csv("data/prior_ens.txt", header = F)
prior.sub = as.vector(prior.sub[,1])

prior.ens = prep_prior(nc.prior)
prior.ens = prior.ens[,,prior.sub]

# slice prior to member e
prior = prior.ens[,,e]

# prep posterior at member e
post = prep_post_time(nc.post, e)

# calculate extremal depth
post.coef = sapply(1:dim(post)[3], function(x) as.vector(post[,,x]))

#### after a million years this returned 0.787  at e =1-- indicating the prior is well within the posterior set
ed_nb = edepth(as.vector(prior), post.coef)

# Plots
# lower = matrix(central[[1]], nrow(prior), ncol(prior))
# upper = matrix(central[[2]], nrow(prior), ncol(prior))
# 
# prior.diff = matrix(prior.dev, nrow(prior), ncol(prior))
# 
# field_plot(lower, nc.prior, main = "lower")
# field_plot(upper, nc.prior, main = "upper")
# field_plot(prior.diff, nc.prior, main="no basis")


### Using coefficients we get an ED of 0
# commonly use 17x30 and 10x15
lat.basis = 30
lon.basis = 30

basis = create_basis(lat.basis, lon.basis, nc.prior)
proj = solve(t(basis) %*% basis) %*% t(basis)

# prep prior
prior.sub = read.csv("data/prior_ens.txt", header = F)
prior.sub = as.vector(prior.sub[,1])

prior.ens = prep_prior(nc.prior)
prior.ens = prior.ens[,,prior.sub]

# slice prior to member e
prior = prior.ens[,,e]

# prep posterior at member e
post = prep_post_time(nc.post, e)

prior.coef = proj %*% as.vector(prior)
post.coef = sapply(1:dim(post)[3], function(x) proj %*% as.vector(post[,,x]))

### return 0.872 at e=1
ed_b = edepth(prior.coef, post.coef)

# Use the x% central regions to measure outlyingness
ed = edepth_set(post.coef)
central = central_region(post.coef, ed, 0.05)
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

summary(as.vector(prior.diff))
