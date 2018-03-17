source('code/setup.R')

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

# compute the central regions and find the ED of each prior ensemble member in the list
# also save the CR difference fields
ens_num = seq(0, 100, by=20)
ens_num[1] = 1

# ens_num = 100
cr_diff_fields = array(0, dim = c(dim(prior.ens)[1], dim(prior.ens)[2], length(ens_num)))
eds = rep(0, length(ens_num))
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
  central = central_region(post.coef, ed, 0.05)
  
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
  eds[which(ens_num == e)] = edepth(prior.coef, post.coef)
}

# plot the difference fields
for(e in 1:length(ens_num)) {
  print(field_plot(cr_diff_fields[,,e], nc.prior, main = paste0("Ensemble ", ens_num[e])))
}




# commonly use 17x30 and 10x15
lat.basis = 17
lon.basis = 30

basis = create_basis(lat.basis, lon.basis, nc.prior)

# prep prior
prior.sub = read.csv("data/prior_ens.txt", header = F)
prior.sub = as.vector(prior.sub[,1])

prior = prep_prior(nc.prior)
prior = prior[,,prior.sub]
prior = prior[,,20]

# prep posterior (time slice)
post = prep_post_time(nc.post, 20)

##### FIT BASIS AND FIND ED #####
proj = solve(t(basis) %*% basis) %*% t(basis)

prior.coef = proj %*% as.vector(prior)
post.coef = sapply(1:dim(post)[3], function(x) proj %*% as.vector(post[,,x]))

# calculate extremal depth
ed = edepth_set(post.coef)
central = central_region(post.coef, ed, 0.05)


par(mfrow=c(2,1))
coef.range = 1:20
plot(post.coef[coef.range,1], type = "l")
for (i in 2:100) {
  lines(post.coef[coef.range,i])
}
lines(prior.coef[coef.range], col = "red")


coef.range = 1:20
plot(post.coef[coef.range,1001-20], type = "l")

for (i in (1001-20):1001) {
  lines(post.coef[coef.range,i])
}
lines(prior.coef[coef.range], col = "red")




