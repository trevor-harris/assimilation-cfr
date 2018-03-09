source('code/setup.R')

t = 400

prior.sub = read.csv("data/prior_ens.txt", header = F)
prior.sub = as.vector(prior.sub[,1])

prior = prep_prior(nc.prior)
prior = prior[,,prior.sub]
prior = apply(prior, 1:2, mean)

# commonly use 17x30 and 10x15
lat.basis = 10
lon.basis = 15

basis = create_basis(lat.basis, lon.basis, nc.prior)
fields = preprocess(t, nc.ens, nc.prior)

ens = fields[["ens"]]

##### FIT BASIS AND FIND ED #####
prior.alpha = coef(fastLmPure(basis, as.vector(prior)))
ens.alpha = sapply(1:dim(ens)[3], function(x) coef(fastLmPure(basis, as.vector(ens[,,x]))))

# calculate extremal depth
ed = sapply(1:ncol(ens.alpha), function(x) edepth(ens.alpha[,x], ens.alpha))
central = central_region(ens.alpha, ed, 0.005)

# Use the x% central regions to measure outlyingness
alpha.dev = rep(0, length(prior.alpha))
for (a in 1:length(alpha.dev)) {
  if (central[[1]][a] > prior.alpha[a]) {
    alpha.dev[a] = prior.alpha[a] - central[[1]][a]
  } 
  if (prior.alpha[a] > central[[2]][a]) {
    alpha.dev[a] = prior.alpha[a] - central[[2]][a]
  }
}
prior.vec = basis %*% alpha.dev

prior.diff = matrix(prior.vec , nrow(prior), ncol(prior))

main = paste0(lat.basis*lon.basis, " basis functions") 
field_plot(prior.diff, nc.prior, main, c(-1, 1))
ggsave("plots/150coef.png", width = 8, height = 5.2)


# commonly use 17x30 and 10x15
lat.basis = 12
lon.basis = 20

basis = create_basis(lat.basis, lon.basis, nc.prior)
fields = preprocess(t, nc.ens, nc.prior)

ens = fields[["ens"]]

##### FIT BASIS AND FIND ED #####
prior.alpha = coef(fastLmPure(basis, as.vector(prior)))
ens.alpha = sapply(1:dim(ens)[3], function(x) coef(fastLmPure(basis, as.vector(ens[,,x]))))

# calculate extremal depth
ed = sapply(1:ncol(ens.alpha), function(x) edepth(ens.alpha[,x], ens.alpha))
central = central_region(ens.alpha, ed, 0.005)

# Use the x% central regions to measure outlyingness
alpha.dev = rep(0, length(prior.alpha))
for (a in 1:length(alpha.dev)) {
  if (central[[1]][a] > prior.alpha[a]) {
    alpha.dev[a] = prior.alpha[a] - central[[1]][a]
  } 
  if (prior.alpha[a] > central[[2]][a]) {
    alpha.dev[a] = prior.alpha[a] - central[[2]][a]
  }
}
prior.vec = basis %*% alpha.dev

prior.diff = matrix(prior.vec , nrow(prior), ncol(prior))

main = paste0(lat.basis*lon.basis, " basis functions") 
field_plot(prior.diff, nc.prior, main, c(-1, 1))
ggsave("plots/240coef.png", width = 8, height = 5.2)

d2 = edepth(prior.alpha, ens.alpha)

# commonly use 17x30 and 10x15
lat.basis = 15
lon.basis = 25

basis = create_basis(lat.basis, lon.basis, nc.prior)
fields = preprocess(t, nc.ens, nc.prior)

ens = fields[["ens"]]

##### FIT BASIS AND FIND ED #####
prior.alpha = coef(fastLmPure(basis, as.vector(prior)))
ens.alpha = sapply(1:dim(ens)[3], function(x) coef(fastLmPure(basis, as.vector(ens[,,x]))))

# calculate extremal depth
ed = sapply(1:ncol(ens.alpha), function(x) edepth(ens.alpha[,x], ens.alpha))
central = central_region(ens.alpha, ed, 0.005)

# Use the x% central regions to measure outlyingness
alpha.dev = rep(0, length(prior.alpha))
for (a in 1:length(alpha.dev)) {
  if (central[[1]][a] > prior.alpha[a]) {
    alpha.dev[a] = prior.alpha[a] - central[[1]][a]
  } 
  if (prior.alpha[a] > central[[2]][a]) {
    alpha.dev[a] = prior.alpha[a] - central[[2]][a]
  }
}
prior.vec = basis %*% alpha.dev

prior.diff = matrix(prior.vec , nrow(prior), ncol(prior))

main = paste0(lat.basis*lon.basis, " basis functions") 
field_plot(prior.diff, nc.prior, main, c(-1, 1))
ggsave("plots/375coef.png", width = 8, height = 5.2)

d3 = edepth(prior.alpha, ens.alpha)

# commonly use 17x30 and 10x15
lat.basis = 17
lon.basis = 30

basis = create_basis(lat.basis, lon.basis, nc.prior)
fields = preprocess(t, nc.ens, nc.prior)

ens = fields[["ens"]]

##### FIT BASIS AND FIND ED #####
prior.alpha = coef(fastLmPure(basis, as.vector(prior)))
ens.alpha = sapply(1:dim(ens)[3], function(x) coef(fastLmPure(basis, as.vector(ens[,,x]))))

# calculate extremal depth
ed = sapply(1:ncol(ens.alpha), function(x) edepth(ens.alpha[,x], ens.alpha))
central = central_region(ens.alpha, ed, 0.005)

# Use the x% central regions to measure outlyingness
alpha.dev = rep(0, length(prior.alpha))
for (a in 1:length(alpha.dev)) {
  if (central[[1]][a] > prior.alpha[a]) {
    alpha.dev[a] = prior.alpha[a] - central[[1]][a]
  } 
  if (prior.alpha[a] > central[[2]][a]) {
    alpha.dev[a] = prior.alpha[a] - central[[2]][a]
  }
}
prior.vec = basis %*% alpha.dev

prior.diff = matrix(prior.vec , nrow(prior), ncol(prior))

main = paste0(lat.basis*lon.basis, " basis functions") 
field_plot(prior.diff, nc.prior, main, c(-1, 1))
ggsave("plots/510coef.png", width = 8, height = 5.2)

d4 = edepth(prior.alpha, ens.alpha)


