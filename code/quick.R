source('code/setup.R')

basis = create_basis(17, 30, nc.prior)
fields = preprocess(100, nc.ens, nc.prior)

ens = fields[["ens"]]
prior = fields[["prior"]]

##### FIT BASIS AND FIND ED #####
prior.alpha = coef(fastLmPure(basis, as.vector(prior)))
ens.alpha = sapply(1:dim(ens)[3], function(x) coef(fastLmPure(basis, as.vector(ens[,,x]))))

# calculate extremal depth
ed = sapply(1:ncol(ens.alpha), function(x) edepth(ens.alpha[,x], ens.alpha))
central = central_region(ens.alpha, ed, 0.05)

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

main = paste0("T = ", ) 
field_plot(prior.diff, nc.prior, "T = 100 with 510 coef")
