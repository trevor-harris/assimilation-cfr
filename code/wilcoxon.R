# Wilcoxon for ED

source('code/setup.R')

# commonly use 17x30 and 10x15
lat.basis = 5
lon.basis = 5

basis = create_basis(lat.basis, lon.basis, nc.prior)
proj = solve(t(basis) %*% basis) %*% t(basis)

# prep prior
prior.sub = read.csv("data/prior_ens.txt", header = F)
prior.sub = as.vector(prior.sub[,1])

prior.ens = prep_prior(nc.prior)

# setup list of fields
z = prior.ens[,,-prior.sub]
x = prior.ens[,,prior.sub]
y = prep_post_ens(nc.post, 1)

# subset for testing
z = z[,,1:20]
x = x[,,1:10]
y = y[,,1:10]

# switch to basis functions
z.coef = sapply(1:dim(z)[3], function(e) proj %*% as.vector(z[,,e]))
x.coef = sapply(1:dim(x)[3], function(e) proj %*% as.vector(x[,,e]))
y.coef = sapply(1:dim(y)[3], function(e) proj %*% as.vector(y[,,e]))

# smooth em
x.poly = predict(loess(x.coef[,1] ~ design[,1]))
plot(x.coef[,1], type = "l")
lines(x.poly, col ="blue")

# calculate ed for each wrt to z.coef
z.ed = edepth_set(z.coef)
x.ed = sapply(1:ncol(x.coef), function(e) edepth(x.coef[,e], z.coef))
y.ed = sapply(1:ncol(y.coef), function(e) edepth(y.coef[,e], z.coef))

plot(x.coef[,1], type = "l")
for (e in 2:10) {lines(x.coef[,e])}
for (e in 1:20) {lines(z.coef[,e], col="red")}
