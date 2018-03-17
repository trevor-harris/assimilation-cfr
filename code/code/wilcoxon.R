# Wilcoxon for ED
source('code/setup.R')

##### CONNECT TO DATA #####
# read ensembles and prior ncdf4 objects
nc.post = nc_open('/Users/Trevor/research/assimilation-cfr/data/tas_ens_da_hydro_r.1000-2000_d.16-Feb-2018.nc')
nc.prior = nc_open('/Users/Trevor/research/assimilation-cfr/data/tas_prior_da_hydro_r.1000-2000_d.16-Feb-2018.nc')


# 28x32 selected by BIC procedure
lat.basis = 28
lon.basis = 32

basis = create_basis(lat.basis, lon.basis, nc.prior)
proj = solve(t(basis) %*% basis) %*% t(basis)

# prep prior
prior.sub = read.csv("data/prior_ens.txt", header = F)
prior.sub = as.vector(prior.sub[,1])

prior.ens = prep_prior(nc.prior)

# setup list of fields
z = prior.ens[,,-prior.sub]
x = prior.ens[,,prior.sub]

# switch to basis functions
z.coef = sapply(1:dim(z)[3], function(e) proj %*% as.vector(z[,,e]))
x.coef = sapply(1:dim(x)[3], function(e) proj %*% as.vector(x[,,e]))

# ready z.coef and z.depths for ED
z.coef = sort_by_ed(z.coef)
z.depths = depth_set(z.coef)

# find x.ed first
x.ed = edepth_multi_fast(x.coef, z.coef, z.depths)

wil = rep(0, 1001)
for (t in 1:1001) {
  y = prep_post_ens(nc.post, t)
  
  # switch to basis functions
  y.coef = sapply(1:dim(y)[3], function(e) proj %*% as.vector(y[,,e]))
  
  # calculate ed for each wrt to z.coef
  y.ed = edepth_multi_fast(y.coef, z.coef, z.depths)
  
  # calculate the ranks and extract the y ranks
  ranks = rank(c(x.ed, y.ed), ties.method = "average")
  ranks = ranks[length(x.ed):length(ranks)]
  
  nr = length(y.ed)
  wil[t] = sum(ranks)
}

# monte carlo distribution and associated p value of the w test
runs = 1000000
wdist = sapply(1:runs, function(t) sum(sample(1:200, 100, replace=FALSE)))
wil_p = sapply(1:1001, function(t) sum(wil[t] <= wdist)) / runs

plot(wil_p, main = "Wilcoxon p values for times 1 to 1001")
abline(0.05, 0)
write.csv(wil_p, "wilcoxon_pvalues")

# multiple comparisons correction
wil_fdr = p.adjust(wil_p, "fdr")

plot(wil_fdr)
abline(0.05, 0)
write.csv(wil_fdr, "wilcoxon_fdr")
