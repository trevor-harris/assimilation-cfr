rm(list = ls())
gc()

library(extdepth)
library(ncdf4)
library(dplyr)
library(reshape2)
library(plotly)

source('code/functions.R')

##### CONNECT TO DATA #####
# read ensembles and prior ncdf4 objects
nc.post = nc_open('/Users/trevh/Research/assimilation-cfr/data/tas_ens_da_hydro_r.1000-2000_d.16-Feb-2018.nc')
nc.prior = nc_open('/Users/trevh/Research/assimilation-cfr/data/tas_prior_da_hydro_r.1000-2000_d.16-Feb-2018.nc')



##### prep prior #####
prior.sub = read.csv("data/prior_ens.txt", header = F)
prior.sub = as.vector(prior.sub[,1])

prior.ens = prep_prior(nc.prior)
prior.ens = prior.ens[,,prior.sub]



##### prep posterior ####
post.ens = prep_post_ens(nc.post, 100)

# split em all
nlat = 24
nlon = 18
regions = 32
ens = dim(prior.ens)[3]

prior.split = vapply(1:ens, function(x) matsplitter(prior.ens[,,x], nlat, nlon),
                     FUN.VALUE = array(0,dim = c(nlat, nlon, regions)))

post.split = vapply(1:ens, function(x) matsplitter(post.ens[,,x], nlat, nlon),
                    FUN.VALUE = array(0,dim = c(nlat, nlon, regions)))


# find the permutation distribution of wilcoxon fields
perms = 100
perm.fields = matrix(0, regions, perms)
for (p in 1:perms) {
  new.fields = permute_fields(prior.split, post.split, seed = p)
  perm.fields[,p] = wilcox.field(new.fields[[1]], new.fields[[2]])
}

# find the observed wilcoxon field
wilco.field = wilcox.field(prior.split, post.split)


plot_ly(showscale = F) %>%
  add_surface(z = ~matrix(wilco.field, 8, 8))


# What if we averaged over the projected values and found a single test statistic
iter=5000
seed = 12345
nlat = dim(prior.split)[1]
nlon = dim(prior.split)[2]
regions = dim(prior.split)[3]
ens = dim(prior.split)[4]

prior.proj = matrix(0, ens, iter)
post.proj = matrix(0, ens, iter)
wilco.field = 1:regions
proj = matrix(rnorm(nlat*nlon*iter), nlat*nlon, iter)

for(r in 1:regions) {
  for(e in 1:ens) {
    prior.proj[e,] = as.vector(as.vector(prior.split[,,r,e]) %*% proj)
    post.proj[e,] = as.vector(as.vector(post.split[,,r,e]) %*% proj)
  }
  wilco.field[r] = wilcox.test(rowMeans(prior.proj), rowMeans(post.proj), alternative = "two.sided",
                               paired = TRUE)$statistic
  
  # wilco = sapply(1:iter, function(x) wilcox.test(prior.proj[,x], 
  #                                                post.proj[,x],
  #                                                alternative = "two.sided", 
  #                                                paired = TRUE)$statistic)
  # wilco.field[r] = mean(wilco)
}
# wilcox.test(prior.proj[,1], post.proj[,1])
wilco.field

plot_ly(showscale = F) %>%
  add_surface(z = ~matrix(wilco.field, 8, 8))







# plot the observed and permuted (as vectors)
plot(perm.fields[,1], type = "l", ylab = "Mean Wilcoxon", xlab = "Region", ylim = c(4400, 5800))
for(p in 2:perms) {
  lines(perm.fields[,p])
}
lines(wilco.field, col = "red", lwd = 2)



# find the central regions
perm.ed = edepth_set(perm.fields)
perm.cr = central_region(perm.fields, perm.ed)

# add them to the plot
lines(perm.cr[[1]], col = "blue", lwd = 2)
lines(perm.cr[[2]], col = "blue", lwd = 2)



# what are the most extreme functions?
perm.out = perm.fields[,perm.ed <= 0.05]
for(p in 1:5) {
  lines(perm.out[,p], col = "blue", lwd = 2)
}

# depth of the observed field
edepth(wilco.field, perm.fields)

