library(ncdf4)
library(tictoc)
library(refund)

source("../research/assimilation-cfr/paper/code/depth_tests.R")

prep_prior = function(nc.prior) {
  
  n.lon = nc.prior$dim$lon$len
  n.lat = nc.prior$dim$lat$len
  n.ens = nc.prior$dim$time2$len
  
  # extract data from the ncdf4 objects
  prior = ncvar_get(nc.prior, attributes(nc.prior$var)$names[1], start = c(1, 1, 1), count = c(-1, -1, -1))
  
  # transpose for intuitive (to me) layout
  prior = aperm(prior, c(2, 1, 3))
  
  # remove lat means
  # prior = vapply(1:n.ens, function(x) prior[,,x] - rowMeans(prior[,,x]), FUN.VALUE = matrix(0, nrow = n.lat, ncol = n.lon))
  
  # normalize
  lats = as.vector(nc.prior$dim$lat$vals)
  latmat = matrix(rep(lats, n.lon), n.lat, n.lon)
  latmat = sqrt(abs(cos(latmat*pi/180)))
  
  prior = vapply(1:n.ens, function(x) prior[,,x]*latmat, FUN.VALUE = matrix(0, nrow = n.lat, ncol = n.lon))
  
  return(prior)
}
prep_post = function(nc.post, t) {
  
  n.lon = nc.post$dim$lon$len
  n.lat = nc.post$dim$lat$len
  n.ens = nc.post$dim$sub_ens$len
  
  # extract data from the ncdf4 objects
  ens = ncvar_get(nc.post, attributes(nc.post$var)$names[1], start = c(1, 1, t, 1), count = c(-1, -1, 1, -1))
  
  # transpose for intuitive (to me) layout
  ens = aperm(ens, c(2, 1, 3))
  
  # remove lat means
  # ens = vapply(1:n.ens, function(x) ens[,,x] - rowMeans(ens[,,x]), FUN.VALUE = matrix(0, nrow = n.lat, ncol = n.lon))
  
  # normalize
  lats = as.vector(nc.post$dim$lat$vals)
  latmat = matrix(rep(lats, n.lon), n.lat, n.lon)
  latmat = sqrt(abs(cos(latmat*pi/180)))
  
  ens = vapply(1:n.ens, function(x) ens[,,x]*latmat, FUN.VALUE = matrix(0, nrow = n.lat, ncol = n.lon))
  
  return(ens)
}
fadtest = function(f, g, k = 35) {
  ### FAD method
  h = cbind(f, g)
  fpc = fpca.face(t(h), knots = k)
  
  fscores = fpc$scores[1:ncol(f),]
  gscores = fpc$scores[(ncol(f) + 1):ncol(h),]
  
  vals = sapply(1:ncol(fscores), function(x) adk.test(fscores[,x], gscores[,x])$adk[2,])
  pvals = vals[2,]
  tvals = vals[1,]
  c(max(tvals), min(1, min(pvals)*ncol(fscores)))
}


# prior
nc.prior = nc_open('../research/proxy/data/tas_prior_da_hydro_r.1000-2000_d.16-Feb-2018.nc')

prior = prep_prior(nc.prior)
prior = prior[,,prior_ind]
prior = flatten(prior)

# post
nc.post = nc_open('../research/proxy/data/tas_ens_da_hydro_r.1000-2000_d.16-Feb-2018.nc')

fad = matrix(0, 998, 2)
for(i in 1:998) {
  tic(paste0("Year ", i))
  
  post = prep_post(nc.post, i)
  post = flatten(post)
  
  fad[i, ] = fadtest(prior, post)
  toc()
}


fadvals = data.frame(time = years, value = fad[,1])
ggplot(data = fadvals, aes(years, value)) +
  geom_point() +
  geom_smooth(color = "red", size = 2, se = F) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5), 
        text = element_text(size=20)) +
  ylab("FAD") +
  xlab("Year")
  # ggtitle("Average difference between prior and posterior means")
ggsave("../research/assimilation-cfr/paper/misc/fad_over_time.png", width = 8, heigh = 6)

plot(fad[,1], xlab = "Time", ylab = "FAD")
plot(fad[,2], xlab = "Time", ylab = "FAD pval")


h = cbind(prior, post)
fpc = fpca.face(t(prior), knots = k)

plot(fpc$evalues)

fscores = fpc$scores[1:ncol(f),]
gscores = fpc$scores[(ncol(f) + 1):ncol(h),]

boxplot(fscores)

f.ind = sample(1:998, size = 499)
f = prior[,f.ind]
g = prior[,-f.ind]
c(fadtest(f, g)[2], bandtest(f, g)[2], kolm(f, g)[2])
# # prior
# nc.prior = nc_open('../research/proxy/data/tas_prior_da_hydro_r.1000-2000_d.16-Feb-2018.nc')
# 
# prior = prep_prior(nc.prior)
# prior = prior[,,sample(1:998, 100)]
# prior = flatten(prior)
# prior.fpca = fpca.face(t(priorsub))
# 
# boxplot(prior.fpca$scores)
# 
# # prior
# nc.post = nc_open('../research/proxy/data/tas_ens_da_hydro_r.1000-2000_d.16-Feb-2018.nc')
# 
# post = prep_post(nc.post, 998)
# post = flatten(post)
# post.fpca = fpca.face(t(post))
# 
# boxplot(post.fpca$scores)
# 
# fadtest(prior, post)
