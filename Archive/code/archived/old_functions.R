# main function for preprocessinfg t
preprocess = function(t, nc.post, nc.prior) {
  
  n.lon = nc.ens$dim$lon$len
  n.lat = nc.ens$dim$lat$len
  n.ens = nc.ens$dim$sub_ens$len
  
  # extract data from the ncdf4 objects
  ens = ncvar_get(nc.ens, attributes(nc.ens$var)$names[1], start = c(1, 1, t, 1), count = c(-1, -1, 1, -1))
  prior = ncvar_get(nc.prior, attributes(nc.prior$var)$names[1], start = c(1, 1, t), count = c(-1, -1, 1))
  
  # transpose for intuitive (to me) layout
  ens = aperm(ens, c(2, 1, 3))
  prior = t(prior)
  
  # remove lat means
  # ens = vapply(1:n.ens, function(x) ens[,,x] - rowMeans(ens[,,x]), FUN.VALUE = matrix(0, nrow = n.lat, ncol = n.lon))
  # prior = prior - rowMeans(prior)
  
  # normalize
  lats = as.vector(nc.ens$dim$lat$vals)
  latmat = matrix(rep(lats, n.lon), n.lat, n.lon)
  latmat = sqrt(abs(cos(latmat*pi/180)))
  
  ens = vapply(1:n.ens, function(x) ens[,,x]*latmat, FUN.VALUE = matrix(0, nrow = n.lat, ncol = n.lon))
  prior = prior * latmat
  
  return(list(ens = ens, prior = prior))
}

#### Repeat for each region #### -- Outer most testing (this is the full procedure)
wilcox.field = function(prior.split, post.split, iter=500, seed = 12345) {
  set.seed(seed)
  
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
    wilco = sapply(1:iter, function(x) wilcox.test(prior.proj[,x], post.proj[,x], alternative = "two.sided",paired = TRUE)$statistic)
    wilco.field[r] = mean(wilco)
  }
  wilco.field
}


ks.field = function(prior.split, post.split, iter=500, seed = 12345) {
  set.seed(seed)
  
  nlat = dim(prior.split)[1]
  nlon = dim(prior.split)[2]
  regions = dim(prior.split)[3]
  ens = dim(prior.split)[4]
  
  prior.proj = matrix(0, ens, iter)
  post.proj = matrix(0, ens, iter)
  kol.field = 1:regions
  proj = matrix(rnorm(nlat*nlon*iter), nlat*nlon, iter)
  
  for(r in 1:regions) {
    for(e in 1:ens) {
      prior.proj[e,] = as.vector(as.vector(prior.split[,,r,e]) %*% proj)
      post.proj[e,] = as.vector(as.vector(post.split[,,r,e]) %*% proj)
    }
    kol = sapply(1:iter, function(x) ks.test(prior.proj[,x], post.proj[,x])$statistic)
    kol.field[r] = mean(kol)
  }
  kol.field
}