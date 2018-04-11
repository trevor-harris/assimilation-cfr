#### DEFINE USEFUL FUNCTIONS ####
create_basis = function(nlat.centers, nlon.centers, nc) {
  
  n.center = nlat.centers * nlon.centers
  
  lats = as.vector(nc$dim$lat$vals)
  lons = as.vector(nc$dim$lon$vals)
  n.lon = nc$dim$lon$len
  n.lat = nc$dim$lat$len
  
  # lats.cen = seq(min(lats), max(lats), length.out = nlat.centers)
  lons.cen = seq(min(lons), max(lons), length.out = nlon.centers)
  
  # cosine normalize to add more basis towards tropics and away from poles
  lats.pi = cos((lats / max(lats)) * (pi/2))
  lats.cen = seq(min(-lats.pi), max(lats.pi), length.out = nlat.centers)
  lats.cen = (acos(lats.cen) / (pi/2) * max(lats)) - 90
  lats.cen = -lats.cen
  
  
  # find the radius (1.5 times the minimal distance)
  # will either be the distance b/w the first point and the next (horizontally or vertically)
  rad = 1.5 * min(abs(lats.cen[1] - lats.cen[2]), abs(lons.cen[1] - lons.cen[2]))
  
  latlon = as.matrix(expand.grid(1:nlat.centers, 1:nlon.centers))
  basis = matrix(0, n.lat*n.lon, n.center)
  
  for (i in 1:n.center) {
    
    # dist = (lat - center.lat)^2 + (lon - center.lon)^2
    dist = c(outer((lats - lats.cen[latlon[i,][1]])^2, (lons - lons.cen[latlon[i,][2]])^2, "+"))
    
    # use the bisquare (radial) basis
    basis[,i] = (1 - dist/rad^2)^2 * (sqrt(dist) < rad)
  }
  
  return(basis)
}

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

prep_post_ens = function(nc.post, t) {
  
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


prep_post_time = function(nc.post, e) {
  
  n.lon = nc.post$dim$lon$len
  n.lat = nc.post$dim$lat$len
  n.time = nc.post$dim$time$len
  
  # extract data from the ncdf4 objects
  ens = ncvar_get(nc.post, attributes(nc.post$var)$names[1], start = c(1, 1, 1, e), count = c(-1, -1, -1, 1))
  
  # transpose for intuitive (to me) layout
  ens = aperm(ens, c(2, 1, 3))
  
  # remove lat means
  # ens = vapply(1:n.time, function(x) ens[,,x] - rowMeans(ens[,,x]), FUN.VALUE = matrix(0, nrow = n.lat, ncol = n.lon))
  
  # normalize
  lats = as.vector(nc.post$dim$lat$vals)
  latmat = matrix(rep(lats, n.lon), n.lat, n.lon)
  latmat = sqrt(abs(cos(latmat*pi/180)))
  
  ens = vapply(1:n.time, function(x) ens[,,x]*latmat, FUN.VALUE = matrix(0, nrow = n.lat, ncol = n.lon))
  
  return(ens)
}

# Plot smooth tilemap with contours
field_plot <- function(field, nc, main = "", zlim = c(-max(abs(field)), max(abs(field)))) {
  
  lats = as.vector(nc$dim$lat$vals)
  lons = as.vector(nc$dim$lon$vals)
  
  n.lon = nc$dim$lon$len
  n.lat = nc$dim$lat$len
  
  # rotates the matrix
  # field = t(apply(field, 2, rev))
  dimnames(field) = list(lats, lons-180)
  
  field.gg = melt(field)
  colnames(field.gg) = c("lat", "lon", "value")
  
  world = map_data("world")
  world = world[world$long <= 178, ]
  
  ggplot() +
    geom_raster(data = field.gg, aes(x=lon, y=lat, fill=value), interpolate = TRUE) +
    geom_contour(data = field.gg, aes(x=lon, y=lat, z=value)) +
    geom_polygon(data = world, aes(x=long, y=lat, group=group), fill = NA, color="black") +
    coord_cartesian() +
    scale_fill_distiller(palette = "RdBu", limits = zlim) +
    theme_void() +
    theme(plot.title = element_text(hjust = 0.5)) +
    labs(title = main)
}


#### regionalize ####
matsplitter<-function(M, r, c) {
  rg <- (row(M)-1)%/%r+1
  cg <- (col(M)-1)%/%c+1
  rci <- (rg-1)*max(cg) + cg
  N <- prod(dim(M))/r/c
  cv <- unlist(lapply(1:N, function(x) M[rci==x]))
  dim(cv)<-c(r,c,N)
  cv
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

#### permutation test ####

# permute_fields = function(prior.split, post.split, seed = 12345) {
#   set.seed(seed)
#   
#   nlat = dim(prior.split)[1]
#   nlon = dim(prior.split)[2]
#   regions = dim(prior.split)[3]
#   ens = dim(prior.split)[4]
#   
#   prior.split.new = prior.split
#   post.split.new = post.split
#   
#   prior.ind = sample.int(nlat * nlon)
#   
#   for(r in 1:regions) {
#     for (e in 1:ens) {
#       both.split = rbind(prior.split[,,r,e], post.split[,,r,e])
#       prior.split.new[,,r,e] = matrix(as.vector(both.split)[prior.ind], nlat, nlon)
#       post.split.new[,,r,e] = matrix(as.vector(both.split)[-prior.ind], nlat, nlon)
#     }
#   }
#   
#   return(list(prior.split.new, post.split.new))
# }

permute_fields = function(prior.split, post.split, seed = 12345) {
  set.seed(seed)
  
  nlat = dim(prior.split)[1]
  nlon = dim(prior.split)[2]
  regions = dim(prior.split)[3]
  ens = dim(prior.split)[4]
  
  prior.split.new = prior.split
  post.split.new = post.split
  
  prior.ind = sample(c(0, 1), size = ens, replace = TRUE)
  
  for (e in 1:ens) {
    if (prior.ind[e] == 1) {
      prior.split.new[,,,e] = post.split[,,,e]
      post.split.new[,,,e] = prior.split[,,,e]
    } else {
      prior.split.new[,,,e] = prior.split[,,,e]
      post.split.new[,,,e] = post.split[,,,e]
    }
  }
  
  return(list(prior.split.new, post.split.new))
}


