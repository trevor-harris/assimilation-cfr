# rm(list = ls())
# gc()

# years = 851:1848
# ens = 100
#
# # library(extdepth)
# library(ncdf4)
# library(dplyr)
# library(reshape2)
# library(ggplot2)
# # library(OpenImageR)
# library(future)
# library(future.apply)
# library(tictoc)
# # library(fields)
#
#
# world = map_data("world")
# world = world[world$long <= 178, ]
#
# # prepare data
# prep_prior = function(nc.prior) {
#
#   n.lon = nc.prior$dim$lon$len
#   n.lat = nc.prior$dim$lat$len
#   n.ens = nc.prior$dim$time2$len
#
#   # extract data from the ncdf4 objects
#   prior = ncvar_get(nc.prior, attributes(nc.prior$var)$names[1], start = c(1, 1, 1), count = c(-1, -1, -1))
#
#   # transpose for intuitive (to me) layout
#   prior = aperm(prior, c(2, 1, 3))
#
#   # remove lat means
#   # prior = vapply(1:n.ens, function(x) prior[,,x] - rowMeans(prior[,,x]), FUN.VALUE = matrix(0, nrow = n.lat, ncol = n.lon))
#
#   # normalize
#   lats = as.vector(nc.prior$dim$lat$vals)
#   latmat = matrix(rep(lats, n.lon), n.lat, n.lon)
#   latmat = sqrt(abs(cos(latmat*pi/180)))
#
#   prior = vapply(1:n.ens, function(x) prior[,,x]*latmat, FUN.VALUE = matrix(0, nrow = n.lat, ncol = n.lon))
#
#   return(prior)
# }
# prep_post = function(nc.post, t) {
#
#   n.lon = nc.post$dim$lon$len
#   n.lat = nc.post$dim$lat$len
#   n.ens = nc.post$dim$sub_ens$len
#
#   # extract data from the ncdf4 objects
#   ens = ncvar_get(nc.post, attributes(nc.post$var)$names[1], start = c(1, 1, t, 1), count = c(-1, -1, 1, -1))
#
#   # transpose for intuitive (to me) layout
#   ens = aperm(ens, c(2, 1, 3))
#
#   # remove lat means
#   # ens = vapply(1:n.ens, function(x) ens[,,x] - rowMeans(ens[,,x]), FUN.VALUE = matrix(0, nrow = n.lat, ncol = n.lon))
#
#   # normalize
#   lats = as.vector(nc.post$dim$lat$vals)
#   latmat = matrix(rep(lats, n.lon), n.lat, n.lon)
#   latmat = sqrt(abs(cos(latmat*pi/180)))
#
#   ens = vapply(1:n.ens, function(x) ens[,,x]*latmat, FUN.VALUE = matrix(0, nrow = n.lat, ncol = n.lon))
#
#   return(ens)
# }
# flatten = function(mat) {
#   matrix(mat, prod(dim(mat)[1:2]), dim(mat)[3])
# }
#
# # plots
# field_plot2 <- function(field, nc) {
#
#   lats = as.vector(nc$dim$lat$vals)
#   lons = as.vector(nc$dim$lon$vals)
#   dimnames(field) = list(lats, ifelse(lons >= 180, lons - 360, lons))
#
#   field.gg = melt(field)
#   colnames(field.gg) = c("lat", "lon", "Temp")
#
#   world = map_data("world")
#   world = world[world$long <= 178, ]
#
#   zlim = c(-max(abs(field)), max(abs(field)))
#
#   plt = ggplot() +
#     geom_raster(data = field.gg, aes(x=lon, y=lat, fill=Temp), interpolate = F) +
#     geom_polygon(data = world, aes(x=long, y=lat, group=group), fill = NA, color="black") +
#     coord_cartesian() +
#     scale_fill_distiller(palette = "RdBu") +
#     theme_void()
#   # labs(fill = "Anomaly")
#   # theme(legend.position="none")
#
#   print(plt)
# }
#
#
# ##### Actual data
# nc.post = nc_open('../research/proxy/data/tas_ens_da_hydro_r.1000-2000_d.16-Feb-2018.nc')
# nc.prior = nc_open('../research/proxy/data/tas_prior_da_hydro_r.1000-2000_d.16-Feb-2018.nc')
# prior_ind = read.csv("../research/proxy/data/prior-ens.txt", header = F)$V1
#
#
# prior = prep_prior(nc.prior)
# prior = prior[,,prior_ind]
#
#
#
# param = matrix(0, 998, 2)
# for(y in 1:998) {
#   tic()
#   post = prep_post(nc.post, y)
#
#   mu1 = matrix(0, nrow(prior), ncol(prior))
#   for(i in 1:nrow(prior)) {
#     for(j in 1:ncol(prior)) {
#       mu1[i, j] = mean(prior[i,j,])
#     }
#   }
#
#   mu2 = matrix(0, nrow(prior), ncol(prior))
#   for(i in 1:nrow(post)) {
#     for(j in 1:ncol(post)) {
#       mu2[i, j] = mean(post[i,j,])
#     }
#   }
#
#   sd1 = matrix(0, nrow(prior), ncol(prior))
#   for(i in 1:nrow(prior)) {
#     for(j in 1:ncol(prior)) {
#       sd1[i, j] = sd(prior[i,j,])
#     }
#   }
#
#   sd2 = matrix(0, nrow(prior), ncol(prior))
#   for(i in 1:nrow(post)) {
#     for(j in 1:ncol(post)) {
#       sd2[i, j] = sd(post[i,j,])
#     }
#   }
#
#
#   param[y,1] = mean((mu1 - mu2)^2)
#   param[y,2] = mean((sd1 / sd2)^2)
#
#   toc()
# }d

plot(param[,1])
plot(param[,2])

mu = data.frame(time = years, value = param[,1])
ggplot(data = mu, aes(years, value)) +
  geom_point() +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5), 
        text = element_text(size=20)) +
  ylab("Average difference squared") +
  xlab("Year")
  # ggtitle("Average difference between background and analysis means")
ggsave("../research/assimilation-cfr/paper/misc/means.png", width = 8, heigh = 6)


sig = data.frame(time = years, value = param[,2])
ggplot(data = sig, aes(years, value)) +
  geom_point() +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5), 
        text = element_text(size=20)) +
  ylab("Average ratio squared") +
  xlab("Year")
  # ggtitle("Average ratio of background over analysis standard deviations")
ggsave("../research/assimilation-cfr/paper/misc/sds.png", width = 8, heigh = 6)



