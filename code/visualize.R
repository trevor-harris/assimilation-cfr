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
field_plot(prior.diff, nc.prior, "T = 100 with 510 coef")


# Use outliers (50% interval * constant) to find outlying regions
alpha.dev = rep(0, length(prior.alpha))
for (a in 1:length(alpha.dev)) {
  if (outlier[[1]][a] > prior.alpha[a]) {
    alpha.dev[a] = prior.alpha[a] - outlier[[1]][a]
  } 
  if (prior.alpha[a] > outlier[[2]][a]) {
    alpha.dev[a] = prior.alpha[a] - outlier[[2]][a]
  }
}
prior.vec = basis %*% alpha.dev


prior.diff = matrix(prior.vec , n.lat, n.lon)
field_plot(prior.diff, nc.prior)



#### image.plot version ####
# # shift longitutde to -180 to 180 format
# lons.re = lons - 180
# 
# # overlay coeff diff field over the earth
# image.plot(x = lons.re, y = lats, prior.plot, zlim = c(-max(abs(prior.plot)), max(abs(prior.plot))))
# map("world", add = T)

#### OLD TRIAL METHODS ####
# # how far does the prior deviate from the central region
# alpha.dev = rep(0, length(prior.alpha))
# for (a in 1:length(alpha.dev)) {
#   if (central[[1]][a] > prior.alpha[a]) {
#     alpha.dev[a] = prior.alpha[a] - central[[1]][a]
#   } 
#   if (prior.alpha[a] > central[[2]][a]) {
#     alpha.dev[a] = prior.alpha[a] - central[[2]][a]
#   }
# }
# 
# # reshape for plotting
# alpha.dev = matrix(alpha.dev, nlat.centers, nlon.centers)
# alpha.plot = t(apply(alpha.dev, 2, rev))
# 
# # create custom color ramp for positive and negative differences
# colormap = colorRampPalette(c("blue", "lightblue", "white", "yellow", "red"))(2001)
# 
# # plot how far coefficients are above or below the central region
# image.plot(alpha.plot,
#            col = colormap,
#            breaks = seq(-1.001, 1, 0.001))
# 
# 
# 
# 
# # Method 1 - smoothed coefficient difference contour map over the original space
# basis.point = expand.grid(floor(seq(1, length(lats), length.out = nlat.centers)),
#                           floor(seq(1, length(lons), length.out = nlon.centers)))
# coef.map = prior*0
# for (p in 1:length(alpha.dev)) {
#   lat = basis.point[,1][p]
#   lon = basis.point[,2][p]
#   
#   coef.map[floor(lat), floor(lon)] = alpha.dev[p]
# }
# 
# coef.fit = fastLmPure(basis, as.vector(coef.map))$fitted.values
# coef.fit = matrix(coef.fit, n.lat, n.lon)
# map.plot = t(apply(coef.fit, 2, rev))
# 
# # plot how far coefficients are above or below the central region
# image.plot(map.plot)
# 
# 
# 
# 
# # Method 2 - plot the difference in the actual (smoothed) fields
# # map central regions back to fields
# lower = basis %*% central[[1]]
# upper = basis %*% central[[2]]
# 
# prior.vec = basis %*% prior.alpha
# 
# # how far does the prior deviate from the central region
# prior.dev = rep(0, length(prior.vec))
# for (a in 1:length(prior.dev)) {
#   if (lower[a] > prior.dev[a]) {
#     prior.dev[a] = prior.vec[a] - lower[a]
#   } 
#   if (prior.dev[a] > upper[a]) {
#     prior.dev[a] = prior.vec[a] - upper[a]
#   }
# }
# 
# # reshape for plotting
# prior.dev = matrix(prior.dev , n.lat, n.lon)
# prior.plot = t(apply(prior.dev, 2, rev))
# 
# image.plot(prior.plot)

# plot_ly() %>% add_surface(z = ~matrix(prior.vec , n.lat, n.lon))
# plot_ly() %>%
#   add_surface(z = ~matrix(lower , n.lat, n.lon)) %>%
#   add_surface(z = ~matrix(upper , n.lat, n.lon))
# plot_ly() %>% add_surface(z = ~matrix(prior.dev , n.lat, n.lon))
