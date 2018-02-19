# how far does the prior deviate from the central region
alpha.dev = rep(0, length(prior.alpha))
for (a in 1:length(alpha.dev)) {
  if (central[[1]][a] > prior.alpha[a]) {
    alpha.dev[a] = prior.alpha[a] - central[[1]][a]
  } 
  if (prior.alpha[a] > central[[2]][a]) {
    alpha.dev[a] = prior.alpha[a] - central[[2]][a]
  }
}

# reshape for plotting
alpha.dev = matrix(alpha.dev, nlat.centers, nlon.centers)
alpha.plot = t(apply(alpha.dev, 2, rev))

# create custom color ramp for positive and negative differences
colormap = colorRampPalette(c("red", "white", "blue"))(2001)

# plot how far coefficients are above or below the central region
image.plot(alpha.plot,
           col = colormap,
           breaks = seq(-1.001, 1, 0.001))

# smoothed coefficient contour map over the original space
basis.point = expand.grid(floor(seq(1, length(lats), length.out = nlat.centers)),
                          floor(seq(1, length(lons), length.out = nlon.centers)))
coef.map = prior*0
for (p in 1:length(alpha.dev)) {
  lat = basis.point[,1][p]
  lon = basis.point[,2][p]
  
  coef.map[floor(lat), floor(lon)] = alpha.dev[p]
}

coef.fit = fastLmPure(basis, as.vector(coef.map))$fitted.values
coef.fit = matrix(coef.fit, n.lat, n.lon)


map.plot = t(apply(coef.fit, 2, rev))

# create custom color ramp for positive and negative differences
colormap = colorRampPalette(c("red", "white", "blue"))(2001)

# plot how far coefficients are above or below the central region
image.plot(map.plot)
