##### BASIS SELECTION #####

source('code/setup.R')

# get data dim
n.lon = nc.prior$dim$lon$len
n.lat = nc.prior$dim$lat$len
n.time = nc.prior$dim$time2$len

##### AVERAGE OVER ALL PRIOR FIELDS #####
prior.avg = matrix(0, n.lat, n.lon)
for (t in 1:n.time) {
  # extract data from the ncdf4 objects
  prior = ncvar_get(nc.prior, attributes(nc.prior$var)$names[1], start = c(1, 1, t), count = c(-1, -1, 1))
  
  # transpose for intuitive (to me) layout
  prior = t(prior)
  
  # remove lat means
  # prior = prior - rowMeans(prior)
  
  # normalize
  lats = as.vector(nc.prior$dim$lat$vals)
  latmat = matrix(rep(lats, n.lon), n.lat, n.lon)
  latmat = sqrt(abs(cos(latmat*pi/180)))
  
  prior = prior * latmat
  prior.avg = prior.avg + prior
}
prior.avg = prior.avg / n.time

plot_ly(showscale = FALSE) %>% 
  add_surface(z = ~prior.avg)

##### FIND OPTIMAL NUMBER OF BASIS FUNCTIONS #####
lon.range = 20:40
lat.range = 20:40

MSE = matrix(0, length(lon.range), length(lat.range))
bic.field = matrix(0, length(lon.range), length(lat.range))
for (nlo in lon.range) {
  for (nla in lat.range) {
    basis = create_basis(nla, nlo, nc.prior)
    
    # check the smoothed approximation
    prior.fit = lm(as.vector(prior.avg) ~ basis - 1)
    prior.smooth = matrix(prior.fit$fitted.values, 96, 144)
    
    MSE[nlo-min(lon.range)+1, nla-min(lat.range)+1] = sum(prior.fit$residuals^2)
    bic.field[nlo-min(lon.range)+1, nla-min(lat.range)+1] = BIC(prior.fit)
  }
}

mse.plot = MSE
bic.plot = bic.field

which(bic.field == min(bic.field), arr.ind = T) + c(19, 19)

# mse
dimnames(mse.plot) = list(lon.range, lat.range)
mse.gg = melt(mse.plot)
colnames(mse.gg) = c("lon", "lat", "value")

ggplot(mse.gg, aes(lon, lat, z = value, fill = value)) + 
  geom_raster(interpolate = T) +
  geom_contour(color = "grey30") +
  geom_point(aes(x=26, y=40), colour="black") +
  scale_fill_distiller(palette = "RdBu") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_blank()) +
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(x = "longitude knots",
       y = "latitude knots",
       title = "MSE", 
       fill="MSE")

ggsave(paste0("paper/figures/mse_select", ".png"), width = 5, height = 3.2)


# bic
dimnames(bic.plot) = list(lon.range, lat.range)
bic.gg = melt(bic.plot)
colnames(bic.gg) = c("lon", "lat", "value")

ggplot(bic.gg, aes(lon, lat, z = value, fill = value)) + 
  geom_raster(interpolate = T) +
  geom_contour(color = "grey30") +
  geom_point(aes(x=26, y=40), colour="black") +
  scale_fill_distiller(palette = "RdBu") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_blank()) +
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(x = "longitude knots",
       y = "latitude knots",
       title = "BIC", 
       fill="BIC")

ggsave(paste0("paper/figures/bic_select", ".png"), width = 5, height = 3.2)
