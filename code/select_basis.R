# source('code/setup.R')

##### CONNECT TO DATA #####
# read ensembles and prior ncdf4 objects
nc.post = nc_open('/Users/Trevor/research/assimilation-cfr/data/tas_ens_da_hydro_r.1000-2000_d.16-Feb-2018.nc')
nc.prior = nc_open('/Users/Trevor/research/assimilation-cfr/data/tas_prior_da_hydro_r.1000-2000_d.16-Feb-2018.nc')


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
  prior = prior - rowMeans(prior)
  
  # normalize
  lats = as.vector(nc.prior$dim$lat$vals)
  latmat = matrix(rep(lats, n.lon), n.lat, n.lon)
  latmat = sqrt(abs(cos(latmat*pi/180)))
  
  prior = prior * latmat
  prior.avg = prior.avg + prior
}
prior.avg = prior.avg / n.time



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


# # MSE over lat conditional on longittude
# plot(MSE[1,], type = "l", ylim = c(20000, 100000))
# for(i in 2:nrow(MSE)) {
#   lines(MSE[i,])
# }


# # reshape for plotting
# rotate <- function(x) t(apply(x, 2, rev))
# mse.plot = rotate(MSE)
# bic.plot = rotate(bic.field)

mse.plot = MSE
bic.plot = bic.field

which(bic.field == min(bic.field), arr.ind = T) + c(19, 19)

# mse
dimnames(mse.plot) = list(lon.range, lat.range)
mse.gg = melt(mse.plot)
colnames(mse.gg) = c("lon", "lat", "value")

ggplot(mse.gg, aes(lon, lat, z = value, fill = value)) + 
  geom_raster(interpolate = TRUE) +
  geom_contour(color = "grey30") +
  # geom_point(aes(x=22, y=20), colour="black") +
  # geom_point(aes(x=24, y=24), colour="black") +
  # geom_point(aes(x=26, y=28), colour="black") +
  geom_point(aes(x=28, y=32), colour="black") +
  scale_fill_distiller(palette = "RdBu") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(title = "MSE")

ggsave(paste0("paper/figures/mse_select", ".png"), width = 5, height = 3.2)



# bic
dimnames(bic.plot) = list(lon.range, lat.range)
bic.gg = melt(bic.plot)
colnames(bic.gg) = c("lon", "lat", "value")

ggplot(bic.gg, aes(lon, lat, z = value, fill = value)) + 
  geom_raster(interpolate = TRUE) +
  geom_contour(color = "grey30") +
  # geom_point(aes(x=22, y=20), colour="black") +
  # geom_point(aes(x=24, y=24), colour="black") +
  # geom_point(aes(x=26, y=28), colour="black") +
  geom_point(aes(x=28, y=32), colour="black") +
  scale_fill_distiller(palette = "RdBu") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(title = "BIC")

ggsave(paste0("paper/figures/bic_select", ".png"), width = 5, height = 3.2)


which(bic.plot == min(bic.plot), arr.ind = T) + c(19, 19)

# MSE surface
plot_ly(showscale = FALSE) %>% 
  add_surface(z = ~MSE) %>%
  layout(
    scene = list(
      xaxis = list(title = "lon"),
      yaxis = list(title = "lat"),
      zaxis = list(title = "MSE")
    )
  )

# runtime surface
plot_ly(showscale = FALSE) %>% 
  add_surface(z = ~bic.field) %>%
  layout(
    scene = list(
      xaxis = list(title = "lon"),
      yaxis = list(title = "lat"),
      zaxis = list(title = "BIC")
    )
  )

# for the top 20 lowest MSE what are the lat and lon?
for(i in 0:20) {
  minN.mse = sort(MSE, decreasing = TRUE)[length(MSE) - i]
  minN.ind = which(MSE == minN.mse, arr.ind = TRUE)
  cat(" nlon = ", minN.ind[,2] + min(lon.range)-1, " nlat = ", minN.ind[,1] + min(lat.range)-1,  minN.mse, "\n")
}

min.mse = which(MSE == min(MSE), arr.ind = TRUE)
cat(" nlon = ", min.mse[,2] + min(lon.range)-1, "\n", "nlat = ", min.mse[,1] + min(lat.range)-1)
