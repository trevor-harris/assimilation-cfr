rm(list = ls())
gc()

library(ncdf4)
library(dplyr)
library(plotly)
library(extdepth)
library(mgcv)

# open connection to TAS file
nc = nc_open('/Users/trevh/Research/assimilation-cfr/data/tas_ens_da_hydro_r.1000-2000_d.30-Nov-2017.nc')

# import 1 ensemble for all time points
clime = ncvar_get(nc, attributes(nc$var)$names[1], count = c(-1, -1, -1, 1))
clime = sapply(seq(dim(clime)[3]), function(x) as.vector(clime[ , , x]))
clime = t(clime)

# import and convert data (all ensembles for 1 time point)
tas = ncvar_get(nc, attributes(nc$var)$names[1], count = c(-1, -1, 1, -1))
tas = lapply(seq(dim(tas)[3]), function(x) as.matrix(tas[ , , x]))


# smooth surfaces
temp = tas[[1]]
df <- data.frame(x = rep(seq_len(ncol(temp)), each = nrow(temp)),
                 y = rep(seq_len(nrow(temp)), times = ncol(temp)),
                 z = c(temp))

mod = gam(z ~ te(x, y, k=15), data = df)$fitted
temp.smooth = matrix(mod, nrow(temp), ncol(temp))
temp.diff = matrix(temp - mod, nrow(temp), ncol(temp))

# raw
plot_ly(showscale = FALSE) %>% 
  add_surface(z = ~temp)

# smooth
plot_ly(showscale = FALSE) %>% 
  add_surface(z = ~temp.smooth)

# difference
plot_ly(showscale = FALSE) %>% 
  add_surface(z = ~temp.diff)

# test
mod2 = gam(z ~ s(x, y, k=200, bs="tp"), data = df)
temp.smooth2 = matrix(mod2$fitted, nrow(temp), ncol(temp))
temp.diff2 = matrix(temp - mod2$fitted, nrow(temp), ncol(temp))


plot_ly(showscale = FALSE) %>% 
  add_surface(z = ~temp.smooth2)

# find basis coefficients

# apply extdepth to basis coef

