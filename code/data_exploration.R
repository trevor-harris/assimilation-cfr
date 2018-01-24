rm(list = ls())
gc()

library(ncdf4)
library(dplyr)
library(plotly)
library(extdepth)
library(mgcv)

# open connection to TAS file
nc = nc_open('/Users/trevh/Research/assimilation-cfr/data/tas_ens_da_hydro_r.1000-2000_d.30-Nov-2017.nc')

# import and convert data
tas = ncvar_get(nc, attributes(nc$var)$names[1], count = c(-1, -1, 1, -1))
tas = lapply(seq(dim(tas)[3]), function(x) as.matrix(tas[ , , x]))


# smooth surfaces
m = tas[[1]]
df <- data.frame(x = rep(seq_len(ncol(m)), each = nrow(m)),
                 y = rep(seq_len(nrow(m)), times = ncol(m)),
                 z = c(m))

# mod = gam(z ~ te(x, y), data = df)
mod <- matrix(gam(z ~ te(x, y), data = df)$fitted,
              nrow(m),
              ncol(m))

plot_ly(showscale = FALSE) %>% 
  add_surface(z = ~tas[[1]])

plot_ly(showscale = FALSE) %>% 
  add_surface(z = ~mod)

# find basis coefficients

# apply extdepth to basis coef

