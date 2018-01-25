rm(list = ls())
gc()

library(ncdf4)
library(dplyr)
library(plotly)
library(extdepth)
library(mgcv)
library(irlba)

# open connection to TAS file
nc = nc_open('/Users/trevh/Research/assimilation-cfr/data/tas_ens_da_hydro_r.1000-2000_d.30-Nov-2017.nc')

# import 1 ensemble for all time points
clime = ncvar_get(nc, attributes(nc$var)$names[1], count = c(50, 50, 50, 1))
clime = sapply(seq(dim(clime)[3]), function(x) as.vector(clime[ , , x]))
clime = t(clime)

clime.cov = cov(clime)
clime.eigen = partial_eigen(clime.cov, n = 10)

clime.eigen2 = svdr(clime, k=5)
clime.eigen2 = ssvd(clime)

clime.eigen3 = svd(clime.cov)$d
clime.eigen3 = eigen(clime.cov)

# partial eigen is much faster
microbenchmark::microbenchmark(partial_eigen(clime.cov), times = 10)
microbenchmark::microbenchmark(eigen(clime.cov), times = 10)

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

