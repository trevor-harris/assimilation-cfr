rm(list = ls())
gc()

library(ncdf4)
library(dplyr)
library(plotly)
library(mgcv)
library(irlba)

# preprocess
ens_preprocess <- function(nc, t) {
  
  # import 1 ensemble for all time points
  field = ncvar_get(nc, attributes(nc$var)$names[1], count = c(-1, -1, t, -1))
  field = lapply(1:dim(field)[3], function(x) t(as.matrix(field[ , , x])))
  
  lats = as.vector(nc$dim$lat$vals)
  nlon = nc$dim$lon$len
  nlat = nc$dim$lat$len
  
  latmat = matrix(rep(lats, nlon), nlat, nlon)
  latmat = sqrt(abs(cos(latmat*pi/180)))
  
  for (c in 1:length(field)) {
    field[[c]] = data.frame(x = rep(seq_len(nlon), each = nlat),
                            y = rep(seq_len(nlat), times = nlon),
                            z = c(field[[c]]))
    # gam smoother
    mod = gam(z ~ te(x, y), data = field[[c]])$fitted
    
    # return to matrix format
    field[[c]] = matrix(mod, nlat, nlon)
    
    # remove lat mean
    field.t = field[[c]] - rowMeans(field[[c]])
    
    # cosine normalize
    field[[c]] = field.t*latmat
  }
  return(field)
}


# open connection to TAS file
nc = nc_open('/Users/trevh/Research/assimilation-cfr/data/tas_ens_da_hydro_r.1000-2000_d.30-Nov-2017.nc')

ens = ens_preprocess(nc, 1)

plot_ly(showscale = FALSE) %>% 
  add_surface(z = ~ens[[1]])


