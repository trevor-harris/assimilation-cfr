rm(list = ls())
gc()

library(ncdf4)
library(extdepth)
library(mgcv)
library(irlba)

# open connection to TAS file
nc = nc_open('/Users/trevh/Research/assimilation-cfr/data/tas_ens_da_hydro_r.1000-2000_d.30-Nov-2017.nc')

# import 1 ensemble for all time points
clime = ncvar_get(nc, attributes(nc$var)$names[1], count = c(-1, -1, -1, 1))
clime = sapply(seq(dim(clime)[3]), function(x) as.vector(clime[ , , x]))
clime = t(clime)

clime.cov = cov(clime)
cl.eig = irlba(clime.cov, nv=5)

# subset down to reduce the dimension of the covariance matrix
clime.sm = clime[,seq(from=1, to=ncol(clime), by=2)]
clime.sm.cov = cov(clime.sm)
clsm.eig0 = irlba(clime.sm.cov, nv=5)

# subset down to reduce the dimension of the covariance matrix
clime.sm = clime[,seq(from=1, to=ncol(clime), by=5)]
clime.sm.cov = cov(clime.sm)
clsm.eig1 = irlba(clime.sm.cov, nv=5)

# subset down to reduce the dimension of the covariance matrix
clime.sm = clime[,seq(from=1, to=ncol(clime), by=10)]
clime.sm.cov = cov(clime.sm)
clsm.eig2 = irlba(clime.sm.cov, nv=5)

# subset down to reduce the dimension of the covariance matrix
clime.sm = clime[,seq(from=1, to=ncol(clime), by=20)]
clime.sm.cov = cov(clime.sm)
clsm.eig3 = irlba(clime.sm.cov, nv=5)

basis_coef <- function(vec, basis) {
  vec = vec - mean(vec)
  nbase = ncol(basis)
  coef = rep(0, nbase)
  for (i in 1:nbase) {
    coef[i] = vec %*% basis[,i]
  }
  return(coef)
}

basis_coef(clime[1,], cl.eig$u)
basis_coef(clime.sm[1,], clsm.eig3$u)


# look about the same
plot(clsm.eig0$u[,1])
plot(clsm.eig1$u[,1])
plot(clsm.eig2$u[,1])
plot(clsm.eig3$u[,1])


