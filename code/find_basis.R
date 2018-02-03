library(ncdf4)
library(dplyr)
library(plotly)
library(extdepth)
library(mgcv)
library(irlba)


### AFTER PRIOR HAS BEEN PROCESSED

# convert to big matrix
prior.mat = matrix(0, length(prior), length(prior[[1]]))
for (t in 1:length(prior)) {
  cat(t, "\n")
  prior.mat[t, ] = as.vector(unlist(prior[[t]]))
}

prior.cov = cov(prior.mat)
prior.eig = irlba(prior.cov, nv=20)

basis_coef <- function(vec, basis) {
  vec = vec - mean(vec)
  nbase = ncol(basis)
  coef = rep(0, nbase)
  for (i in 1:nbase) {
    coef[i] = vec %*% basis[,i]
  }
  return(coef)
}

basis_coef(as.vector(ens[[1]]), prior.eig$u)
prior.eig$d



#### TESTING SUBSAMPLING

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


# look about the same
plot(clsm.eig0$u[,1])
plot(clsm.eig1$u[,1])
plot(clsm.eig2$u[,1])
plot(clsm.eig3$u[,1])


