library(MASS)
library(fields)
library(mvtnorm)

#### SIMULATION
# each of these functions will generate realizations from the named process
# fields: number of fields (functions) to generate
# pts: number of sample points in the field / function
# mu: mean parameter, Can be a scalar or a function the size of the field
# sd: standard deviation parameter in the matern covariance
# range: range parameter in the matern covariance
# nu: smoothness parameter from the matern covariance function


## 1d gaussian process
gp1d = function(fields = 100, mu = 0, sd = 1, pts = 50, range = 0.3, nu = 1) {
  
  grid = seq(0, 1, length.out = pts)
  distmat = as.matrix(dist(grid))
  
  sigma = fields::Matern(distmat, range = range, nu = nu)
  sigma.cho = t(chol(sigma))
  
  gps = matrix(0, pts, fields)
  for(f in 1:fields) {
    gps[,f] = (sigma.cho %*% rnorm(pts, sd = sd)) + mu
  }
  return(gps)
}

## 2d gaussian process
gp2d = function(fields = 100, mu = 0, sd = 1, l = 1, pts = 25, 
                    range = 1, nu = 0.5) {
  
  grid = seq(0, 1, length.out = pts)
  grid = expand.grid(grid, grid)
  distmat = as.matrix(dist(grid))
  
  sigma = fields::Matern(distmat, range = range, nu = nu)
  sigma.cho = t(chol(sigma))
  
  gps = array(0, dim=c(pts, pts, fields))
  for(f in 1:fields) {
    gps[,,f] = (sigma.cho %*% rnorm(pts^2, sd = sd)) + mu
  }
  return(gps)
}

flatten = function(mat) {
  matrix(mat, prod(dim(mat)[1:2]), dim(mat)[3])
}

## 1d t-process
tp1d = function(fields = 10, mu = 0, df = 1, sd = 1, l = 1, pts = 50, 
                range = 0.1, nu = 1) {
  
  grid = seq(0, 1, length.out = pts)
  distmat = as.matrix(dist(grid))
  
  sigma = fields::Matern(distmat, range = range, nu = nu)
  
  t(rmvt(fields, sigma, df)) + mu
}

## 2d t-process
tp2d = function(fields = 100, mu = 0, df = 2, sd = 1, l = 1, pts = 25, 
                range = 0.1, nu = 1) {
  
  grid = seq(0, 1, length.out = pts)
  grid = expand.grid(grid, grid)
  distmat = as.matrix(dist(grid))
  
  sigma = fields::Matern(distmat, range = range, nu = nu) / 3
  sigma.cho = t(chol(sigma))
  
  gps = array(0, dim=c(pts, pts, fields))
  for(f in 1:fields) {
    chi = sqrt(df / rchisq(1, df = df))
    td = as.vector(sigma.cho %*% rnorm(pts^2, sd = sd))
    
    gps[,,f] =  td * chi + mu
  }
  
  return(gps)
}