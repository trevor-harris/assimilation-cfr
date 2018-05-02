# helper function
isbetween = function(x, a, b) {
  if(x < a) return(0)
  if(x > b) return(0)
  else return(1)
}

#### SIMULATION ####
sim_gp = function(fields = 100, mu = 0, l = 1, pts = 30, ker = "exp") {
  # kernel function
  exp_kernel <- function(X,l){
    D <- plgp::distance(X)
    Sigma <- exp(-D/l)
  }
  
  gaus_kernel <- function(X, l) {
    D <- plgp::distance(X)
    Sigma <- exp(-(D/l)^2)
  }
  
  # spatial locations to estimate
  grid = seq(-1, 1, length = pts)
  grid = expand.grid(grid, grid)
  
  # calc sigma with cov kernel
  if (ker == "exp") sigma = exp_kernel(grid, l = l)
  if (ker == "gauss") sigma = gaus_kernel(grid, l = l)
  # calc mu
  if(length(mu) == 1) {
    mu = rep(mu, dim(sigma)[1])
  } else {
    mu = mu
  }
  
  # sample fields from the GP
  gps = MASS::mvrnorm(fields, mu, sigma)
  
  if(fields > 1) {
    gps = vapply(1:fields, function(i) matrix(gps[i,], pts, pts) 
                 ,FUN.VALUE = array(0, dim = c(pts, pts)))
  } else {
    gps = array(gps, dim=c(pts, pts, 1))
  }
  return(gps)
}



#### SIMULATION ####
sim_gp.1d = function(fields = 100, mu = 0, l = 1, pts = 30) {
  # kernel function
  exp_kernel <- function(X,l){
    D <- plgp::distance(X)
    Sigma <- exp(-D/l)
  }
  
  # spatial locations to estimate
  grid = seq(-1, 1, length = pts)
  
  # calc sigma with cov kernel
  sigma = exp_kernel(grid, l = l)
  
  # calc mu
  if(length(mu) == 1) {
    mu = rep(mu, dim(sigma)[1])
  } else {
    mu = mu
  }
  
  # sample fields from the GP
  gps = MASS::mvrnorm(fields, mu, sigma)
  return(gps)
}

