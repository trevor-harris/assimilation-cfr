library(refund)
library(roahd)
library(fChange)

source("../research/proxy/code/depths.R")
source("../research/proxy/code/depth_tests.R")
source("../research/proxy/code/simulation.R")


gp1d = function(fields = 20, mu = 0, sd = 1, l = 1, pts = 50, 
                range = 0.4, nu = 1.5, phi = 1) {
  
  grid = seq(0, 1, length.out = pts)
  distmat = as.matrix(dist(grid))
  
  sigma = fields::Matern(distmat, range = range, nu = nu, phi = phi)
  sigma.cho = t(chol(sigma))
  
  gps = matrix(0, pts, fields)
  for(f in 1:fields) {
    gps[,f] = (sigma.cho %*% rnorm(pts, sd = sd)) + mu
  }
  return(gps)
}
approxFourier = function(basis, pts = 50) {
  t = seq(0, 2*pi, length.out = 50)
  fmat = matrix(0, pts, ncol(basis))
  
  for(f in 1:ncol(fmat)) {
    
    fmat[,f] = fmat[,f] + basis[1, f]
    
    for(i in 1:((nrow(basis)-1)/2)) {
      fmat[,f] = fmat[,f] + basis[2*i, f] * sin(i*t)
    }
    for(i in 1:((nrow(basis)-1)/2)) {
      fmat[,f] = fmat[,f] + basis[2*i + 1, f] * cos(i*t)
    }
  }
  fmat
}

p = 50
f = gp1d(50, sd = 1, pts = p)
g = gp1d(50, sd = 1.5, pts = p)

plt_funs(f, g)

fadtest(f, g, p/2)
bandtest(f, g)[2]
kolm(f, g)[2]



nb = 50
fcoef = vapply(1:nb, function(x) rnorm(100, sd = 1/x^(1/3)), FUN.VALUE = rep(0, 100))
gcoef = vapply(1:nb, function(x) rnorm(100, sd = 1/x^(1/3)), FUN.VALUE = rep(0, 100))

fcoef = apply(fcoef, 2, function(x) x - mean(x))
gcoef = apply(gcoef, 2, function(x) x - mean(x))

boxplot(fcoef)

gcoef = sapply(1:nb, function(x) gcoef[,x] + rnorm(1, sd = 0.3) * (x > nb/2) * sample(c(0, 1), 1))

boxplot(gcoef)

f = approxFourier(t(fcoef))
g = approxFourier(t(gcoef))

fadtest(f, g)
bandtest(f, g)[2]
kolm(f, g)[2]

nb = 120
sims = 1000
out = matrix(0, nrow = sims, ncol = 3)
for(i in 1:sims) {
  fcoef = vapply(1:nb, function(x) rnorm(100, sd = 1/sqrt(x)), FUN.VALUE = rep(0, 100))
  gcoef = vapply(1:nb, function(x) rnorm(100, sd = 1/sqrt(x)), FUN.VALUE = rep(0, 100))
  
  fcoef = apply(fcoef, 2, function(x) x - mean(x))
  gcoef = apply(gcoef, 2, function(x) x - mean(x))
  
  gcoef = sapply(1:nb, function(x) gcoef[,x] + rnorm(1, sd = 0.15) * (x > 10) * sample(c(0, 1), 1))
  # gcoef = sapply(1:nb, function(x) gcoef[,x] + rnorm(1, sd = 0.45) * as.integer((1:nb %% 5) == 0))
  
  f = approxFourier(t(fcoef))
  g = approxFourier(t(gcoef))
  
  out[i, ] = c(fadtest(f, g)[2], bandtest(f, g)[2], kolm(f, g)[2])
  
}
colMeans(out < 0.05)

boxplot(fcoef)
boxplot(gcoef)
plt_funs(f, g)


