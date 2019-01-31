rm(list = ls())

library(tictoc)

#### SIMULATION
gp1d = function(fields = 100, mu = 0, sd = 1, l = 50, pts = 50) {
  grid = 1:pts
  distmat = as.matrix(dist(grid))
  
  # calc sigma with cov kernel
  sigma = exp(-distmat / l)
  
  sigma.eig = eigen(sigma)
  sigma.half = sigma.eig$vectors %*% diag(sqrt(sigma.eig$values)) %*% t(sigma.eig$vectors)
  
  gps = matrix(0, pts, fields)
  for(f in 1:fields) {
    gps[,f] = (sigma.half %*% rnorm(pts, sd = sd)) + mu
  }
  return(gps)
}

# Expected Depth (XD) pointwise
depth = function(g, fmat) {
  
  # Computes the depth values of a function with respect to a set of functions (fmat)
  fn = ncol(fmat)
  depth = rep(0, length(g))
  
  for (row in 1:nrow(fmat)) {
    diff = abs(sum(sign(g[row] - fmat[row,])))
    depth[row] = 1 - (diff / fn)
  }
  
  return(depth)
}

# expected depth 
xdepth = function(gmat, fmat) {
  apply(gmat, 2, function(x) mean(depth(x, fmat)))
}

# true depth pointwise
gpdepth = function(f) {
  ds = sapply(1:length(f), function(x) pnorm(f[x]))
  mean(1 - abs(1 - 2*ds))
}

# true depth
gdepth = function(f) {
  apply(f, 2, gpdepth)
}


sims = 1000
# ns = c(25, 50, 100, 200, 400, 800, 1200)
ns = c(25, 50, 200, 400, 800, 1400, 2000, 2600)
diff = matrix(0, sims, length(ns))

for(j in 1:length(ns)) {
  tic()
  for(i in 1:sims) {
    f = gp1d(ns[j], pts = 50, l = 10)
    xf = xdepth(f, f)
    gf = gdepth(f)
    
    # diff[i, j] = mean(xf < 0.3) - mean(gf < 0.3)
    diff[i, j] = mean(as.numeric(xf < 0.3) - as.numeric(gf < 0.3))
  }
  toc()
}
boxplot(sapply(1:length(ns), function(x) diff[,x] * sqrt(ns[x])))
abline(h = 0)

scaled_diff = sapply(1:length(ns), function(x) diff[,x] * sqrt(ns[x]))
plot(apply(scaled_diff, 2, sd))


f = gp1d(100, pts = 20)
xf = mean(depth(f[,1], f))
gf = gdepth(f)[1]


sims = 1000
ns = c(25, 50, 100, 200, 400)
diff = matrix(0, sims, length(ns))

for(j in 1:length(ns)) {
  tic()
  for(i in 1:sims) {
    f = gp1d(ns[j], pts = 20)

    diff[i, j] = abs(as.numeric(mean(depth(f[,1], f)) < 0.3) - as.numeric(gdepth(f)[1] < 0.3))
  }
  toc()
}
boxplot(sapply(1:length(ns), function(x) diff[,x] * sqrt(ns[x])))
abline(h = 0)

1-apply(diff, 2, mean)




sims = 1000
ns = c(25, 50, 100, 200, 400)
diff = matrix(0, sims, length(ns))

for(j in 1:length(ns)) {
  tic()
  for(i in 1:sims) {
    f1 = gp1d(ns[j], pts = 20)
    f2 = gp1d(ns[j], pts = 20)
    
    f1d = xdepth(f, f)
    gf = gdepth(f)
    
    # diff[i, j] = mean(xf < 0.3) - mean(gf < 0.3)
    diff[i, j] = mean(as.numeric(xf < 0.3) - as.numeric(gf < 0.3))
  }
  toc()
}
boxplot(sapply(1:length(ns), function(x) diff[,x] * sqrt(ns[x])))
abline(h = 0)





xf[1]
gpdepth(f[,1])

sims = 100
diff = rep(0, sims)
for(i in 1:sims) {
  f = gp1d(100, pts = 20)
  xf = xdepth(f, f)
  
  diff[i] = sqrt() * (xf[1] - gpdepth(f[,1]))
}
boxplot(diff)



# ok so the depths converge to normals with a smaller sd
sims = 1000
ns = c(25, 50, 100, 200, 400, 800)
diff = matrix(0, sims, length(ns))

x1 = gp1d(1, pts = 20)
for(j in 1:length(ns)) {
  for(i in 1:sims) {
    f = cbind(x1, gp1d(ns[j], pts = 20))
    diff[i, j] = sqrt(ns[j]) * (mean(depth(f[,1], f)) - gpdepth(f[,1]))
  }
}
boxplot(diff)
abline(h = 0)

diff = melt(diff)
diff[["Var2"]] = as.factor(diff[["Var2"]])
levels(diff[["Var2"]]) = ns

ggplot(diff, aes(x = Var2, y = value)) + 
  geom_boxplot() + 
  geom_hline(yintercept = 0) +
  theme_classic() +
  xlab("Number of functions (n)") +
  ylab("sqrt(n) * (D(x, P_n) - D(x, P))")


# True Depth CDFs converge to Normal
#####
sims = 100
ns = c(25, 50, 100, 200, 400, 800)
ns = c(1000)
diff = matrix(0, sims, length(ns))

for(j in 1:length(ns)) {
  for(i in 1:sims) {
    f = gp1d(ns[j], pts = 20)
    g = gp1d(ns[j], pts = 20)
    
    gpf = apply(f, 2, gpdepth)
    gpg = apply(g, 2, gpdepth)
    
    diff[i, j] = sqrt(ns[j]/2) * (mean(gpf < 0.5) - mean(gpg < 0.5) )
  }
}
x = seq(-4, 4, length.out = 10000)
hist(diff[,1] * 2, breaks = 50, probability = T)
lines(density(diff[,1] * 2))
lines(x, sapply(x, dnorm))
####




# Difference between true and empirical depth of G wrt F goes to Normal
#####
sims = 10000
# ns = c(25, 50, 100, 200, 400, 800)
ns = c(1000)
diff = matrix(0, sims, length(ns))

for(j in 1:length(ns)) {
  for(i in 1:sims) {
    f = gp1d(ns[j], pts = 20)
    g = gp1d(ns[j], pts = 20)
    
    gpg = apply(g, 2, gpdepth)
    xg = xdepth(g, f)
    
    diff[i, j] = sqrt(ns[j]/2) * (mean(gpg < 0.5) - mean(xg < 0.5))
  }
}
boxplot(diff)

x = seq(-4, 4, length.out = 10000)
hist(diff[,1] * 2, breaks = 50, probability = T)
lines(density(diff[,1] * 2))
lines(x, sapply(x, dnorm))
####



#####
sims = 100
ns = c(25, 50, 100, 200, 400, 800)
ns = c(1000)
diff = matrix(0, sims, length(ns))

for(j in 1:length(ns)) {
  for(i in 1:sims) {
    f = gp1d(ns[j], pts = 20)
    g = gp1d(ns[j]+2, pts = 20)
    
    gpf = apply(f, 2, gpdepth)
    xf = xdepth(f, f)
    
    diff[i, j] = sqrt(ns[j]/2) * (mean(gpf < 0.5) - mean(xf < 0.5))
  }
}
boxplot(diff)

x = seq(-4, 4, length.out = 10000)
hist(diff[,1] * 4, breaks = 50, probability = T)
lines(density(diff[,1] * 4))
lines(x, sapply(x, dnorm))
####


diff1 = (diff[,1] / sqrt(ns[j]/2)) * (sqrt(ns[j]))
hist(diff1 * 4, breaks = 50, probability = T)
lines(density(diff1 * 4))
lines(x, sapply(x, dnorm))






# ok so the depths converge to normals with a smaller sd
sims = 1000
ns = c(25, 50, 100, 200, 400, 800)
diff = matrix(0, sims, length(ns))

for(j in 1:length(ns)) {
  for(i in 1:sims) {
    f = gp1d(ns[j], pts = 20)
    g = gp1d(ns[j], pts = 20)
    diff[i, j] = sqrt(ns[j]) * (mean(depth(g[,1], f)) - gpdepth(g[,1]))
  }
}
boxplot(diff)
qqnorm(diff[,6] * 2)



# so do the difference of the cdfs
sims = 1000
ns = c(25, 50, 100, 200, 400, 800)
diff = matrix(0, sims, length(ns))

for(j in 1:length(ns)) {
  for(i in 1:sims) {
    f = gp1d(ns[j], pts = 20)
    
    xf = xdepth(f, f)
    xt = apply(f, 2, gpdepth)
    
    diff[i, j] = sqrt(ns[j]) * (mean(xf < 0.5) - mean(xt < 0.5))
  }
}
boxplot(diff)
qqnorm(diff[,6] * 2)




# so do the difference of the cdfs
sims = 1000
ns = c(25, 50, 100, 200, 400, 800, 1600, 3200)
diff = matrix(0, sims, length(ns))

for(j in 1:length(ns)) {
  for(i in 1:sims) {
    f = gp1d(ns[j], pts = 20)
    g = gp1d(ns[j], pts = 20)
    
    xf = xdepth(f, f)
    xg = xdepth(g, f)
    
    gpf = apply(f, 2, gpdepth)
    gpg = apply(g, 2, gpdepth)
    
    diff[i, j] = sqrt(ns[j]/2) * (((mean(xf < 0.5) - mean(gpf < 0.5))) - ((mean(xg < 0.5) - mean(gpg < 0.5))))
  }
}
boxplot(diff)
# qqnorm(diff[,1])

diff = melt(diff)
diff[["Var2"]] = as.factor(diff[["Var2"]])
levels(diff[["Var2"]]) = ns

ggplot(diff, aes(x = Var2, y = value)) + 
  geom_boxplot() + 
  geom_hline(yintercept = 0) +
  theme_classic() +
  xlab("Number of functions (n)") +
  ylab("four part diff")


##### quadratic convergence hopefully
# so do the difference of the cdfs
sims = 200
ns = c(25, 50, 100, 200, 400, 800, 1600)
diff2 = matrix(0, sims, length(ns))

for(j in 1:length(ns)) {
  tic(paste0("N = ", ns[j]))
  for(i in 1:sims) {
    f = gp1d(ns[j], pts = 20)
    g = gp1d(ns[j], pts = 20)
    
    xf = xdepth(f, f)
    xg = xdepth(g, f)
    
    gpf = apply(f, 2, gpdepth)
    gpg = apply(g, 2, gpdepth)
    
    f_hat = mean(xf < 0.5)
    g_hat = mean(xg < 0.5)
    
    f_til = mean(gpf < 0.5)
    g_til = mean(gpg  < 0.5)
    
    diff2[i, j] = (ns[j]/2) * ((f_hat - f_til)^2 - (f_hat - f_til)*(g_hat - g_til))
  }
  toc()
}
boxplot(diff2)

plot(apply(diff2, 2, mean))




#### just need the covariance to decay fast enough
# which is the same as the probabilities decaying in this case

f = gp1d(100, pts = 20)
xf = rev(sort(xdepth(f, f)))
# plot(xf)

t = seq(0, 1, length.out = 10000)
mean(xf[1] <= t) * (mean(xf[300] > t))

covtotal = 0
for(i in 1:(ncol(f)-1)) {
  for(j in (i+1):ncol(f)) {
    covtotal = covtotal + (mean(xf[i] <= t) * (mean(xf[j] > t)))
  }
}
covtotal / ncol(f)


t = 0.4
sims = 200
ns = 50
cov12 = 0
for(k1 in 1:(ns-1)) {
  tic(paste0("k1 = ", k1))
  
  for(k2 in k1:ns) {
    
    p1 = rep(0, sims)
    p2 = rep(0, sims)
    
    for(i in 1:sims) {
      f = gp1d(ns, pts = 20)
      xf = rev(sort(xdepth(f, f)))
      
      p1[i] = as.integer(xf[k1] < t)
      p2[i] = as.integer(xf[k2] < t)
      
    }
    cov12 = cov12 + cov(p1, p2)
  }
  toc()
}

cov12 / ns


##### what about just the variance of the sums
t = 0.5
sims = 500
ns = 500

j = 1
v1 = rep(0, sims)
v2 = rep(0, sims)
for(i in 1:sims) {
  f = gp1d(ns[j], pts = 20)
  g = gp1d(ns[j], pts = 20)
  
  xf = xdepth(f, f)
  xg = xdepth(g, f)
  
  gpf = apply(f, 2, gpdepth)
  gpg = apply(g, 2, gpdepth)
  
  v1[i] = (sum(xf < t) - sum(gpf < t)) / ns[j]
  v2[i] = (sum(xg < t) - sum(gpg < t)) / ns[j]
}

# diff of variances
ns[j]*(var(v1) - var(v2))

# OG
ns[j]*(mean(v1^2) + mean(v2^2) - 2*mean(v1 * v2))




##### small test
##### what about just the variance of the sums
t = 0.3
sims = 500
ns = as.integer(seq(25, 200, length.out = 30))

varest = rep(0, length(ns))
cvest1 = rep(0, length(ns))
cvest2 = rep(0, length(ns))

v1 = rep(0, sims)
v2 = rep(0, sims)
for(j in 1:length(ns)) {
  tic(paste0("N = ", ns[j]))
  
  for(i in 1:sims) {
    f = gp1d(ns[j], pts = 20)
    g = gp1d(ns[j], pts = 20)
    
    xf = xdepth(f, f)
    xg = xdepth(g, f)
    
    gpf = apply(f, 2, gpdepth)
    gpg = apply(g, 2, gpdepth)
    
    v1[i] = (sum(xf < t) - sum(gpf < t)) / ns[j]
    v2[i] = (sum(xg < t) - sum(gpg < t)) / ns[j]
  }
  varest[j] = ns[j]*(var(v1) - var(v2))
  cvest1[j] = ns[j]*(var(v1) - mean(v1*v2))
  cvest2[j] = ns[j]*(var(v2) - mean(v1*v2))
  
  toc()
}

plot(abs(varest))
plot(abs(cvest1))
plot(abs(cvest2))



#### Just the variances
##### small test
##### what about just the variance of the sums
t = 0.7
sims = 1000
ns = as.integer(seq(10, 1000, length.out = 50))

varest = rep(0, length(ns))
cvest1 = rep(0, length(ns))
# cvest2 = rep(0, length(ns))

v1 = rep(0, sims)
v2 = rep(0, sims)
tic("Total")
for(j in 1:length(ns)) {
  tic(paste0("N = ", ns[j]))
  
  for(i in 1:sims) {
    f = gp1d(ns[j], pts = 20)
    # g = gp1d(ns[j], pts = 20)
    
    xf = xdepth(f, f)
    # xg = xdepth(g, f)
    
    gpf = apply(f, 2, gpdepth)
    # gpg = apply(g, 2, gpdepth)
    
    v1[i] = mean(xf < t) - mean(gpf < t)
    # v2[i] = mean(xg < t) - mean(gpg < t)
  }
  
  cvest1[j] = ns[j]*var(v1) / 2
  # cvest2[j] = ns[j]*var(v2) / 2
  
  toc()
}
toc()

plot(abs(cvest1), main = "l2 convergence of fn - fn")
# plot(abs(cvest2), main = "l2 convergence of gn - gn")




##### what about just the variance of the sums
t = 0.4
sims = 1000
ns = as.integer(seq(25, 700, length.out = 10))

varest = rep(0, length(ns))
cvest1 = rep(0, length(ns))
cvest2 = rep(0, length(ns))

v1 = rep(0, sims)
tic("Total")
for(j in 1:length(ns)) {
  tic(paste0("N = ", ns[j]))
  
  for(i in 1:sims) {
    f = gp1d(ns[j], pts = 20)
    
    xf = xdepth(f, f)
    gpf = apply(f, 2, gpdepth)
    
    # a = as.integer(xf < t)
    # b = as.integer(gpf < t)
    
    v1[i] = mean(as.integer(xf < t) - as.integer(gpf < t))
  }
  
  cvest1[j] = ns[j]*var(v1)
  
  toc()
}
toc()

plot(abs(cvest1), main = "l2 convergence of 1 - fn")



# # so do the difference of the cdfs
# sims = 500
# ns = c(25, 100, 400, 1600)
# diff = matrix(0, sims, length(ns))
# 
# for(j in 1:length(ns)) {
#   for(i in 1:sims) {
#     f = gp1d(ns[j], pts = 20)
#     g = gp1d(ns[j], pts = 20)
#     
#     xf = xdepth(f, f)
#     xg = xdepth(g, f)
#     
#     f_til = mean(gpf < 0.5)
#     g_til = mean(gpg  < 0.5)
#     
#     gpf = apply(f, 2, gpdepth)
#     gpg = apply(g, 2, gpdepth)
#     
#     diff[i, j] = sqrt(ns[j]^2 / (2*(ns[j]))) * (((mean(xf < 0.5))) - ((mean(xg < 0.5))))
#   }
# }
# boxplot(diff)
# # qqnorm(diff[,1])
# 
# diff = melt(diff)
# diff[["Var2"]] = as.factor(diff[["Var2"]])
# levels(diff[["Var2"]]) = ns
# 
# ggplot(diff, aes(x = Var2, y = value)) + 
#   geom_boxplot() + 
#   geom_hline(yintercept = 0) +
#   theme_classic() +
#   xlab("Number of functions (n)") +
#   ylab("four part diff")
# 
# 
# 
# 
# 
# # check them correlations
# sims = 1000
# ns = 400
# cdfs = matrix(0, sims, 2)
# 
# j = 1
# for(i in 1:sims) {
#   f = gp1d(ns[j], pts = 20)
#   g = gp1d(ns[j], pts = 20)
#   
#   xf = xdepth(f, f)
#   xg = xdepth(g, f)
#   
#   gpf = apply(f, 2, gpdepth)
#   gpg = apply(g, 2, gpdepth)
#   
#   cdfs[i, 1] = mean(xf < 0.5) - mean(gpf < 0.5)
#   cdfs[i, 2] = mean(xg < 0.5) - mean(gpg < 0.5)
# }
# 
# # plot(cdfs)
# # plot((cdfs[,1] - cdfs[,2]))
# 
# max((cdfs[,1] - cdfs[,2]))
# 1 / sqrt(ns/2)
# 
# mean((cdfs[,1] - cdfs[,2]) < 1 / sqrt(ns/2))
# 
>>>>>>> 559e76f7f14af028c340bd9fa130d8a12858685f
