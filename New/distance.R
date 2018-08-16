rm(list = ls()); gc()

library(extdepth)
library(tictoc)
library(fdasrvf)

# sims
gp1d = function(fields = 100, mu = 0, l = 50, pts = 50) {
  
  grid = 1:pts
  distmat = as.matrix(dist(grid))
  
  # calc sigma with cov kernel
  sigma = exp(-distmat / l)
  
  sigma.eig = eigen(sigma)
  sigma.half = sigma.eig$vectors %*% diag(sqrt(sigma.eig$values)) %*% t(sigma.eig$vectors)
  
  gps = matrix(0, pts, fields)
  for(f in 1:fields) {
    gps[,f] = (sigma.half %*% rnorm(pts)) + mu
  }
  return(gps)
}
plt_funs = function(f, g) {
  f = as.matrix(f)
  g = as.matrix(g)
  plot(f[,1], type = "l", ylim = c(min(cbind(f, g)), max(cbind(f, g))), col = "red")
  for(i in 2:ncol(f)) {
    lines(f[,i], col = "red")
  }
  for(i in 1:ncol(g)) {
    lines(g[,i], col = "blue")
  }
}

# distance
num.int = function(y, x = seq_len(length(y)) / length(y)) {
  sum(diff(x) * (head(y, -1) + tail(y, -1))) / 2
}
l2 = function(f, g) {
  sqrt(num.int((f-g)^2))
}
ks.dist = function(x, y) {
  nx = ncol(x)
  ny = ncol(y)
  
  x_dist = rep(0, nx)
  y_dist = rep(0, ny)
  
  x_ed = edepth_set(x)
  x_med = x[,which.max(x_ed)]
  
  for(f in 1:nx) {
    x_dist[f] = l2(x[,f], x_med)
  }
  for(f in 1:ny) {
    y_dist[f] = l2(y[,f], x_med)
  }
  
  x_dist = x_dist[x_dist > 0]
  
  ks.test(x_dist, y_dist)$p.value
  # wilcox.test(x_dist, y_dist)$p.value
  # t.test(x_dist, y_dist)$p.value
}

# depth
contained = function(l, u, x) {
  min((l < x) & (u > x))
}

# cosine examples
cos_mu = cos(seq(-pi, pi, length.out = 50))
plot(cos_mu)

gp1 = gp1d()
gp2 = gp1d()
plt_funs(gp1, gp2)

gp1.cr = central_region(gp1, edepth_set(gp1))
plt_funs(cbind(gp1.cr$lower, gp1.cr$upper), gp2)

mean(sapply(1:ncol(gp2), function(x) contained(gp1.cr$lower, gp1.cr$upper, gp2[,x])))

gp2.cr = central_region(gp2, edepth_set(gp2))
plt_funs(cbind(gp2.cr$lower, gp2.cr$upper), gp1)

mean(sapply(1:ncol(gp1), function(x) contained(gp2.cr$lower, gp2.cr$upper, gp1[,x])))

ks.dist(gp1, gp2)

time = 1:50/50
qp1 = f_to_srvf(gp1, time)
qp2 = f_to_srvf(gp2, time)
plt_funs(qp1, qp2)

ks.dist(qp1, qp2)

# trend example
up_mu = seq(-1, 1, length.out = 50) / 2
dn_mu = seq(1, -1, length.out = 50) / 2

gp1 = gp1d(mu = up_mu, l = 50)
gp2 = gp1d(mu = dn_mu, l = 50)
plt_funs(gp1, gp2)

ks.dist(gp2, gp1)


# transformation idea
gp1.med = gp1[,which.max(edepth_set(gp1))]
gp2.med = gp2[,which.max(edepth_set(gp2))]

qp1 = gp1 - gp2.med
qp2 = gp2 - gp1.med
plt_funs(qp1, qp2)

ks.dist(qp2, qp1)


depth_test = function(f, g) {
  fmed = f[,which.max(edepth_set(f))]
  gmed = g[,which.max(edepth_set(g))]
  
  f = f - gmed
  g = g - fmed
  
  ks.dist(f, g)
}

# type 1 error is going to be way too high
f = gp1d(l = 100)
g = gp1d(l = 100) + 0
plt_funs(f, g)
ks.dist(f, g)

fmed = f[,which.max(edepth_set(f))]
gmed = g[,which.max(edepth_set(g))]

f[,1] %*% fmed
num.int(f[,1] * fmed)

f = f - gmed
g = g - fmed
plt_funs(f, g)
ks.dist(f, g)

