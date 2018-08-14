rm(list = ls())

library(extdepth)
library(tictoc)
library(fdasrvf)

#### FUNCTIONS ####
ks.dist = function(x, y) {
  nx = ncol(x)
  ny = ncol(y)
  
  x_dist = rep(0, nx)
  y_dist = rep(0, ny)
  
  x_ed = edepth_set(x)
  x_med = x[,which.max(x_ed)]
  
  for(f in 1:nx) {
    x_dist[f] = l2(x[,f], x_med)
    y_dist[f] = l2(y[,f], x_med)
  }
  
  x_dist = x_dist[x_dist > 0]
  
  ks.test(x_dist, y_dist)$p.value
}
num.int = function(y, x = seq_len(length(y)) / length(y)) {
  sum(diff(x) * (head(y, -1) + tail(y, -1))) / 2
}
l2 = function(f, g) {
  sqrt(num.int((f-g)^2))
}


#### SIMULATION ####
t = seq(-1, 1, length.out = 100)
f1 = t^2
f1 = outer(f1, seq(-2, 2, length.out = 10), FUN = "+") + rnorm(10, 0, 0.1)
f1 = f1 - mean(f1)

f2 = 1 - t^2
f2 = outer(f2, seq(-2, 2, length.out = 10), FUN = "+") + rnorm(10, 0, 0.1)
f2 = f2 - mean(f2)

f = cbind(f1, f2)




#### SRSF ####
q = f_to_srvf(f, t)

plot(f[,1], type = "l", ylim = c(-3, 3))
for(i in 2:20) {
  lines(f[,i])
}

plot(q[,1], type = "l", ylim = c(-3, 3))
for(i in 2:20) {
  lines(q[,i])
}

q1 = q[,1:10]
q2 = q[,11:20]

ks.dist(f1, f2)
ks.dist(q1, q2)
