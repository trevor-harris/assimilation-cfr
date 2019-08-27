<<<<<<< HEAD
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
gp2d = function(fields = 100, mu = 0, sd = 1, l = 30, pts = 30) {
  grid = 1:pts
  grid = expand.grid(grid, grid)
  distmat = as.matrix(dist(grid))
  
  # calc sigma with cov kernel
  sigma = exp(-distmat / l)
  
  sigma.eig = eigen(sigma)
  sigma.half = sigma.eig$vectors %*% diag(sqrt(sigma.eig$values)) %*% t(sigma.eig$vectors)
  
  gps = array(0, dim=c(pts, pts, fields))
  for(f in 1:fields) {
    gps[,,f] = (sigma.half %*% rnorm(pts^2, sd = sd)) + mu
  }
  return(gps)
}
flatten = function(mat) {
  matrix(mat, prod(dim(mat)[1:2]), dim(mat)[3])
}

### matern based gp
gp1d.mat = function(fields = 100, mu = 0, sd = 1, l = 1, pts = 50, 
                    range = 1, alpha = 1/range, smoothness = 0.5, nu = smoothness, 
                    phi = 1) {
  
  grid = seq(0, 1, length.out = pts)
  distmat = as.matrix(dist(grid))
  
  sigma = fields::Matern(distmat, range = range, alpha = alpha, 
                         smoothness = smoothness, nu = nu, 
                         phi = phi)
  
  sigma.eig = eigen(sigma)
  sigma.half = sigma.eig$vectors %*% diag(sqrt(sigma.eig$values)) %*% t(sigma.eig$vectors)
  
  gps = matrix(0, pts, fields)
  for(f in 1:fields) {
    gps[,f] = (sigma.half %*% rnorm(pts, sd = sd)) + mu
  }
  return(gps)
}

gp2d.mat = function(fields = 100, mu = 0, sd = 1, l = 1, pts = 25, 
                    range = 1, alpha = 1/range, smoothness = 0.5, nu = smoothness, 
                    phi = 1) {
  
  grid = seq(0, 1, length.out = pts)
  grid = expand.grid(grid, grid)
  distmat = as.matrix(dist(grid))
  
  sigma = fields::Matern(distmat, range = range, alpha = alpha, 
                         smoothness = smoothness, nu = nu, 
                         phi = phi)
  
  sigma.eig = eigen(sigma)
  sigma.half = sigma.eig$vectors %*% diag(sqrt(sigma.eig$values)) %*% t(sigma.eig$vectors)
  
  gps = array(0, dim=c(pts, pts, fields))
  for(f in 1:fields) {
    gps[,,f] = (sigma.half %*% rnorm(pts^2, sd = sd)) + mu
  }
  return(gps)
}


plt_funs = function(f, g, yl = "Value", fcol = "red", gcol = "blue", domain = seq(0, 1, length.out = nrow(f))) {
  domain = domain
  
  f = melt(f)
  f$Var1 = domain
  
  plt = ggplot() +
    geom_line(data = f, aes(x = Var1, y = value, group = Var2, alpha = 0.5), color = fcol)
  
  if(!missing(g)) {
    g = melt(g)
    g$Var1 = domain
    
    plt = plt + 
      geom_line(data = g, aes(x = Var1, y = value, group = Var2, alpha = 0.5), color = gcol)
  }
  plt = plt + 
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(legend.position="none") +
    ylab(yl) +
    xlab("Time")
  
  print(plt)
}

# plt_funs = function(f, g, main = "Functions") {
#   
#   if(missing(g)) {
#     f = as.matrix(f)
#     plot(f[,1], type = "l", ylim = c(min(f), max(f)), col = "red", main = main)
#     for(i in 2:ncol(f)) {
#       lines(f[,i], col = "red")
#     }
#   }
#   else {
#     f = as.matrix(f)
#     g = as.matrix(g)
#     plot(f[,1], type = "l", ylim = c(min(cbind(f, g)), max(cbind(f, g))), col = "red", main = main)
#     for(i in 2:ncol(f)) {
#       lines(f[,i], col = "red")
#     }
#     for(i in 1:ncol(g)) {
#       lines(g[,i], col = "blue")
#     }
#   }
=======
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
gp2d = function(fields = 100, mu = 0, sd = 1, l = 30, pts = 30) {
  grid = 1:pts
  grid = expand.grid(grid, grid)
  distmat = as.matrix(dist(grid))
  
  # calc sigma with cov kernel
  sigma = exp(-distmat / l)
  
  sigma.eig = eigen(sigma)
  sigma.half = sigma.eig$vectors %*% diag(sqrt(sigma.eig$values)) %*% t(sigma.eig$vectors)
  
  gps = array(0, dim=c(pts, pts, fields))
  for(f in 1:fields) {
    gps[,,f] = (sigma.half %*% rnorm(pts^2, sd = sd)) + mu
  }
  return(gps)
}
flatten = function(mat) {
  matrix(mat, prod(dim(mat)[1:2]), dim(mat)[3])
}

### matern based gp
gp1d.mat = function(fields = 100, mu = 0, sd = 1, l = 1, pts = 50, 
                    range = 1, alpha = 1/range, smoothness = 0.5, nu = smoothness, 
                    phi = 1) {
  
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

gp2d.mat = function(fields = 100, mu = 0, sd = 1, l = 1, pts = 25, 
                        range = 1,nu = 0.5, phi = 1) {
  
  grid = seq(0, 1, length.out = pts)
  grid = expand.grid(grid, grid)
  distmat = as.matrix(dist(grid))
  
  sigma = fields::Matern(distmat, range = range, nu = nu, phi = phi)
  sigma.cho = t(chol(sigma))
  
  gps = array(0, dim=c(pts, pts, fields))
  for(f in 1:fields) {
    gps[,,f] = (sigma.cho %*% rnorm(pts^2, sd = sd)) + mu
  }
  return(gps)
}

plt_funs = function(f, g, yl = "Value", fcol = "red", gcol = "blue", domain = seq(0, 1, length.out = nrow(f))) {
  domain = domain
  
  f = melt(f)
  f$Var1 = domain
  
  plt = ggplot() +
    geom_line(data = f, aes(x = Var1, y = value, group = Var2, alpha = 0.5), color = fcol)
  
  if(!missing(g)) {
    g = melt(g)
    g$Var1 = domain
    
    plt = plt + 
      geom_line(data = g, aes(x = Var1, y = value, group = Var2, alpha = 0.5), color = gcol)
  }
  plt = plt + 
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(legend.position="none") +
    ylab(yl) +
    xlab("Time")
  
  print(plt)
}

# plt_funs = function(f, g, main = "Functions") {
#   
#   if(missing(g)) {
#     f = as.matrix(f)
#     plot(f[,1], type = "l", ylim = c(min(f), max(f)), col = "red", main = main)
#     for(i in 2:ncol(f)) {
#       lines(f[,i], col = "red")
#     }
#   }
#   else {
#     f = as.matrix(f)
#     g = as.matrix(g)
#     plot(f[,1], type = "l", ylim = c(min(cbind(f, g)), max(cbind(f, g))), col = "red", main = main)
#     for(i in 2:ncol(f)) {
#       lines(f[,i], col = "red")
#     }
#     for(i in 1:ncol(g)) {
#       lines(g[,i], col = "blue")
#     }
#   }
>>>>>>> cc1e9d5fb30fc062b83ab08c31f08af177e46b7a
# }