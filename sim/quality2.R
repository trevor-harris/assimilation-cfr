# T TEST STUFF

rm(list = ls()); gc()
library(extdepth)
library(tictoc)

# sims
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
plt_funs = function(f, g, main = "Plot") {
  
  if(missing(g)) {
    f = as.matrix(f)
    plot(f[,1], type = "l", ylim = c(min(f), max(f)), col = "red", main = main)
    for(i in 2:ncol(f)) {
      lines(f[,i], col = "red")
    }
  }
  else {
    f = as.matrix(f)
    g = as.matrix(g)
    plot(f[,1], type = "l", ylim = c(min(cbind(f, g)), max(cbind(f, g))), col = "red", main = main)
    for(i in 2:ncol(f)) {
      lines(f[,i], col = "red")
    }
    for(i in 1:ncol(g)) {
      lines(g[,i], col = "blue")
    }
  }
}

# depth CDF
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
xdepth = function(gmat, fmat) {
  apply(gmat, 2, function(x) mean(depth(x, fmat)))
}

# SYM KS
ks_cdf = function(x, n = 10) {
  1 - 2*(sum(sapply(1:n, function(k) ((-1)^(k-1)) * exp(-2*(k^2)*(x^2)))))
}

coverage = function(f, g) {
  ffxd = xdepth(f, f)
  gfxd = xdepth(g, f)
  fgxd = xdepth(f, g)
  ggxd = xdepth(g, g)
  
  tf = seq(0, 1, length.out = max(1000, 3*length(ffxd)))  
  ffr = sapply(tf, function(y) mean(ffxd <= y))
  gfr = sapply(tf, function(y) mean(gfxd <= y))
  
  tg = seq(0, 1, length.out = max(1000, 3*length(ggxd)))
  fgr = sapply(tg, function(y) mean(fgxd <= y))
  ggr = sapply(tg, function(y) mean(ggxd <= y))
  
  rate = sqrt((ncol(g)*ncol(f)) / (ncol(g) + ncol(f)))
  
  ksf = rate*max(abs(ffr - gfr))
  ksg = rate*max(abs(fgr - ggr))
  
  1 - ks_cdf(max(ksf, ksg))^2
}
quality = function(f, g) {
  fxd = xdepth(f, f)
  gxd = xdepth(g, f)
  
  r = sapply(gxd, function(y) mean(y >= fxd))
  1 - pnorm((mean(r) - 0.5) / sqrt((1/ncol(f) + 1/ncol(g))/12))
}


set.seed(1023)
sims = 1000
regina = rep(0, sims)
trevor = rep(0, sims)
for(s in 1:sims) {
  tic("Total")
  cat("Simulation ", s, "\n")
  
  gp1 = gp1d(500, l = 20)
  gp2 = gp1d(500, l = 20)
  
  regina[s] = quality(gp1, gp2)
  trevor[s] = coverage(gp1, gp2)
  
  toc()
  cat("Regina: ", mean(regina[1:s] < 0.05), "\n")
  cat("Trevor: ", mean(trevor[1:s] < 0.05), "\n")
  cat("\n")
}

regina.boot20 = sapply(1:100, function(x) mean(sample(regina, size = length(regina), TRUE) <= 0.05))
trevor.boot20 = sapply(1:100, function(x) mean(sample(trevor, size = length(trevor), TRUE) <= 0.05))

boxplot(regina.boot, trevor.boot)

N = 250

set.seed(5)
sims = 1000
regina = rep(0, sims)
trevor = rep(0, sims)
for(s in 1:sims) {
  tic("Total")
  cat("Simulation ", s, "\n")
  
  gp1 = gp1d(N, l = 5)
  gp2 = gp1d(N, l = 5)
  
  regina[s] = quality(gp1, gp2)
  trevor[s] = coverage(gp1, gp2)
  
  toc()
  cat("Regina: ", mean(regina[1:s] < 0.05), "\n")
  cat("Trevor: ", mean(trevor[1:s] < 0.05), "\n")
  cat("\n")
}

regina.boot5 = sapply(1:100, function(x) mean(sample(regina, size = length(regina), TRUE) <= 0.05))
trevor.boot5 = sapply(1:100, function(x) mean(sample(trevor, size = length(trevor), TRUE) <= 0.05))

# boxplot(regina.boot5, trevor.boot5, main = "Smoothness = 5")



set.seed(10)
sims = 1000
regina = rep(0, sims)
trevor = rep(0, sims)
for(s in 1:sims) {
  tic("Total")
  cat("Simulation ", s, "\n")
  
  gp1 = gp1d(N, l = 10)
  gp2 = gp1d(N, l = 10)
  
  regina[s] = quality(gp1, gp2)
  trevor[s] = coverage(gp1, gp2)
  
  toc()
  cat("Regina: ", mean(regina[1:s] < 0.05), "\n")
  cat("Trevor: ", mean(trevor[1:s] < 0.05), "\n")
  cat("\n")
}

regina.boot10 = sapply(1:100, function(x) mean(sample(regina, size = length(regina), TRUE) <= 0.05))
trevor.boot10 = sapply(1:100, function(x) mean(sample(trevor, size = length(trevor), TRUE) <= 0.05))

# boxplot(regina.boot10, trevor.boot10, main = "Smoothness = 10")



set.seed(20)
sims = 1000
regina = rep(0, sims)
trevor = rep(0, sims)
for(s in 1:sims) {
  tic("Total")
  cat("Simulation ", s, "\n")
  
  gp1 = gp1d(N, l = 20)
  gp2 = gp1d(N, l = 20)
  
  regina[s] = quality(gp1, gp2)
  trevor[s] = coverage(gp1, gp2)
  
  toc()
  cat("Regina: ", mean(regina[1:s] < 0.05), "\n")
  cat("Trevor: ", mean(trevor[1:s] < 0.05), "\n")
  cat("\n")
}

regina.boot20 = sapply(1:100, function(x) mean(sample(regina, size = length(regina), TRUE) <= 0.05))
trevor.boot20 = sapply(1:100, function(x) mean(sample(trevor, size = length(trevor), TRUE) <= 0.05))

# boxplot(regina.boot20, trevor.boot20, main = "Smoothness = 20")



boots = data.frame(Method = c(rep("Quality Index", 300), rep("KSD", 300)),
                   Smoothness = c(rep("5", 100), rep("10", 100), rep("20", 100)),
                   Size = c(regina.boot5, regina.boot10, regina.boot20,
                            trevor.boot5, trevor.boot10, trevor.boot20))

boots$Smoothness = factor(boots$Smoothness, levels = c("5", "10", "20"), ordered = TRUE)

ggplot(boots, aes(x = Smoothness, y = Size, color = Method)) +
  geom_boxplot() +
  geom_abline(intercept = 0.05, slope = 0) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5)) +
  ggtitle("Size by Method and Smoothness")
                   








N = 75
M = 250

set.seed(5)
sims = 1000
regina = rep(0, sims)
trevor = rep(0, sims)
for(s in 1:sims) {
  tic("Total")
  cat("Simulation ", s, "\n")
  
  gp1 = gp1d(N, l = 5)
  gp2 = gp1d(M, l = 5)
  
  regina[s] = quality(gp1, gp2)
  trevor[s] = coverage(gp1, gp2)
  
  toc()
  cat("Regina: ", mean(regina[1:s] < 0.05), "\n")
  cat("Trevor: ", mean(trevor[1:s] < 0.05), "\n")
  cat("\n")
}

regina.boot5 = sapply(1:100, function(x) mean(sample(regina, size = length(regina), TRUE) <= 0.05))
trevor.boot5 = sapply(1:100, function(x) mean(sample(trevor, size = length(trevor), TRUE) <= 0.05))

# boxplot(regina.boot5, trevor.boot5, main = "Smoothness = 5")



set.seed(10)
sims = 1000
regina = rep(0, sims)
trevor = rep(0, sims)
for(s in 1:sims) {
  tic("Total")
  cat("Simulation ", s, "\n")
  
  gp1 = gp1d(N, l = 10)
  gp2 = gp1d(M, l = 10)
  
  regina[s] = quality(gp1, gp2)
  trevor[s] = coverage(gp1, gp2)
  
  toc()
  cat("Regina: ", mean(regina[1:s] < 0.05), "\n")
  cat("Trevor: ", mean(trevor[1:s] < 0.05), "\n")
  cat("\n")
}

regina.boot10 = sapply(1:100, function(x) mean(sample(regina, size = length(regina), TRUE) <= 0.05))
trevor.boot10 = sapply(1:100, function(x) mean(sample(trevor, size = length(trevor), TRUE) <= 0.05))

# boxplot(regina.boot10, trevor.boot10, main = "Smoothness = 10")



set.seed(20)
sims = 1000
regina = rep(0, sims)
trevor = rep(0, sims)
for(s in 1:sims) {
  tic("Total")
  cat("Simulation ", s, "\n")
  
  gp1 = gp1d(N, l = 20)
  gp2 = gp1d(M, l = 20)
  
  regina[s] = quality(gp1, gp2)
  trevor[s] = coverage(gp1, gp2)
  
  toc()
  cat("Regina: ", mean(regina[1:s] < 0.05), "\n")
  cat("Trevor: ", mean(trevor[1:s] < 0.05), "\n")
  cat("\n")
}

regina.boot20 = sapply(1:100, function(x) mean(sample(regina, size = length(regina), TRUE) <= 0.05))
trevor.boot20 = sapply(1:100, function(x) mean(sample(trevor, size = length(trevor), TRUE) <= 0.05))

# boxplot(regina.boot20, trevor.boot20, main = "Smoothness = 20")



boots = data.frame(Method = c(rep("Quality Index", 300), rep("KSD", 300)),
                   Smoothness = c(rep("5", 100), rep("10", 100), rep("20", 100)),
                   Size = c(regina.boot5, regina.boot10, regina.boot20,
                            trevor.boot5, trevor.boot10, trevor.boot20))

boots$Smoothness = factor(boots$Smoothness, levels = c("5", "10", "20"), ordered = TRUE)

ggplot(boots, aes(x = Smoothness, y = Size, color = Method)) +
  geom_boxplot() +
  geom_abline(intercept = 0.05, slope = 0) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5)) +
  ggtitle("Size by Method and Smoothness (N < M)")



