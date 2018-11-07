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
# meandepth = function(gmat, fmat) {
#   apply(gmat, 2, function(x) mean(depth(x, fmat)))
# }
ks_cdf = function(x, n = 10) {
  1 - 2*(sum(sapply(1:n, function(k) ((-1)^(k-1)) * exp(-2*(k^2)*(x^2)))))
}

# ks.mean = function(f, g) {
#   # how well does g fit in with f?
#   ff.ed = meandepth(f, f)
#   gf.ed = meandepth(g, f)
#   
#   # how well does f fit in with g?
#   fg.ed = meandepth(f, g)
#   gg.ed = meandepth(g, g)
#   
#   f.surv = rev(c(0, sort(ff.ed)))
#   gf.cdf = sapply(f.surv, function(x) mean(gf.ed > x))
#   
#   g.surv = rev(c(0, sort(gg.ed)))
#   fg.cdf = sapply(g.surv, function(x) mean(fg.ed > x))
#   
#   # how well do they fit with themselves?
#   uni = seq(0, 1, length.out = length(gf.cdf))
#   
#   rate = sqrt((ncol(g)*ncol(f)) / (ncol(g) + ncol(f)))
#   
#   ksf = rate * max(abs(uni - gf.cdf))
#   ksg = rate * max(abs(uni - fg.cdf))
#   
#   # 1 - ks_cdf(ksf)*ks_cdf(ksg)
#   1 - ks_cdf(max(ksf, ksg))^2
# }


# TESTING USING THE ACTUAL CR
xdepth = function(gmat, fmat) {
  apply(gmat, 2, function(x) mean(depth(x, fmat)))
}
central_region <- function(fmat, edepths, alpha = 0.05) {
  # Takes a list of extremal depths and an alpha level and returns the functions corresponding to
  # the lower and upper bounds. Also returns the median function
  
  # fmat = f
  # edepths = ed
  # alpha = 1000/1000
  
  # filter out functions outside the alpha level
  fset = as.matrix(fmat[,edepths > alpha])
  if (ncol(fset) == 0) fset = as.matrix(fmat[,which.max(edepths)])
  
  # lower
  lower = sapply(1:nrow(fset), function(x) min(fset[x,]))
  
  # upper
  upper = sapply(1:nrow(fset), function(x) max(fset[x,]))
  
  return(list(lower = lower,
              upper = upper))
}

# ks_cdf = function(t, n = 10) {
#   1 - 2*(sum(sapply(1:n, function(x) (-1)^(x-1) * exp(-2*(x^2)*t^2))))
# }
# coverage = function(f, g, s = 2*max(ncol(f), ncol(g))) {
#   f = gp1d(50)
#   g = gp1d(50)
#   s = 2*max(ncol(f), ncol(g))
#   
#   f.n = ncol(f)
#   g.n = ncol(g)
#   
#   f.xd = xdepth(f, f)
#   g.xd = xdepth(g, g)
#   
#   f.xd = edepth_set(f)
#   g.xd = edepth_set(g)
#   
#   f.fcr = rep(0, s)
#   g.fcr = rep(0, s)
#   f.gcr = rep(0, s)
#   g.gcr = rep(0, s)
#   
#   for (a in 1:s) {
#     f.cr = central_region(f, f.xd, 1-(a/s))
#     
#     f.fcr[a] = mean(sapply(1:f.n, function(x) min((f[,x] >= f.cr$lower) & (f[,x] <= f.cr$upper))))
#     g.fcr[a] = mean(sapply(1:g.n, function(x) min((g[,x] >= f.cr$lower) & (g[,x] <= f.cr$upper))))
#     
#     g.cr = central_region(g, g.xd, 1-(a/s))
#     
#     f.gcr[a] = mean(sapply(1:f.n, function(x) min((f[,x] >= g.cr$lower) & (f[,x] <= g.cr$upper))))
#     g.gcr[a] = mean(sapply(1:g.n, function(x) min((g[,x] >= g.cr$lower) & (g[,x] <= g.cr$upper))))
#   }
# 
#   plot(f.fcr, type = "l")
#   lines(g.fcr)
# 
#   plot(g.gcr, type = "l")
#   lines(f.gcr)
#   # 
#   # plt_funs(f, g)
# 
#   rate = sqrt((ncol(g)*ncol(f)) / (ncol(g) + ncol(f)))
#   
#   f.ks = rate * max(abs(f.fcr - g.fcr))
#   g.ks = rate * max(abs(f.gcr - g.gcr))
# 
#   # 1 - ks_cdf(f.ks)
#     
#   1 - ks_cdf(max(f.ks, g.ks))^2
# }
# coverage = function(f, g, s = 2*max(ncol(f), ncol(g))) {
#   f = gp1d()
#   g = gp1d()
#   s = 2*max(ncol(f), ncol(g))
#   
#   f.n = ncol(f)
#   g.n = ncol(g)
#   
#   ff.xd = xdepth(f, f)
#   gf.xd = xdepth(f, g)
#   
#   fg.xd = xdepth(g, f)
#   gg.xd = xdepth(g, g)
# 
#   f.gcdf = rep(0, s)
#   g.fcdf = rep(0, s)
#   
#   for (a in 1:s) {
#     alpha = 1-a/s
#     f.gcdf[a] = mean(sapply(1:f.n, function(x) (ff.xd[x] > alpha) & (gf.xd[x] > alpha)))
#     g.fcdf[a] = mean(sapply(1:g.n, function(x) (gg.xd[x] > alpha) & (fg.xd[x] > alpha)))
#   }
#   
#   # 1 - ks_cdf(sqrt(ncol(f)) * max(abs(seq(0, 1, length.out = s) - f.fcr)))
#   
#   rate = sqrt((ncol(g)*ncol(f)) / (ncol(g) + ncol(f)))
#   
#   f.ks = rate * max(abs(f.gcdf - seq(0, 1, length.out = s)))
#   g.ks = rate * max(abs(g.fcdf - seq(0, 1, length.out = s)))
#   
#   1 - ks_cdf(max(f.ks, g.ks))^2
#   
#   # plot(f.gcdf, type = "l")
#   # lines(seq(0, 1, length.out = s))
#   # 
#   # plot(g.fcdf, type = "l")
#   # lines(seq(0, 1, length.out = s))
#   # 
#   # plt_funs(f, g)
#   # 
#   # rate = sqrt((ncol(g)*ncol(f)) / (ncol(g) + ncol(f)))
#   # 
#   # f.ks = rate * max(abs(f.fcr - g.fcr))
#   # g.ks = rate * max(abs(f.gcr - g.gcr))
#   # 
#   # 1 - ks_cdf(f.ks)
#   # 1 - ks_cdf(max(f.ks, g.ks))^2
# }
# 
# coverage = function(f, g, s = 2) {
#   # f = gp1d()
#   # g = gp1d()
#   s = s*max(ncol(f), ncol(g))
#   
#   f.n = ncol(f)
#   g.n = ncol(g)
#   
#   ff.xd = rdepth(f, f)
#   gf.xd = rdepth(g, f)
#   
#   gg.xd = rdepth(g, g)
#   fg.xd = rdepth(f, g)
#   
#   f.gcdf = sapply(sort(ff.xd), function(a) mean(gf.xd > a))
#   g.fcdf = sapply(sort(gg.xd), function(a) mean(fg.xd > a))
#   
#   rate = sqrt((g.n*f.n) / (f.n + g.n))
#   f.ks = rate * max(abs(f.gcdf - seq(1, 0, length.out = f.n)))
#   g.ks = rate * max(abs(g.fcdf - seq(1, 0, length.out = g.n)))
#   
#   1 - ks_cdf(max(f.ks, g.ks))^2
# }
# 
# # Ks test version
# coverage = function(f, g, s = 2) {
#   f = gp1d()
#   g = gp1d()
#   s = 2*max(ncol(f), ncol(g))
#   
#   f.n = ncol(f)
#   g.n = ncol(g)
#   
#   ff.xd = rdepth(f, f)
#   gf.xd = rdepth(g, f)
#   
#   gg.xd = rdepth(g, g)
#   fg.xd = rdepth(f, g)
#   
#   t.test(ff.xd, gf.xd)
#   
#   alpha = seq(1, 0, length.out = s)
#   
#   ff.cdf = sapply(alpha, function(a) mean(ff.xd > a))
#   gf.cdf = sapply(alpha, function(a) mean(gf.xd > a))
#   fg.cdf = sapply(alpha, function(a) mean(fg.xd > a))
#   gg.cdf = sapply(alpha, function(a) mean(gg.xd > a))
#   
#   rate = sqrt((g.n*f.n) / (f.n + g.n))
#   f.ks = sqrt(f.n) * max(abs(ff.cdf - gf.cdf))
#   g.ks = sqrt(g.n) * max(abs(gg.cdf - fg.cdf))
#   
#   plot(ff.cdf, type = "l")
#   lines(gf.cdf, col = "red")
# 
#   plot(gg.cdf, type = "l")
#   lines(fg.cdf, col = "red")
#   
#   1 - ks_cdf(max(f.ks, g.ks))^2
# }

# coverage = function(f, g, s = 2) {
#   # f = gp1d()
#   # g = gp1d()
#   s = s*max(ncol(f), ncol(g))
#   
#   f.n = ncol(f)
#   g.n = ncol(g)
#   
#   ff.xd = xdepth(f, f)
#   gf.xd = xdepth(f, g)
#   
#   fg.xd = xdepth(g, f)
#   gg.xd = xdepth(g, g)
#   
#   f.gcdf = rep(0, s)
#   g.fcdf = rep(0, s)
#   
#   for (a in 1:s) {
#     alpha = 1-a/s
#     f.gcdf[a] = mean(sapply(1:f.n, function(x) (ff.xd[x] > alpha) & (gf.xd[x] > alpha)))
#     g.fcdf[a] = mean(sapply(1:g.n, function(x) (gg.xd[x] > alpha) & (fg.xd[x] > alpha)))
#   }
#   
#   # 1 - ks_cdf(sqrt(ncol(f)) * max(abs(seq(0, 1, length.out = s) - f.fcr)))
#   
#   rate = sqrt((ncol(g)*ncol(f)) / (ncol(g) + ncol(f)))
#   
#   f.ks = sqrt((f.n + g.n)) * max(abs(f.gcdf - seq(0, 1, length.out = s)))
#   g.ks = sqrt((f.n + g.n)) * max(abs(g.fcdf - seq(0, 1, length.out = s)))
#   
#   1 - ks_cdf(max(f.ks, g.ks))^2
# }

# 1way t test version
coverage1 = function(f, g, s = 2) {
  ff.xd = xdepth(f, f)
  gf.xd = xdepth(g, f)
  
  t.test(ff.xd, gf.xd)$p.value
}

coverageks = function(f, g, s = 2) {
  ff.xd = xdepth(f, f)
  gf.xd = xdepth(g, f)

  alpha = seq(1, 0, length.out = s*ncol(f))
  ff.cdf = sapply(alpha, function(a) mean(ff.xd > a))
  gf.cdf = sapply(alpha, function(a) mean(gf.xd > a))

  rate = sqrt((ncol(g)*ncol(f)) / (ncol(g) + ncol(f)))
  stat = max(abs(ff.cdf - gf.cdf))
  1 - ks_cdf(rate * stat)
}

# set.seed(1023)
# sims = 1000
# cov1 = rep(0, sims)
# covk = rep(0, sims)
# for(s in 1:sims) {
#   tic("Total")
#   cat("Simulation ", s, "\n")
#   
#   # gp1 = gp1d(200, mu = 0, pts = 50, l = 10)
#   # gp2 = gp1d(50, mu = 0, pts = 50, l = 10)
#   
#   gp1 = gp1d()
#   gp2 = gp1d()
#   
#   cov1[s] = coverage1(gp1, gp2, 2)
#   covk[s] = coverageks(gp1, gp2, 2)
#   
#   toc()
#   cat(mean(cov1[1:s] < 0.05), "\n")
#   cat(mean(covk[1:s] < 0.05), "\n")
#   cat("\n")
# }
# mean(cov1 < 0.05)
# mean(covk < 0.05)
# plot(cov1)
# plot(covk)

nscore = function(xd1, xd2) {
  xd = c(xd1, xd2)
  nxd = qnorm(rank(xd) / (length(xd) + 1))
  return(list(xd1 = nxd[1:length(xd1)],
              xd2 = nxd[(length(xd1) + 1):length(xd)]))
}

coverage2 = function(f, g, s = 2) {
  ff.xd = xdepth(f, f)
  gf.xd = xdepth(g, f)

  gg.xd = xdepth(g, g)
  fg.xd = xdepth(f, g)
  
  stat = max(t.test(ff.xd, gf.xd)$statistic, t.test(gg.xd, fg.xd)$statistic)
  1 - pnorm(stat)^2
}


ncov = function(f, g, s = 2) {
  ff.xd = xdepth(f, f)
  gf.xd = xdepth(g, f)

  gg.xd = xdepth(g, g)
  fg.xd = xdepth(f, g)
  
  fscores = nscore(ff.xd, gf.xd)
  gscores = nscore(gg.xd, fg.xd)
    
  stat = max(t.test(fscores$xd1, fscores$xd2)$statistic,
             t.test(gscores$xd1, gscores$xd2)$statistic)
  1 - pnorm(stat)^2
}

coverage3 = function(f, g, s = 2) {
  ff.xd = xdepth(f, f)
  gf.xd = xdepth(g, f)
  
  ff.xd = xdepth(g, g)
  gf.xd = xdepth(f, g)
  
  f.stat = t.test(ff.xd, gf.xd)$statistic
  g.stat = t.test(ff.xd, gf.xd)$statistic
  
  f.weight = ncol(f) / (ncol(f) + ncol(g))
  g.weight = ncol(g) / (ncol(f) + ncol(g))
  
  stat = (f.weight * f.stat) + (g.weight * g.stat)
  1 - pnorm(stat)
}


gp1 = gp1d()^2
gp2 = gp1d()

plt_funs(gp1, gp2)

coverage1(gp1d(), gp1d())
coverage2(gp1d(), gp1d())

set.seed(1023)
sims = 500
cov1 = rep(0, sims)
cov2 = rep(0, sims)
cov3 = rep(0, sims)
for(s in 1:sims) {
  tic("Total")
  cat("Simulation ", s, "\n")
  
  # gp1 = gp1d(200, mu = 0, pts = 50, l = 10)
  # gp2 = gp1d(50, mu = 0, pts = 50, l = 10)
  
  gp1 = gp1d(fields = 100, l = 10)
  gp2 = gp1d(fields = 100, l = 10)^2
  
  cov1[s] = coverage1(gp1, gp2, 2) 
  cov2[s] = coverage2(gp1, gp2, 2)
  cov3[s] = ncov(gp1, gp2, 2) 
  
  # cov1[s] = coverage1(gp1d(), gp1d(25), 2) 
  # cov2[s] = coverage2(gp1d(), gp1d(25), 2) 
  
  toc()
  cat(mean(cov1[1:s] < 0.05), "\n")
  cat(mean(cov2[1:s] < 0.05), "\n")
  cat(mean(cov3[1:s] < 0.05), "\n")
  cat("\n")
}
mean(cov1 < 0.05)
mean(cov2 < 0.05)
mean(cov3 < 0.05)

plot(cov1)
plot(cov2)
plot(cov3)


plot(cov1, cov2)


# f = gp1d(50)
# g = gp1d(50)
# ed = edepth_set(f)
# 
# central_region <- function(fmat, edepths, alpha = 0.05) {
#   # Takes a list of extremal depths and an alpha level and returns the functions corresponding to
#   # the lower and upper bounds. Also returns the median function
#   
#   # fmat = f
#   # edepths = ed
#   # alpha = 1000/1000
#   
#   # filter out functions outside the alpha level
#   fset = as.matrix(fmat[,edepths > alpha])
#   if (ncol(fset) == 0) fset = as.matrix(fmat[,which.max(edepths)])
#   
#   # lower
#   lower = sapply(1:nrow(fset), function(x) min(fset[x,]))
#   
#   # upper
#   upper = sapply(1:nrow(fset), function(x) max(fset[x,]))
#   
#   return(list(lower = lower,
#               upper = upper))
# }
# 
# f.incr = rep(0, 1000)
# g.incr = rep(0, 1000)
# for (a in 1:1000) {
#   cr = central_region(f, ed, a/1000)
#   
#   f.incr[a] = 1-mean(sapply(1:50, function(x) min((f[,x] >= cr$lower) & (f[,x] <= cr$upper))))
#   g.incr[a] = 1-mean(sapply(1:50, function(x) min((g[,x] >= cr$lower) & (g[,x] <= cr$upper))))
# }
# 
# plot(f.incr)
# 
# ks = max(abs(f.incr - g.incr))

# plt_funs(gp1, gp2)
