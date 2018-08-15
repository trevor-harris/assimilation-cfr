rm(list = ls()); gc()

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

# median depth
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
mdepth = function(gmat, fmat) {
  apply(gmat, 2, function(x) median(depth(x, fmat)))
}
ks.med = function(x, y) {
  mx = mdepth(x, x)
  my = mdepth(y, x)
  
  ks.test(mx, my)$p.value
}

# ED
edepth_multi = function(gmat, fmat_sorted) {
  fmat = fmat_sorted
  fdepths = depth_set(fmat)
  gdepths = apply(gmat, 2, function(x) depth(x, fmat))
  
  # get the allowed r values (for calculating the dCDF)
  rvals = sort(unique(c(gdepths, fdepths)))
  
  ged = rep(1, ncol(gmat))
  for(g in 1:ncol(gmat)) {
    for(f in 1:ncol(fmat)) {
      for(r in rvals) {
        dg = sum(gdepths[,g] <= r)
        df = sum(fdepths[,f] <= r)
        if(dg != df) {
          break;
        }
      }
      if (dg > df) {
        ged[g] = (f-1) / ncol(gmat)
        break;
      }
    }
  }
  ged
}
ks_pval = function(t, n = 10) {
  2*(sum(sapply(1:n, function(x) (-1)^(x-1) * exp(-2*(x^2)*t))))
}
ks.depth = function(f, g) {
  
  f = edepth_sort(f)
  fed = (1:ncol(f)) / ncol(f)
  ged = edepth_multi(g, f)
  
  f.surv = rev(c(0, sort(fed)))
  g.cdf = sapply(1:length(f.surv), function(x) mean(ged >= f.surv[x]))
  f.cdf = rev(f.surv)
  
  ks = max(abs(f.cdf - g.cdf))
  ks_pval(sqrt(ncol(g))*ks)
}

# median depth
ks.mcdf = function(f, g) {
  fed = mdepth(f, f)
  ged = mdepth(g, f)
  
  f.surv = rev(c(0, sort(fed)))
  g.cdf = sapply(1:length(f.surv), function(x) mean(ged >= f.surv[x]))
  f.cdf = sapply(1:length(f.surv), function(x) mean(fed >= f.surv[x]))
  
  ks = max(abs(f.cdf - g.cdf))
  ks_pval(sqrt(ncol(g))*ks)
}


sims = 100
kdist = rep(0, sims)
kmedi = rep(0, sims)
kdept = rep(0, sims)
kmcdf = rep(0, sims)
for(s in 1:sims) {
  tic("Total")
  cat("Simulation ", s, "\n")
  
  gp1 = gp1d()
  gp2 = gp1d()
  
  kdist[s] = ks.dist(gp1, gp2)
  kmedi[s] = ks.med(gp1, gp2)
  kdept[s] = ks.depth(gp1, gp2)
  kmcdf[s] = ks.mcdf(gp1, gp2)
  
  toc()
  cat("\n")
}

mean(kdist < 0.05)
mean(kmedi < 0.05)
mean(kdept < 0.05)
mean(kmcdf < 0.05)

plot(kdist)
plot(kmedi)
plot(kdept)
plot(kmcdf)

boxplot(kdist, kmedi, kdept, kmcdf)


