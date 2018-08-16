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
ks.depth = function(f, g) {
  
  f = edepth_sort(f)
  fed = (1:ncol(f)) / ncol(f)
  ged = edepth_multi(g, f)
  
  f.surv = rev(c(0, sort(fed)))
  g.cdf = sapply(1:length(f.surv), function(x) mean(ged >= f.surv[x]))
  f.cdf = rev(f.surv)
  
  ks = max(abs(f.cdf - g.cdf))
  rate = sqrt((ncol(g)*ncol(f)) / (ncol(g) + ncol(f)))
  ks_pval(rate*ks)
}

meandepth = function(gmat, fmat) {
  apply(gmat, 2, function(x) mean(depth(x, fmat)))
}
ks.mean = function(f, g) {
  fed = meandepth(f, f)
  ged = meandepth(g, f)
  
  f.surv = rev(c(0, sort(fed)))
  g.cdf = sapply(1:length(f.surv), function(x) mean(ged >= f.surv[x]))
  f.cdf = sapply(1:length(f.surv), function(x) mean(fed >= f.surv[x]))
  
  # plot(g.cdf, type = "l", col = "blue")
  # lines(f.cdf, col = "red")
  
  ks = max(abs(f.cdf - g.cdf))
  rate = sqrt((ncol(g)*ncol(f)) / (ncol(g) + ncol(f)))
  ks_pval(rate*ks)
}

cos_mu = cos(seq(-pi, pi, length.out = 50))/2
sin_mu = sin(seq(-pi, pi, length.out = 50))/2

gp1 = gp1d(mu = sin_mu, pts = 100, l = 500)
gp2 = gp1d(mu = cos_mu, pts = 100, l = 500)

plt_funs(gp1, gp2)
ks.dist(gp1, gp2)



f = gp1d(500, mu = sin_mu, pts = 100, l = 10)
g = gp1d(500, mu = sin_mu, pts = 100, l = 10)

plt_funs(f, g)
ks.dist(f, g)
ks.med(f, g)


cos_mu = cos(seq(-pi, pi, length.out = 50)) / 4
sin_mu = sin(seq(-pi, pi, length.out = 50)) / 4

sims = 1000
kdist = rep(0, sims)
for(s in 1:sims) {
  tic("Total")
  cat("Simulation ", s, "\n")
  
  gp1 = gp1d(500, mu = cos_mu, pts = 50, l = 10)
  gp2 = gp1d(500, mu = sin_mu, pts = 50, l = 10)
  
  kdist[s] = ks.mean(gp1, gp2) 
  
  toc()
  cat("\n")
}
mean(kdist < 0.05)
plot(kdist)

plt_funs(gp1, gp2)


