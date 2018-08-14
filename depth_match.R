rm(list = ls()); gc()

library(extdepth)
library(fda)

weather = CanadianWeather

temp = as.matrix(weather$dailyAv[,,1])

atlantic = temp[,weather$region == "Atlantic"]
arctic = temp[,weather$region == "Arctic"]
continental = temp[,weather$region == "Continental"]

plt_2funs(continental, atlantic)
ks.dist(atlantic, continental)


ks_pval = function(t, n = 10) {
  2*(sum(sapply(1:n, function(x) (-1)^(x-1) * exp(-2*(x^2)*t))))
}

dks = function(f, g) {

  fed = edepth_set(f)
  ged = rep(0, ncol(g))
  for(i in ncol(g)) {
    ged[i] = edepth(g[,i], cbind(g[,i],f))
  }
  
  f.surv = rev(c(0, sort(ed)))
  g.cdf = sapply(1:length(f.surv), function(x) mean(ged >= f.surv[x]))
  f.cdf = rev(f.surv)
  
  ks = max(abs(f.cdf - g.cdf))
  ks_pval(sqrt(ncol(g))*ks)
}

p1 = depth_ks(arctic, continental)
p2 = depth_ks(arctic, atlantic)
p3 = depth_ks(atlantic, continental)
depth_ks(arctic, continental)
depth_ks(continental, arctic)
