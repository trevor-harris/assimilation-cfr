set.seed(1023)
sims = 500
real = rep(0, sims)
fake = rep(0, sims)
for(s in 1:sims) {
  tic("Total")
  cat("Simulation ", s, "\n")
  
  f = gp1d(400)
  g = gp1d(400)
  
  ffxd = rescale(xdepth(f, f))
  gfxd = rescale(xdepth(g, f))
  fgxd = rescale(xdepth(f, g))
  ggxd = rescale(xdepth(g, g))
  
  ffr = sapply(sort(ffxd), function(y) mean(ffxd >= y))
  gfr = sapply(sort(ffxd), function(y) mean(gfxd >= y))
  fgr = sapply(sort(ggxd), function(y) mean(fgxd >= y))
  ggr = sapply(sort(ggxd), function(y) mean(ggxd >= y))
  
  rate = sqrt((ncol(g)*ncol(f)) / (ncol(g) + ncol(f)))
  
  ksf1 = rate*max(abs(ffr - gfr))
  ksg1 = rate*max(abs(fgr - ggr))
  
  real[s] = max(ksf1, ksg1)
  
  
  tf = 1:length(ffxd) / length(ffxd)
  ffr = sapply(tf, function(y) mean(ffxd >= y))
  gfr = sapply(tf, function(y) mean(gfxd >= y))
  
  tg = 1:length(ggxd) / length(ggxd)
  fgr = sapply(tg, function(y) mean(fgxd >= y))
  ggr = sapply(tg, function(y) mean(ggxd >= y))
  
  ksf = rate*max(abs(ffr - gfr))
  ksg = rate*max(abs(fgr - ggr))
  
  fake[s] = max(ksf, ksg)
  
  toc()
  cat("Diff: ", mean((real[1:s] - fake[1:s])^2), "\n")
  cat("\n")
}