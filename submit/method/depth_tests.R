#### K

# kolmogorov pvals
ks_cdf = function(x, n = 100) {
  if(x < 0.05) return(0)
  1 - 2*(sum(sapply(1:n, function(k) ((-1)^(k-1)) * exp(-2*(k^2)*(x^2)))))
}

kolm = function(f, g) {
  ff.xd = xdepth(f, f)
  fg.xd = xdepth(f, g)
  
  gg.xd = xdepth(g, g)
  gf.xd = xdepth(g, f)
  
  ff.cdf = sapply(ff.xd, function(y) mean(ff.xd <= y))
  gf.cdf = sapply(ff.xd, function(y) mean(gf.xd <= y))
  fg.cdf = sapply(gg.xd, function(y) mean(fg.xd <= y))
  gg.cdf = sapply(gg.xd, function(y) mean(gg.xd <= y))
  
  rate = sqrt((ncol(g)*ncol(f)) / (ncol(g) + ncol(f)))
  
  ksf = rate*max(abs(ff.cdf - gf.cdf))
  ksg = rate*max(abs(gg.cdf - fg.cdf))
  
  ks = max(ksf, ksg)
  c(ks/rate, 1-ks_cdf(ks))
}


#### QI

# Quality Index with Expected Depth
quality = function(f, g) {
  fxd = xdepth(f, f)
  gxd = xdepth(g, f)
  
  q = mean(sapply(gxd, function(y) mean(fxd <= y)))
  pval = 1 - pnorm((0.5 - q) / sqrt(((1/ncol(f)) + (1/ncol(g)))/12))
  
  c(q, pval)
}


# # K test with Extremal Depth
# kolm.ed = function(f, g) {
#   ff.xd = edepth_set(f)
#   fg.xd = apply(f, 2, function(x) edepth(x, g))
#   
#   gg.xd = edepth_set(g)
#   gf.xd = apply(g, 2, function(x) edepth(x, f))
#   
#   tf = seq(0, 1, length.out = max(1000, 3*length(ff.xd)))  
#   ff.cdf = sapply(tf, function(y) mean(ff.xd <= y))
#   gf.cdf = sapply(tf, function(y) mean(gf.xd <= y))
#   
#   tg = seq(0, 1, length.out = max(1000, 3*length(gg.xd)))
#   fg.cdf = sapply(tg, function(y) mean(fg.xd <= y))
#   gg.cdf = sapply(tg, function(y) mean(gg.xd <= y))
#   
#   rate = sqrt((ncol(g)*ncol(f)) / (ncol(g) + ncol(f)))
#   
#   ksf = rate*max(abs(ff.cdf - gf.cdf))
#   ksg = rate*max(abs(gg.cdf - fg.cdf))
#   
#   ks = max(ksf, ksg)
#   c(ks, 1-ks_cdf(ks))
# }
# 
# # K test with Projection Depth
# kolm.pd = function(f, g) {
#   ff.xd = pdepth(f, f)
#   fg.xd = pdepth(f, g)
#   
#   gg.xd = pdepth(g, g)
#   gf.xd = pdepth(g, f)
#   
#   tf = seq(0, 1, length.out = max(1000, 3*length(ff.xd)))  
#   ff.cdf = sapply(tf, function(y) mean(ff.xd <= y))
#   gf.cdf = sapply(tf, function(y) mean(gf.xd <= y))
#   
#   tg = seq(0, 1, length.out = max(1000, 3*length(gg.xd)))
#   fg.cdf = sapply(tg, function(y) mean(fg.xd <= y))
#   gg.cdf = sapply(tg, function(y) mean(gg.xd <= y))
#   
#   rate = sqrt((ncol(g)*ncol(f)) / (ncol(g) + ncol(f)))
#   
#   ksf = rate*max(abs(ff.cdf - gf.cdf))
#   ksg = rate*max(abs(gg.cdf - fg.cdf))
#   
#   ks = max(ksf, ksg)
#   c(ks, 1-ks_cdf(ks))
# }

# # K test with Expected Depth
# kolm.xd = function(f, g) {
#   ff.xd = xdepth(f, f)
#   gf.xd = xdepth(g, f)
#   
#   fg.xd = xdepth(f, g)
#   gg.xd = xdepth(g, g)
#   
#   tf = seq(0, 1, length.out = max(1000, 3*length(ff.xd)))  
#   ff.cdf = sapply(tf, function(y) mean(ff.xd <= y))
#   gf.cdf = sapply(tf, function(y) mean(gf.xd <= y))
#   
#   tg = seq(0, 1, length.out = max(1000, 3*length(gg.xd)))
#   fg.cdf = sapply(tg, function(y) mean(fg.xd <= y))
#   gg.cdf = sapply(tg, function(y) mean(gg.xd <= y))
#   
#   rate = sqrt((ncol(g)*ncol(f)) / (ncol(g) + ncol(f)))
#   
#   ksf = rate*max(abs(ff.cdf - gf.cdf))
#   ksg = rate*max(abs(gg.cdf - fg.cdf))
#   
#   ks = max(ksf, ksg)
#   c(ks, 1-ks_cdf(ks))
# }


# # Quality Index with Extremal Depth
# quality.ed = function(f, g) {
#   fxd = edepth_set(f)
#   gxd = apply(g, 2, function(x) edepth(x, f))
#   
#   r = sapply(gxd, function(y) mean(y >= fxd))
#   pval = 1 - pnorm((mean(r) - 0.5) / sqrt((1/ncol(f) + 1/ncol(g))/12))
#   
#   c(mean(r), pval)
# }
# 
# # Quality Index with Projection Depth
# quality.pd = function(f, g) {
#   fxd = pdepth(f, f)
#   gxd = pdepth(g, f)
#   
#   r = sapply(gxd, function(y) mean(y >= fxd))
#   pval = 1 - pnorm((mean(r) - 0.5) / sqrt((1/ncol(f) + 1/ncol(g))/12))
#   
#   c(mean(r), pval)
# }
