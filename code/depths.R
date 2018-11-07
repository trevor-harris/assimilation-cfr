#### DEPTHS
# Expected Depth (XD)
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

# Projection Depth (PD)
pdepth = function(gmat, fmat, proj = 1000) {
  u = matrix(runif(nrow(gmat) * proj, -1, 1), nrow(gmat), proj)
  u = apply(u, 2, function(x) x / sqrt(sum(x^2)))
  
  Fu = t(u) %*% fmat
  Gu = t(u) %*% gmat
  
  Fu.mu = apply(Fu, 1, median)
  Fu.sig = apply(Fu, 1, mad)
  
  out = t(abs(t(Gu) - Fu.mu) / Fu.sig)
  
  return(1 / (1 + apply(out, 2, max)))
}

# extremal Depth(ED)
library(extdepth)