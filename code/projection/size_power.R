isbetween = function(x, a, b) {
  if(x < a) return(0)
  if(x > b) return(0)
  else return(1)
}

sim_gp = function(fields = 100, mu = 0, scale = 1, pts = 30) {
  # kernel function
  exp_kernel <- function(X,l){
    D <- plgp::distance(X)
    Sigma <- exp(-D/l)^2
  }
  
  # spatial locations to estimate
  grid = seq(-1, 1, length = pts)
  grid = expand.grid(grid, grid)
  
  # calc sigma with cov kernel
  sigma = exp_kernel(grid, l = scale)
  
  # sample fields from the GP
  gps = MASS::mvrnorm(fields, rep(mu, dim(sigma)[1]), sigma)
  gps = vapply(1:fields, function(i) matrix(gps[i,], pts, pts) 
               ,FUN.VALUE = array(0, dim = c(pts, pts)))
  return(gps)
}

diffs = 0
samples = 100
post_mu = 0
for (i in 1:samples) {
  prior.gp = sim_gp(mu = 0, scale = 1)
  post.gp = sim_gp(mu = post_mu, scale = 1)
  
  # split em all
  sim.prior.split = vapply(1:100, function(x) matsplitter(prior.gp[,,x], 10, 10),
                           FUN.VALUE = array(0, dim = c(10, 10, 9)))
  
  sim.post.split = vapply(1:100, function(x) matsplitter(post.gp[,,x], 10, 10),
                          FUN.VALUE = array(0, dim = c(10, 10, 9)))
  
  
  # find the observed wilco field
  kol.field = ks.field(sim.prior.split, sim.post.split)
  
  # find the permutation distribution
  perms = 100
  perm.fields = matrix(0, 9, perms)
  for (p in 1:perms) {
    new.fields = permute_fields(sim.prior.split, sim.post.split, seed = p)
    perm.fields[,p] = ks.field(new.fields[[1]], new.fields[[2]])
  }
  
  # find the central regions
  perm.ed = edepth_set(perm.fields)
  perm.cr = central_region(perm.fields, perm.ed)
  
  diffs = diffs + sum(sapply(1:length(kol.field), function(x) 1-isbetween(kol.field[x], perm.cr[[1]][x], perm.cr[[2]][x])))
}

diffs / (samples*9)
