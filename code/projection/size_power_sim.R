
samples = 200
diffs = matrix(0, 9, samples)
post_mu = 0
set.seed(100)
for (i in 1:samples) {
  prior.gp = sim_gp(mu = 0, scale = 1)
  post.gp = sim_gp(mu = post_mu, scale = 1)
  
  # split em all
  sim.prior.split = vapply(1:100, function(x) matsplitter(prior.gp[,,x], 10, 10),
                           FUN.VALUE = array(0, dim = c(10, 10, 9)))
  
  sim.post.split = vapply(1:100, function(x) matsplitter(post.gp[,,x], 10, 10),
                          FUN.VALUE = array(0, dim = c(10, 10, 9)))
  
  
  # find the observed wilco field
  kol.field = ks.field.ez(sim.prior.split, sim.post.split)
  
  # find the permutation distribution
  perms = 100
  perm.fields = matrix(0, 9, perms)
  for (p in 1:perms) {
    new.fields = permute_fields(sim.prior.split, sim.post.split, seed = p)
    perm.fields[,p] = ks.field.ez(new.fields[[1]], new.fields[[2]], seed = p)
  }
  
  # find the central regions
  perm.ed = edepth_set(perm.fields)
  perm.cr = central_region(perm.fields, perm.ed)
  
  diffs[,i] = sapply(1:length(kol.field), function(x) 1-isbetween(kol.field[x], perm.cr[[1]][x], perm.cr[[2]][x]))
  cat("iter ", i, "\n")
  # diffs = diffs + sum(sapply(1:length(kol.field), function(x) 1-isbetween(kol.field[x], perm.cr[[1]][x], perm.cr[[2]][x])))
}

sum(diffs) / (9*samples)



samples = 200
diffs_alt = matrix(0, 9, samples)
post_mu = 0.5
for (i in 1:samples) {
  prior.gp = sim_gp(mu = 0, scale = 1)
  post.gp = sim_gp(mu = post_mu, scale = 1)
  
  # split em all
  sim.prior.split = vapply(1:100, function(x) matsplitter(prior.gp[,,x], 10, 10),
                           FUN.VALUE = array(0, dim = c(10, 10, 9)))
  
  sim.post.split = vapply(1:100, function(x) matsplitter(post.gp[,,x], 10, 10),
                          FUN.VALUE = array(0, dim = c(10, 10, 9)))
  
  
  # find the observed wilco field
  kol.field = ks.field.ez(sim.prior.split, sim.post.split)
  
  # find the permutation distribution
  perms = 100
  perm.fields = matrix(0, 9, perms)
  for (p in 1:perms) {
    new.fields = permute_fields(sim.prior.split, sim.post.split, seed = p)
    perm.fields[,p] = ks.field.ez(new.fields[[1]], new.fields[[2]])
  }
  
  # find the central regions
  perm.ed = edepth_set(perm.fields)
  perm.cr = central_region(perm.fields, perm.ed)
  
  diffs_alt[,i] = sapply(1:length(kol.field), function(x) 1-isbetween(kol.field[x], perm.cr[[1]][x], perm.cr[[2]][x]))
  cat("iter ", i, "\n")
  # diffs = diffs + sum(sapply(1:length(kol.field), function(x) 1-isbetween(kol.field[x], perm.cr[[1]][x], perm.cr[[2]][x])))
}

sum(diffs_alt) / (9*samples)
