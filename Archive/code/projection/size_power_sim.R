
set.seed(100)
sim_size_power = function(post_mu, samples = 100) {
  diffs = matrix(0, 9, samples)
  for (i in 1:samples) {
    prior.gp = sim_gp(mu = 0, l = 1)
    post.gp = sim_gp(mu = post_mu, l = 1)
    
    # split em all
    sim.prior.split = vapply(1:100, function(x) matsplitter(prior.gp[,,x], 10, 10),
                             FUN.VALUE = array(0, dim = c(10, 10, 9)))
    
    sim.post.split = vapply(1:100, function(x) matsplitter(post.gp[,,x], 10, 10),
                            FUN.VALUE = array(0, dim = c(10, 10, 9)))
    
    
    # find the observed kst field
    kol.field = kst.field(sim.prior.split, sim.post.split)

    # find the permutation distribution
    perm.fields = kst.permute(sim.prior.split, sim.post.split, 100, 1)
    
    # find the central regions
    perm.ed = edepth_set(perm.fields)
    perm.cr = central_region(perm.fields, perm.ed)
    
    diffs[,i] = sapply(1:length(kol.field), function(x) 1-isbetween(kol.field[x], perm.cr[[1]][x], perm.cr[[2]][x]))
    cat("iter ", i, "\n")
  }
  
  return(diffs)
}

# under the null
diff_null = sim_size_power(0)
sum(diff_null) / (9*100)


# under a single alternative
diff_alt = sim_size_power(0.5)
sum(diff_alt) / (9*100)

# under the alternatives
diff_alt = array(0, dim = c(9, 200, 10))

for (j in 1:10) {
  cat("Alt: ", j/10, "\n")
  diff_alt[,,j] = sim_size_power(j/10)
}

for (j in 1:10) {
  cat(sum(diff_alt[,,j])/(200*9), "\n")
}

save.image(file = "sims.RData")

