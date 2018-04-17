prior.gp = sim_gp(mu = 0, scale = 1)
post.gp = sim_gp(mu = 0, scale = 1)

# split em all
sim.prior.split = vapply(1:100, function(x) matsplitter(prior.gp[,,x], 10, 10),
                         FUN.VALUE = array(0, dim = c(10, 10, 9)))

sim.post.split = vapply(1:100, function(x) matsplitter(post.gp[,,x], 10, 10),
                        FUN.VALUE = array(0, dim = c(10, 10, 9)))


# find the observed wilco field
wilco.field = wilcox.field(sim.prior.split, sim.post.split)

# find the permutation distribution
perms = 100
perm.fields = matrix(0, 9, perms)
for (p in 1:perms) {
  new.fields = permute_fields(sim.prior.split, sim.post.split, seed = p)
  perm.fields[,p] = wilcox.field(new.fields[[1]], new.fields[[2]])
}

# plot the observed and permuted (as vectors)
pltmin = min(c(perm.fields, wilco.field))
pltmax = max(c(perm.fields, wilco.field))

plot(perm.fields[,1], type = "l", main=paste0("Under Null"), ylab = "Mean WILCOX", xlab = "Region", ylim = c(pltmin, pltmax))
for(p in 2:perms) {
  lines(perm.fields[,p])
}
lines(wilco.field, col = "red", lwd = 2)

# find the central regions
perm.ed = edepth_set(perm.fields)
perm.cr = central_region(perm.fields, perm.ed)

# add them to the plot
lines(perm.cr[[1]], col = "blue", lwd = 2)
lines(perm.cr[[2]], col = "blue", lwd = 2)

# add median
lines(perm.fields[perm.ed == 1], col = "blue")




# find the observed wilco field
kolm.field = ks.field(sim.prior.split, sim.post.split)

# find the permutation distribution
perms = 100
perm.fields = matrix(0, 9, perms)
for (p in 1:perms) {
  new.fields = permute_fields(sim.prior.split, sim.post.split, seed = p)
  perm.fields[,p] = ks.field(new.fields[[1]], new.fields[[2]])
}

# plot the observed and permuted (as vectors)
pltmin = min(c(perm.fields, kolm.field))
pltmax = max(c(perm.fields, kolm.field))

plot(perm.fields[,1], type = "l", main=paste0("Under Null"), ylab = "Mean K-S", xlab = "Region", ylim = c(pltmin, pltmax))
for(p in 2:perms) {
  lines(perm.fields[,p])
}
lines(kolm.field, col = "red", lwd = 2)

# find the central regions
perm.ed = edepth_set(perm.fields)
perm.cr = central_region(perm.fields, perm.ed)

# add them to the plot
lines(perm.cr[[1]], col = "blue", lwd = 2)
lines(perm.cr[[2]], col = "blue", lwd = 2)

# add median
lines(perm.fields[perm.ed == 1], col = "blue")