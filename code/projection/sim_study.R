#### SIMULATION STUDY --- PHASE 1 SAME FUNCTIONS ####

# generate spatial ensemble fields from the same distributions
sim.prior = array(0, dim = c(100, 100, 100))
for (i in 1:100) {
  sim.prior[,,i] = runif(1, 0.75, 1.25) * c(outer(dnorm((1:100), 25, 25), dnorm((1:100), 75, 25), "*")) + 
                   runif(1, 0.75, 1.25) * c(outer(dnorm((1:100), 75, 25), dnorm((1:100), 25, 25), "*"))
  sim.prior[,,i] = sim.prior[,,i] * 10000 * matrix(rnorm(100*100, 25, 1), 100, 100)
}

sim.post = array(0, dim = c(100, 100, 100))
for (i in 1:100) {
  sim.post[,,i] = runif(1, 0.75, 1.25) * c(outer(dnorm((1:100), 25, 25), dnorm((1:100), 75, 25), "*")) + 
                  runif(1, 0.75, 1.25) * c(outer(dnorm((1:100), 75, 25), dnorm((1:100), 25, 25), "*"))
  sim.post[,,i] = sim.post[,,i] * 10000 * matrix(rnorm(100*100, 25, 1), 100, 100)
}


# split em all
sim.prior.split = vapply(1:100, function(x) matsplitter(sim.prior[,,x], 20, 20),
                     FUN.VALUE = array(0, dim = c(20, 20, 25)))

sim.post.split = vapply(1:100, function(x) matsplitter(sim.post[,,x], 20, 20),
                     FUN.VALUE = array(0, dim = c(20, 20, 25)))


# find the observed wilco field
wilco.field = wilcox.field(sim.prior.split, sim.post.split)

# find the permutation distribution
perms = 200
perm.fields = matrix(0, 25, perms)
for (p in 1:perms) {
  new.fields = permute_fields(sim.prior.split, sim.post.split, seed = p)
  perm.fields[,p] = wilcox.field(new.fields[[1]], new.fields[[2]])
}

# plot the observed and permuted (as vectors)
plot(perm.fields[,1], type = "l", main="Regular Test - Same", ylab = "Mean Wilcoxon", xlab = "Region", ylim = c(min(perm.fields), max(perm.fields)))
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



# fields.perm = permute_fields(sim.prior.split, sim.post.split, seed = 1)
# 
# plot_ly(showscale = FALSE) %>%
#   add_surface(z = ~sim.prior.split[,,1,3]) %>%
#   add_surface(z = ~fields.perm[[1]][,,1,3])


# plot_ly(showscale = FALSE) %>%
#   add_surface(z = ~sim.prior[,,1]) %>%
#   add_surface(z = ~sim.prior[,,10]) %>%
#   add_surface(z = ~sim.prior[,,50]) %>%
#   add_surface(z = ~sim.prior[,,80]) %>%
#   add_surface(z = ~sim.post[,,1]) %>%
#   add_surface(z = ~sim.post[,,10]) %>%
#   add_surface(z = ~sim.post[,,50]) %>%
#   add_surface(z = ~sim.post[,,80])
  



#### SIMULATION STUDY --- PHASE 2 DIFFERENT FUNCTIONS ####

# generate spatial ensemble fields from the same distributions
sim.prior = array(0, dim = c(100, 100, 100))
for (i in 1:100) {
  sim.prior[,,i] = runif(1, 0.75, 1.25) * c(outer(dnorm((1:100), 25, 25), dnorm((1:100), 75, 25), "*")) + 
                   runif(1, 0.75, 1.25) * c(outer(dnorm((1:100), 75, 25), dnorm((1:100), 25, 25), "*"))
  sim.prior[,,i] = sim.prior[,,i] * 10000 * matrix(rnorm(100*100, 25, 1), 100, 100)
}

sim.post = array(0, dim = c(100, 100, 100))
for (i in 1:100) {
  sim.post[,,i] = runif(1, 0.75, 1.25) * c(outer(dnorm((1:100), 50, 25), dnorm((1:100), 50, 25), "*")) + 
                  runif(1, 0.75, 1.25) * c(outer(dnorm((1:100), 75, 25), dnorm((1:100), 25, 25), "*"))
  sim.post[,,i] = sim.post[,,i] * 10000 * matrix(rnorm(100*100, 25, 1), 100, 100)
}

# matrix(rnorm(100*100, 25, 1), 100, 100)

# split em all
sim.prior.split = vapply(1:100, function(x) matsplitter(sim.prior[,,x], 20, 20),
                         FUN.VALUE = array(0, dim = c(20, 20, 25)))

sim.post.split = vapply(1:100, function(x) matsplitter(sim.post[,,x], 20, 20),
                        FUN.VALUE = array(0, dim = c(20, 20, 25)))


# find the observed wilco field
wilco.field = wilcox.field(sim.prior.split, sim.post.split)

# find the permutation distribution
perms = 200
perm.fields = matrix(0, 25, perms)
for (p in 1:perms) {
  new.fields = permute_fields(sim.prior.split, sim.post.split, seed = p)
  perm.fields[,p] = wilcox.field(new.fields[[1]], new.fields[[2]])
}

# plot the observed and permuted (as vectors)
plot(perm.fields[,1], type = "l", main="Regular Test - Different", ylab = "Mean Wilcoxon", xlab = "Region", ylim = c(min(perm.fields), max(perm.fields)))
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



# fields.perm = permute_fields(sim.prior.split, sim.post.split, seed = 1)
# 
# plot_ly(showscale = FALSE) %>%
#   add_surface(z = ~sim.prior.split[,,1,3]) %>%
#   add_surface(z = ~fields.perm[[1]][,,1,3])


# plot_ly(showscale = FALSE) %>%
#   add_surface(z = ~sim.prior[,,1]) %>%
#   add_surface(z = ~sim.prior[,,10]) %>%
#   add_surface(z = ~sim.prior[,,50]) %>%
#   add_surface(z = ~sim.prior[,,80]) %>%
#   add_surface(z = ~sim.post[,,1]) %>%
#   add_surface(z = ~sim.post[,,10]) %>%
#   add_surface(z = ~sim.post[,,50]) %>%
#   add_surface(z = ~sim.post[,,80])


