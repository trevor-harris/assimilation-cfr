prior.gp = sim_gp(mu = 0, scale = 1)
post.gp = sim_gp(mu = 0, scale = 0.01)

# split em all
sim.prior.split = vapply(1:100, function(x) matsplitter(prior.gp[,,x], 10, 10),
                         FUN.VALUE = array(0, dim = c(10, 10, 9)))

sim.post.split = vapply(1:100, function(x) matsplitter(post.gp[,,x], 10, 10),
                        FUN.VALUE = array(0, dim = c(10, 10, 9)))


prior.split = sim.prior.split
post.split = sim.post.split

iter = 1000

nlat = dim(prior.split)[1]
nlon = dim(prior.split)[2]
regions = dim(prior.split)[3]
ens = dim(prior.split)[4]

prior.proj = matrix(0, ens, iter)
post.proj = matrix(0, ens, iter)
wilco.field = 1:regions
proj = matrix(rnorm(nlat*nlon*iter), nlat*nlon, iter)

for(r in 1:regions) {
  for(e in 1:ens) {
    prior.proj[e,] = as.vector(as.vector(prior.split[,,r,e]) %*% proj)
    post.proj[e,] = as.vector(as.vector(post.split[,,r,e]) %*% proj)
  }
}

ad.test(prior.proj[,1], post.proj[,1])
ks.test(prior.proj[,1], post.proj[,1], alternative = "two.sided")



projected = data.frame(value = c(prior.proj[,1], post.proj[,1]),
                      label = c(rep("prior", times = 100), rep("post", times = 100)))

ggplot(data = projected, aes(value, fill = label)) +
  geom_density(alpha = 0.2) +
  ggtitle("Prior vs Posterior under Alternative with scale = 0.01") +
  theme(plot.title = element_text(hjust = 0.5))
  

ggplot(data = projected, aes(value, color = label)) +
  stat_ecdf() +
  ggtitle("Prior vs Posterior under Alternative with scale = 0.01") +
  theme(plot.title = element_text(hjust = 0.5))

# plot(prior.proj[,1], type = "l", col = "red")
# lines(post.proj[,1], col = "blue")

