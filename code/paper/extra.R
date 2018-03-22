
#  example of approximation
prior = basis %*% (proj %*% as.vector(prior.ens[,,1]))
prior = matrix(prior, 96, 144)

plot_ly(showscale = FALSE) %>% 
   add_surface(z = ~prior.ens[,,2])

plot_ly(showscale = FALSE) %>% 
  add_surface(z = ~prior)


post = prep_post_ens(nc.post, 100)
post.coef = sapply(1:dim(post)[3], function(x) proj %*% as.vector(post[,,x]))

coef.rng = 550:575
plot(post.coef[coef.rng, 1], type = "l", xlab = "Coefficient", ylab = "Value")
for (i in 2:100) {
  lines(post.coef[coef.rng,i])
}

# USEING THE REGULAR FIELDS
prior.coef = proj %*% as.vector(prior.ens[,,10])
post.coef = sapply(1:dim(post)[3], function(x) proj %*% as.vector(post[,,x]))

coef.rng = 1:500
plot(post.coef[coef.rng, 1], type = "l", xlab = "Coefficient", ylab = "Value", main ="Regular")
for (i in 2:100) {
  lines(post.coef[coef.rng,i])
}
lines(prior.coef[coef.rng],col="red")

# look at the cumulative sum of the coefficients
plot(cumsum(post.coef[coef.rng, 1]), type = "l")
for (i in 2:100) {
  lines(cumsum(post.coef[coef.rng, i]))
}
lines(cumsum(prior.coef[coef.rng]),col="red")

# compute the depths of the regular and the integrated coefficients
prior.coef.sum = cumsum(prior.coef)
post.coef.sum = sapply(1:dim(post)[3], function(x) cumsum(post.coef[,x]))

edepth(prior.coef.sum, post.coef.sum)
edepth(prior.coef, post.coef)

# USING THE TRANSPOSED FIELDS
post = prep_post_time(nc.post, 10)
prior.coef = proj %*% as.vector(t(prior.ens[,,10]))
post.coef = sapply(1:dim(post)[3], function(x) proj %*% as.vector(t(post[,,x])))

coef.rng = 1:896
plot(post.coef[coef.rng, 1], type = "l", xlab = "Coefficient", ylab = "Value", main = "transpose")
for (i in 2:100) {
  lines(post.coef[coef.rng,i])
}
lines(prior.coef[coef.rng],col="red")

# look at the cumulative sum of the coefficients
plot(cumsum(post.coef[coef.rng, 1]), type = "l")
for (i in 2:100) {
  lines(cumsum(post.coef[coef.rng, i]))
}
lines(cumsum(prior.coef[coef.rng]),col="red")

# compute the depths of the regular and the integrated coefficients
prior.coef.sum = cumsum(prior.coef)
post.coef.sum = sapply(1:dim(post)[3], function(x) cumsum(post.coef[,x]))

edepth(prior.coef.sum, post.coef.sum)
edepth(prior.coef, post.coef)








# What if we used the SRFS??
# the magnitudes are quite close but the directions are not.
prior.srsf = sign(diff(prior.coef, 1)) * sqrt(abs(diff(prior.coef, 1)))
post.srsf = sign(diff(post.coef, 1)) * sqrt(abs(diff(post.coef, 1)))

# prior.srsf = sqrt(abs(diff(prior.coef, 1)))
# post.srsf = sqrt(abs(diff(post.coef, 1)))


plot(post.srsf[coef.rng, 1], type = "l", xlab = "Coefficient", ylab = "Value")
for (i in 2:100) {
  lines(post.srsf[coef.rng,i])
}
lines(prior.srsf[coef.rng],col="red")

