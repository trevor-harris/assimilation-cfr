
#  example of approximation
prior = basis %*% (proj %*% as.vector(prior.ens[,,1]))
prior = matrix(prior, 96, 144)
post = prep_post_ens(nc.post, 20)


plot_ly(showscale = FALSE) %>%
   add_surface(z = ~prior.ens[,,1])

plot_ly(showscale = FALSE) %>%
  add_surface(z = ~prior)


basis.fun = matrix(basis[,450]+basis[,660], 96, 144)
plot_ly(showscale = FALSE) %>%
  add_surface(z = ~basis.fun)

# post = prep_post_ens(nc.post, 20)
# post.coef = sapply(1:dim(post)[3], function(x) proj %*% as.vector(post[,,x]))
# 
# coef.rng = 550:575
# plot(post.coef[coef.rng, 1], type = "l", xlab = "Coefficient", ylab = "Value")
# for (i in 2:100) {
#   lines(post.coef[coef.rng,i])
# }

# USEING THE REGULAR FIELDS
prior.coef = proj %*% as.vector(t(prior.ens[,,20]))

post = prep_post_ens(nc.post, 20)
post.coef = sapply(1:dim(post)[3], function(x) proj %*% as.vector(t(post[,,x])))

coef.rng = 1:40
plot(post.coef[coef.rng, 1], type = "l", xlab = "Coefficient", ylab = "Value", main ="Regular")
for (i in 2:100) {
  lines(post.coef[coef.rng,i])
}
lines(prior.coef[coef.rng],col="red")

# compute the depths of the regular and the integrated coefficients
prior.int = prior.coef - mean(prior.coef)
post.int = post.coef - colMeans(post.coef)
prior.coef.sum = cumsum(prior.int)
post.coef.sum = sapply(1:dim(post)[3], function(x) cumsum(post.int[,x]))

# edepth(prior.coef.sum, post.coef.sum)
# edepth(prior.coef, post.coef)

plot(post.coef.sum[coef.rng, 1], type = "l", xlab = "Coefficient", ylab = "Value", main ="Integrated")
for (i in 2:100) {
  lines(post.coef.sum[coef.rng,i])
}
lines(prior.coef.sum[coef.rng],col="red")

# USING THE TRANSPOSED FIELDS
e_num = 76
post = prep_post_time(nc.post, e_num)
prior.coef = proj %*% as.vector(t(prior.ens[,,e_num]))
post.coef = sapply(1:dim(post)[3], function(x) proj %*% as.vector(t(post[,,x])))

coef.rng = 1:1040
# plot(post.coef[coef.rng, 1], type = "l", xlab = "Coefficient", ylab = "Value", main = "transpose")
# for (i in 2:100) {
#   lines(post.coef[coef.rng,i])
# }
# lines(prior.coef[coef.rng],col="red")

# look at the cumulative sum of the coefficients
plot(cumsum(post.coef[coef.rng, 1]), type = "l")
for (i in 2:100) {
  lines(cumsum(post.coef[coef.rng, i]))
}
lines(cumsum(prior.coef[coef.rng]),col="red")

# compute the depths of the regular and the integrated coefficients
prior.coef.sum = cumsum(prior.coef)
post.coef.sum = sapply(1:dim(post)[3], function(x) cumsum(post.coef[,x]))

# look at the fft of the coefficients
prior.coef.fft = fft(prior.coef, inverse = F)
plot(abs(prior.coef.fft), type = "l")

edepth(prior.coef.sum, post.coef.sum)
edepth(prior.coef, post.coef)








# What if we used the SRFS??
# the magnitudes are quite close but the directions are not.
prior.srsf = sign(diff(prior.coef, 1)) * sqrt(abs(diff(prior.coef, 1)))
post.srsf = sign(diff(post.coef, 1)) * sqrt(abs(diff(post.coef, 1)))

prior.srsf = abs(diff(prior.coef, 1))
post.srsf = abs(diff(post.coef, 1))

coef.rng = 1:500
plot(post.srsf[coef.rng, 1], type = "l", xlab = "Coefficient", ylab = "Value")
for (i in 2:100) {
  lines(post.srsf[coef.rng,i])
}
lines(prior.srsf[coef.rng],col="red")




library("fdasrvf")
post.align = align_fPCA(post.coef, 1:1040)

coef.rng = 1:15
plot(post.align$fn[coef.rng, 1], type = "l", xlab = "Coefficient", ylab = "Value", main = "aligned")
for (i in 2:100) {
  lines(post.align$fn[coef.rng,i])
}
lines(prior.coef[coef.rng],col="red")
