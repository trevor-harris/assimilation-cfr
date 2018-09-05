##### VERTICAL #####

# define data and import fuctions
source('code/setup.R')

##### PREP THE X FUNCTIONS #####
# setup list of fields
x = prior.ens[,,prior.sub]
x.coef = sapply(1:dim(x)[3], function(e) proj %*% as.vector(x[,,e]))

y = prep_post_ens(nc.post, 900)
y.coef = sapply(1:dim(y)[3], function(e) proj %*% as.vector(y[,,e]))

z.coef = cbind(x.coef, y.coef)
z.coef = z.coef - rowMeans(z.coef)

x.coef = z.coef[,1:100]
y.coef = z.coef[,101:200]

plot_ly(showscale = FALSE) %>%
  add_surface(z = ~matrix(x.coef[,2], 40, 26))

coef.rng = 1:40
plot(x.coef[coef.rng, 1], type = "l", xlab = "Coefficient", ylab = "Value", main ="Regular")
for (i in 2:100) {
  lines(x.coef[coef.rng,i])
  # lines(y.coef[coef.rng,i], col = "red")
}
lines(y.coef[,1], col = "red")


