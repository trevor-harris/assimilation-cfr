library(fdasrvf)
library(extdepth)


#### Some useful functions
num.int = function(y, x = seq_len(length(y)) / length(y)) {
  sum(diff(x) * (head(y, -1) + tail(y, -1))) / 2
}
exp_map = function(v, mu) {
  norm.v = sqrt(num.int(v^2))
  cos(norm.v)*mu + sin(norm.v) * v / norm.v
} 
inv_psi = function(psi) {
  x = seq_len(length(psi)) / length(psi)
  sapply(1:length(psi), function(i) num.int(psi[1:i]^2, x[1:i]))
}


# generate some gamma functions5
gams = t(rgam(100, 5, 100))

# convert to psi and tangent space
gams.srsf = SqrtMean(gams)

# convenience
gams.psi = gams.srsf$psi
gams.vec = gams.srsf$vec

# proof the functions work
plot(inv_psi(exp_map(gams.vec[,2], gams.srsf$mu)), col = "red")
lines(gams[,2])

# compute ED on tangent vectors and the original functions. Compare.
vec.ed = edepth_set(gams.vec)
gam.ed = edepth_set(gams)


# plot tangent functions.
plot(gams.vec[,1], type = "l", col = "grey", ylim = c(min(gams.vec), max(gams.vec)))
for (i in 2:ncol(gams)) {
  lines(gams.vec[,i], col = "grey")
}
# add the cr
vec.cr = central_region(gams.vec, vec.ed, alpha = 0.5)
lines(gams.vec[,which(vec.ed == 1)], col = "red")
lines(vec.cr$upper, col = "red")
lines(vec.cr$lower, col = "red")

# convert back to gamma functions
gam.lower = inv_psi(exp_map(vec.cr$lower, gams.srsf$mu))
gam.upper = inv_psi(exp_map(vec.cr$upper, gams.srsf$mu))

# plot original functions. Add tangent median and CR (red) and original median (blue)
plot(gams[,1], type = "l", col = "grey")
for (i in 2:ncol(gams)) {
  lines(gams[,i], col = "grey")
}
# tangent based ED functions 
lines(gams[,which(vec.ed == 1)], col = "red", lwd = 2)
lines(gam.upper, col = "red", lwd = 1.5)
lines(gam.lower, col = "red", lwd = 1.5)

# original ED functions
gam.cr = central_region(gams, gam.ed, alpha = 0.5)
lines(gams[,which(gam.ed == 1)], col = "blue", lwd = 2)
lines(gam.cr$upper, col = "blue", lwd = 1.5)
lines(gam.cr$lower, col = "blue", lwd = 1.5)