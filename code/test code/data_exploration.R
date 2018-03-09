library(extdepth)
load("basis")

Rcpp::sourceCpp("/Users/trevh/research/extdepth/src/depth.cpp")

# exdepth works so long as g is included in fmat. Does that make sense? No
full_basis = cbind(prior.basis, ens.basis)
Rcpp::sourceCpp("/Users/trevh/research/extdepth/src/depthtest.cpp")


extdepth(prior.basis, cbind(prior.basis, ens.basis))
edepth(prior.basis, ens.basis)

edepth(ens.basis[,3], cbind(prior.basis, ens.basis))

extdepth(ens.basis[,2], cbind(prior.basis, ens.basis))
edepth(ens.basis[,2], cbind(prior.basis, ens.basis))

edepth(cbind(prior.basis, prior.basis), cbind(prior.basis, ens.basis))

# single comparison
plot(prior.basis, type = "l")
lines(ens.basis[,1])


# what do the ens coeff look like
plot(ens.basis[,1], type = "l")
for (i in 2:100) {
  lines(ens.basis[,i])
}
# plus prior
lines(prior.basis, col = "red")

# same comparison except on dCDFs
plot(sapply(1:100, function(x) sum(depth(prior.basis, ens.basis) <= x/100))
     ,type = "l"
     ,col = "red"
     ,ylim=c(0, 5)
     ,xlim=c(1,2))


for (i in 1:100) {
  lines(sapply(1:100, function(x) sum(depth(ens.basis[,i], ens.basis) <= x/100)))
}
# lines(sapply(1:100, function(x) sum(depth(prior.basis, ens.basis) <= x/100)), col = "red")
