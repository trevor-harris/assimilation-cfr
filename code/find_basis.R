
# AFTER PRIOR HAS BEEN PROCESSED

# document this function omg
basis_coef <- function(vec, basis) {
  vec = vec - mean(vec)
  nbase = ncol(basis)
  coef = rep(0, nbase)
  for (i in 1:nbase) {
    coef[i] = vec %*% basis[,i]
  }
  return(coef)
}

prior.basis = prior.eig$d

ens.basis = matrix(0, length(prior.basis), length(ens))
for(e in 1:length(ens)) {
  ens.basis[,e] = basis_coef(as.vector(ens[[e]]), prior.eig$u)
}

# should match prior.basis
prior.basis = basis_coef(as.vector(prior[[1]]), prior.eig$u)

save(ens.basis, prior.basis, file="data/basis")
save(ens.basis, prior.basis, file="basis")
