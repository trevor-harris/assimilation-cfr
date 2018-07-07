


betas = matrix(0, 30, 30)
for(i in 1:30) {
  for(j in 1:30) {
    betas[i,j] = mean(posterior[i,j,]-prior[i,j,])
    # betas[i,j] = t.test(prior[i,j,], posterior[i,j,], paired = T)$statistic
  }
}

permuted_betas = array(0, c(30, 30, perms))

for(p in 1:perms) {
  permuted = permute_fields(prior, posterior)
  prior.p = permuted[[1]]
  posterior.p = permuted[[2]]
  
  betas.p = matrix(0, 30, 30)
  for(i in 1:30) {
    for(j in 1:30) {
      betas.p[i,j] = mean(posterior.p[i,j,]-prior.p[i,j,])
      # betas[i,j] = t.test(prior.p[i,j,], posterior.p[i,j,], paired = T)$statistic
    }
  }
  
  permuted_betas[,,p] = betas.p
}

# observed betas
beta_obs[,,s] = betas

# ED based upper and lower
betas_flat = flatten(permuted_betas)
betas_ed = edepth_set(betas_flat)
betas_cr = central_region(betas_flat, betas_ed)

ed_lower[,,s] = matrix(betas_cr[[1]], 30, 30)
ed_upper[,,s] = matrix(betas_cr[[2]], 30, 30)