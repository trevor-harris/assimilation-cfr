permute_fields = function(prior, posterior) {
  # permutes 2 samples of 2D regionalized functions 
  # prior.split = 4D array. dim1 = lat, dim2 = lon, dim3 = region no, dim4 = sample number
  # post.split = same as prior.split
  # returns a list = [permuted.prior, permuted.post]
  
  nlat = dim(prior)[1]
  nlon = dim(prior)[2]
  ens = dim(prior)[3]
  
  prior.new = prior
  post.new = posterior
  
  prior.ind = sample(c(0, 1), size = ens, replace = TRUE)
  
  for (e in 1:ens) {
    if (prior.ind[e] == 1) {
      prior.new[,,e] = posterior[,,e]
      post.new[,,e] = prior[,,e]
    } else {
      prior.new[,,e] = prior[,,e]
      post.new[,,e] = posterior[,,e]
    }
  }
  
  return(list(prior.new, post.new))
}

# reformat 
flatten = function(mat) {
  matrix(mat, prod(dim(mat)[1:2]), dim(mat)[3])
}
