#### permutation test ####
permute_fields = function(prior.split, post.split) {
  # permutes 2 samples of 2D regionalized functions 
  # prior.split = 4D array. dim1 = lat, dim2 = lon, dim3 = region no, dim4 = sample number
  # post.split = same as prior.split
  # returns a list = [permuted.prior, permuted.post]
  
  nlat = dim(prior.split)[1]
  nlon = dim(prior.split)[2]
  regions = dim(prior.split)[3]
  ens = dim(prior.split)[4]
  
  prior.split.new = prior.split
  post.split.new = post.split
  
  prior.ind = sample(c(0, 1), size = ens, replace = TRUE)
  
  for (e in 1:ens) {
    if (prior.ind[e] == 1) {
      prior.split.new[,,,e] = post.split[,,,e]
      post.split.new[,,,e] = prior.split[,,,e]
    } else {
      prior.split.new[,,,e] = prior.split[,,,e]
      post.split.new[,,,e] = post.split[,,,e]
    }
  }
  
  return(list(prior.split.new, post.split.new))
}

#### KST ####
kst.fast = function(x, y) {
  # computes the Kolmogorov-Smirnov 2 sample test between X and Y
  # no p values, just the statistic

  n <- length(x)
  n.x <- as.double(n)
  n.y <- length(y)
  w <- c(x, y)
  z <- cumsum(ifelse(order(w) <= n.x, 1/n.x, -1/n.y))
  max(abs(z))
}

kst.field = function(prior.split, post.split, iter=100) {
  # computes the KST field statistic between two samples of X and Y regionalized 2D functions
  # prior.split and post.split same as in permute_fields
  
  nlat = dim(prior.split)[1]
  nlon = dim(prior.split)[2]
  regions = dim(prior.split)[3]
  ens = dim(prior.split)[4]
  
  prior.proj = matrix(0, ens, iter)
  post.proj = matrix(0, ens, iter)
  kol.field = 1:regions
  proj = matrix(rnorm(nlat*nlon*iter), nlat*nlon, iter)
  
  for(r in 1:regions) {
    for(e in 1:ens) {
      prior.proj[e,] = as.vector(as.vector(prior.split[,,r,e]) %*% proj)
      post.proj[e,] = as.vector(as.vector(post.split[,,r,e]) %*% proj)
    }
    kol = sapply(1:iter, function(x) kst.fast(prior.proj[,x], post.proj[,x]))
    kol.field[r] = mean(kol)
  }
  kol.field
}

kst.permute2 = function(prior.split, post.split, perms = 100, proj = 100) {
  
  # get dim
  nlat = dim(prior.split)[1]
  nlon = dim(prior.split)[2]
  regions = dim(prior.split)[3]
  ens = dim(prior.split)[4]
  
  # setup intermediate variables
  prior.proj = matrix(0, ens, proj)
  post.proj = matrix(0, ens, proj)
  kol.field = 1:regions
  proj.mat = matrix(rnorm(nlat*nlon*proj), nlat*nlon, proj)
  
  # setup output
  observed.field = rep(0, regions)
  permuted.fields = matrix(0, regions, perms)
  
  for(r in 1:regions) {
    for(e in 1:ens) {
      prior.proj[e,] = as.vector(as.vector(prior.split[,,r,e]) %*% proj.mat)
      post.proj[e,] = as.vector(as.vector(post.split[,,r,e]) %*% proj.mat)
    }
    # first find the observed KST statistic
    kst = sapply(1:proj, function(x) kst.fast(prior.proj[,x], post.proj[,x]))
    observed.field[r] = mean(kst)
    
    # Second find the permutation distribution
    for(p in 1:perms) {
      swap = matrix(rbinom(ens*proj, 1, 0.5), ens, proj)
      prior.perm = prior.proj * (1 - swap) + post.proj * swap
      post.perm = prior.proj * swap + post.proj * (1 - swap)
      
      kst = sapply(1:proj, function(x) kst.fast(prior.perm[,x], post.perm[,x]))
      permuted.fields[r,p] = mean(kst)
    }
  }
  return(list(obs = observed.field,
              perm = permuted.fields))
}

kst.permute = function(prior.split, post.split, perms = 100, iter = 100) {
  # generates the permutation distribution for the kst field.
  
  regions = dim(prior.split)[3]
  perm.fields = matrix(0, regions, perms)
  
  for (p in 1:perms) {
    new.fields = permute_fields(prior.split, post.split)
    prior = new.fields[[1]]
    post = new.fields[[2]]
    perm.fields[,p] = kst.field(prior, post, iter = iter)
  }
  
  return(perm.fields)
}


#### regionalize ####
matsplitter<-function(M, r, c) {
  # splits 1 matrix into c MxR matricies
  # I have no idea how this works
  
  rg <- (row(M)-1)%/%r+1
  cg <- (col(M)-1)%/%c+1
  rci <- (rg-1)*max(cg) + cg
  N <- prod(dim(M))/r/c
  cv <- unlist(lapply(1:N, function(x) M[rci==x]))
  dim(cv)<-c(r,c,N)
  cv
}





























