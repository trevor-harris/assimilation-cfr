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

library(abind)
eigen.proj = function(prior.split, post.split, eig=5) {
  
  nlat = dim(prior.split)[1]
  nlon = dim(prior.split)[2]
  regions = dim(prior.split)[3]
  ens = dim(prior.split)[4]
  
  all.split = array(abind(prior.split, post.split, along=4), dim = c(nlat, nlon, regions, (2*ens)))
  
  all.eigen = array(0, dim = c(nlat*nlon, eig, regions))
  for (r in 1:regions) {
    all.split.vec = t(sapply(1:(2*ens), function(x) as.vector(all.split[,,r,x])))
    all.eigen[,,r] = princomp(all.split.vec)$scores[1:eig]
  }
  return(all.eigen)
}


kst.field.eig = function(prior.split, post.split, eigen=5, proj_mat = NULL) {
  # computes the KST field statistic between two samples of X and Y regionalized 2D functions
  # prior.split and post.split same as in permute_fields
  
  nlat = dim(prior.split)[1]
  nlon = dim(prior.split)[2]
  regions = dim(prior.split)[3]
  ens = dim(prior.split)[4]
  
  if (is.null(proj_mat)) {
    proj_mat = eigen.proj(prior.split, post.split, eigen)
  } else {
    proj_mat = proj_mat
    eigen = dim(proj_mat)[2]
  }
  
  prior.proj = matrix(0, ens, eigen)
  post.proj = matrix(0, ens, eigen)
  kol.field = 1:regions
  
  for(r in 1:regions) {
    for(e in 1:ens) {
      prior.proj[e,] = as.vector(as.vector(prior.split[,,r,e]) %*% proj_mat[,,r])
      post.proj[e,] = as.vector(as.vector(post.split[,,r,e]) %*% proj_mat[,,r])
    }
    kol = sapply(1:eigen, function(x) kst.fast(prior.proj[,x], post.proj[,x]))
    kol.field[r] = mean(kol)
  }
  kol.field
}


kst.permute.eig = function(prior.split, post.split, eigen = 5, perms = 100) {
  # generates the permutation distribution for the kst field.
  
  regions = dim(prior.split)[3]
  perm.fields = matrix(0, regions, perms)
  
  proj_mat = eigen.proj(prior.split, post.split, eigen)
  
  # dont just naively use the field test. Eigenfunctions are INVARIANT to permutation
  for (p in 1:perms) {
    new.fields = permute_fields(prior.split, post.split)
    prior = new.fields[[1]]
    post = new.fields[[2]]
    perm.fields[,p] = kst.field.eig(prior, post, proj_mat = proj_mat)
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





l2.field = function(prior.split, post.split) {
  # computes the KST field statistic between two samples of X and Y regionalized 2D functions
  # prior.split and post.split same as in permute_fields

  nlat = dim(prior.split)[1]
  nlon = dim(prior.split)[2]
  regions = dim(prior.split)[3]
  ens = dim(prior.split)[4]

  prior.proj = rep(0, ens)
  post.proj = rep(0, ens)
  kol.field = 1:regions

  for(r in 1:regions) {
    for(e in 1:ens) {
      # prior.proj[e] = sum((prior.split[,,r,e])^2)
      # post.proj[e] = sum((post.split[,,r,e])^2)

      prior.proj[e] = max((prior.split[,,r,e]))
      post.proj[e] = max((post.split[,,r,e]))
    }
    kol.field[r] = kst.fast(prior.proj, post.proj)
  }
  kol.field
}

l2.permute = function(prior.split, post.split, perms = 100) {
  # generates the permutation distribution for the kst field.
  regions = dim(prior.split)[3]
  perm.fields = matrix(0, regions, perms)

  # dont just naively use the field test. Eigenfunctions are INVARIANT to permutation
  for (p in 1:perms) {
    new.fields = permute_fields(prior.split, post.split)
    prior = new.fields[[1]]
    post = new.fields[[2]]
    perm.fields[,p] = l2.field(prior, post)
  }

  return(perm.fields)
}




























