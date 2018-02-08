# edepth = function(g, fmat) {
#   obs = nrow(fmat)
#   fns = ncol(fmat)
#   
#   # create the depths
#   g = as.matrix(g)
#   g_depth = matrix(0, nrow(g), ncol(g))
#   for (j in 1:ncol(g)) {
#     g_depth[,j] = depth(g[,j], fmat)
#   }
#   g_depth = as.vector(g_depth)
#   
#   dmat = matrix(0, obs, fns)
#   for (j in 1:fns) {
#     dmat[,j] = depth(fmat[,j], fmat)
#   }
#   
#   # create the dcdfs
#   r = sort(unique(c(g_depth, dmat)))
#   
#   g_dcdf = sapply(r, function(x) sum(g_depth <= x)) / length(g_depth)
#   
#   dcdfmat = matrix(0, length(r), ncol(dmat))
#   for (j in 1:fns) {
#     dcdfmat[,j] = sapply(r, function(x) sum(dmat[,j] <= x))
#   }
#   dcdfmat = dcdfmat / nrow(dmat)
#   
#   # compare the dcdfs
#   ed = 0
#   for (j in 1:fns) {
#     for (i in 1:nrow(dcdfmat)) {
#       if(g_dcdf[i] != dcdfmat[i,j]) {
#         ed = ed + (dcdfmat[i,j] > g_dcdf[i])
#         break
#       }
#     }
#   }
#   return(ed / ncol(dcdfmat))
# }
# 

edepth = function(g, fmat) {
  
  # add g to the functional cloud fmat
  fmat = cbind(g, fmat)
  
  # save the dimensions for convenience
  obs = nrow(fmat)
  fns = ncol(fmat)
  
  # find the depths of g
  g_depth = depth(g, fmat)
  
  # find the depths of each f in fmat
  fmat_depth = matrix(0, obs, fns)
  for (j in 1:fns) {
    fmat_depth[,j] = depth(fmat[,j], fmat)
  }
  
  # get the allowed r values (for calculating the dCDF)
  r = sort(unique(c(g_depth, fmat_depth)))
  
  # find dCDF of g
  g_dcdf = sapply(r, function(x) sum(g_depth <= x))
  
  # find the dCDFs of each f in fmat
  fmat_dcdf = matrix(0, length(r), ncol(fmat_depth))
  for (j in 1:fns) {
    fmat_dcdf[,j] = sapply(r, function(x) sum(fmat_depth[,j] <= x))
  }
  
  # compare the dCDF of g against each dCDF of fmat
  ed = 0
  for (j in 1:fns) {
    for (i in 1:nrow(fmat_dcdf)) {
      if(g_dcdf[i] != fmat_dcdf[i,j]) {
        ed = ed + (fmat_dcdf[i,j] > g_dcdf[i])
        break
      }
    }
  }
  return(ed / fns)
}

