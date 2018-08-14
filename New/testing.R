rm(list = ls()); gc()

library(extdepth)
library(tictoc)
library(fda)

#### FUNS ####
create_basis = function(field, basis.pts = 10, r = 1.5) {
  
  x.dim = dim(field)[1]
  y.dim = dim(field)[2]
  
  x.loc = as.integer(seq(1, x.dim, length.out = basis.pts))
  y.loc = as.integer(seq(1, y.dim, length.out = basis.pts))
  
  # find the radius (1.5 times the minimal distance)
  # will either be the distance b/w the first point and the next (horizontally or vertically)
  rad = r * min(abs(x.loc[1] - x.loc[2]), abs(y.loc[1] - y.loc[2]))
  
  n.loc = basis.pts * basis.pts
  grid = expand.grid(x.loc, y.loc)
  basis = matrix(0, x.dim*y.dim, n.loc)
  
  for (i in 1:n.loc) {
    
    # dist = (lat - center.lat)^2 + (lon - center.lon)^2
    dist = c(outer((1:x.dim - grid[i,][[1]])^2, (1:y.dim - grid[i,][[2]])^2, "+"))
    
    # use the bisquare (radial) basis
    basis[,i] = (1 - dist/rad^2)^2 * (sqrt(dist) < rad)
  }
  
  return(basis)
}
plt_2funs = function(f, g) {
  f = as.matrix(f)
  g = as.matrix(g)
  plot(f[,1], type = "l", ylim = c(min(cbind(f, g)), max(cbind(f, g))), col = "red")
  for(i in 2:ncol(f)) {
    lines(f[,i], col = "red")
  }
  for(i in 1:ncol(g)) {
    lines(g[,i], col = "blue")
  }
}
plt_funs = function(fmat) {
  plot(fmat[,1], type = "l", ylim = c(min(fmat), max(fmat)))
  for(i in 2:ncol(fmat)) {
    lines(fmat[,i])
  }
}
sim_gp = function(fields = 100, mu = 0, l = 30, pts = 30) {
  
  grid = 1:pts
  grid = expand.grid(grid, grid)
  distmat = as.matrix(dist(grid))
  
  # calc sigma with cov kernel
  sigma = exp(-distmat / l)
  
  sigma.eig = eigen(sigma)
  sigma.half = sigma.eig$vectors %*% diag(sqrt(sigma.eig$values)) %*% t(sigma.eig$vectors)
  
  
  gps = array(0, dim=c(pts, pts, fields))
  for(f in 1:fields) {
    gps[,,f] = (sigma.half %*% rnorm(pts^2)) + mu
  }
  return(gps)
}

num.int = function(y, x = seq_len(length(y)) / length(y)) {
  sum(diff(x) * (head(y, -1) + tail(y, -1))) / 2
}
l2 = function(f, g) {
  sqrt(num.int((f-g)^2))
}
elastic_depth = function(fmat) {
  # fmat = prior_flat
  
  n = ncol(fmat)
  
  dist_mat = matrix(0, n, n)
  for(i in 1:n) {
    for(j in i:n) {
      dist = l2(fmat[,i], fmat[,j])
      dist_mat[i, j] = dist
      dist_mat[j, i] = dist
    }
  }
  
  1 / (1 + apply(dist_mat, 2, median))
  # depth2 = exp(-apply(dist_mat, 2, median)^2/2)
}

flatten = function(mat) {
  matrix(mat, prod(dim(mat)[1:2]), dim(mat)[3])
}

# distance
ks.dist = function(x, y) {
  nx = ncol(x)
  ny = ncol(y)
  
  x_dist = rep(0, nx)
  y_dist = rep(0, ny)
  
  x_ed = edepth_set(x)
  x_med = x[,which.max(x_ed)]
  
  for(f in 1:nx) {
    x_dist[f] = l2(x[,f], x_med)
  }
  for(f in 1:ny) {
    y_dist[f] = l2(y[,f], x_med)
  }
  
  x_dist = x_dist[x_dist > 0]
  
  ks.test(x_dist, y_dist)$p.value
  # wilcox.test(x_dist, y_dist)$p.value
  # t.test(x_dist, y_dist)$p.value
}

# ED
depth <- function(g, fmat) {
  
  # Computes the depth values of a function with respect to a set of functions (fmat)
  fn = ncol(fmat)
  depth = rep(0, length(g))
  
  for (row in 1:nrow(fmat)) {
    diff = abs(sum(sign(g[row] - fmat[row,])))
    depth[row] = 1 - (diff / fn)
  }
  
  return(depth)
}
edepth_multi = function(gmat, fmat_sorted) {
  fmat = fmat_sorted
  fdepths = depth_set(fmat)
  gdepths = apply(gmat, 2, function(x) depth(x, fmat))
  
  # get the allowed r values (for calculating the dCDF)
  rvals = sort(unique(c(gdepths, fdepths)))
  
  ged = rep(1, ncol(gmat))
  for(g in 1:ncol(gmat)) {
    for(f in 1:ncol(fmat)) {
      for(r in rvals) {
        dg = sum(gdepths[,g] <= r)
        df = sum(fdepths[,f] <= r)
        if(dg != df) {
          break;
        }
      }
      if (dg > df) {
        ged[g] = (f-1) / ncol(gmat)
        break;
      }
    }
  }
  ged
}
ks_pval = function(t, n = 10) {
  2*(sum(sapply(1:n, function(x) (-1)^(x-1) * exp(-2*(x^2)*t))))
}
ks.depth = function(f, g) {
  
  f = edepth_sort(f)
  fed = (1:ncol(f)) / ncol(f)
  ged = edepth_multi(g, f)
  
  f.surv = rev(c(0, sort(fed)))
  g.cdf = sapply(1:length(f.surv), function(x) mean(ged >= f.surv[x]))
  f.cdf = rev(f.surv)
  
  ks = max(abs(f.cdf - g.cdf))
  ks_pval(sqrt(ncol(g))*ks)
}

# integrated ED
depth = function(g, fmat) {
  
  # Computes the depth values of a function with respect to a set of functions (fmat)
  fn = ncol(fmat)
  depth = rep(0, length(g))
  
  for (row in 1:nrow(fmat)) {
    diff = abs(sum(sign(g[row] - fmat[row,])))
    depth[row] = 1 - (diff / fn)
  }
  
  return(depth)
}
cdf = function(x) {
  sapply(1:length(x), function(z) mean(x <= z/length(x)))
}
tdepth = function(gmat, fmat) {
  gdep = apply(gmat, 2, function(x) depth(x, fmat))
  gcdf = apply(gdep, 2, function(x) cdf(x))
}
ks.int = function(x, y) {
  dx = tdepth(x, x)
  dy = tdepth(y, x)

  nx = ncol(x)
  ny = ncol(y)
  
  x_dist = rep(0, nx)
  y_dist = rep(0, ny)
  
  for(f in 1:nx) {
    x_dist[f] = l2(dx[,f], 0)
  }
  for(f in 1:ny) {
    y_dist[f] = l2(dy[,f], 0)
  }
  
  # ks.test(x_dist, y_dist)$p.value
  # wilcox.test(x_dist, y_dist)$p.value
  t.test(x_dist, y_dist)$p.value
}

# elastic
elastic_depth = function(fmat) {
  # fmat = prior_flat
  
  n = ncol(fmat)
  
  dist_mat = matrix(0, n, n)
  for(i in 1:n) {
    for(j in i:n) {
      dist = l2(fmat[,i], fmat[,j])
      dist_mat[i, j] = dist
      dist_mat[j, i] = dist
    }
  }
  
  1 / (1 + apply(dist_mat, 2, median))
  # depth2 = exp(-apply(dist_mat, 2, median)^2/2)
}
elastic_depth_multi = function(gmat, fmat) {
  ng = ncol(gmat)
  nf = ncol(fmat)
  
  dist_mat = matrix(0, ng, nf)
  for(i in 1:ng) {
    for(j in 1:nf) {
      dist = l2(gmat[,i], fmat[,j])
      dist_mat[i, j] = dist
      dist_mat[j, i] = dist
    }
  }
  
  dist_mat = cbind(rep(0, nrow(dist_mat)), dist_mat)
  
  1 / (1 + apply(dist_mat, 1, median))
  # depth2 = exp(-apply(dist_mat, 2, median)^2/2)
}
elas.depth = function(f, g) {
  
  f = prior_flat
  g = post_flat
  
  fed = elastic_depth(f)
  ged = elastic_depth_multi(g, f)
  
  f.surv = rev(c(0, sort(fed)))
  g.cdf = sapply(1:length(f.surv), function(x) mean(ged >= f.surv[x]))
  f.cdf = sapply(1:length(f.surv), function(x) mean(fed >= f.surv[x]))
  
  ks = max(abs(f.cdf - g.cdf))
  ks_pval(sqrt(ncol(g))*ks)
}

# median depth
mdepth = function(gmat, fmat) {
  apply(gmat, 2, function(x) median(depth(x, fmat)))
}
ks.med = function(x, y) {
  mx = mdepth(x, x)
  my = mdepth(y, x)
  
  ks.test(mx, my)$p.value
}


# get args
batch = 30
sims = 40
set.seed(batch + 1023)

kdist = rep(0, sims)
kdept = rep(0, sims)
kinte = rep(0, sims)
kelas = rep(0, sims)
kmedi = rep(0, sims)

for(s in 1:sims) {
  tic("Total")
  cat("Simulation ", s, "\n")
  
  ### Generate Fields
  tic("Generating fields")
  
  ## generate
  nfield = 200
  fields = sim_gp(nfield, mu = 0, l = 45)
  
  ## smoothe
  basis = create_basis(matrix(0, 30, 30), basis.pts = 10)
  proj = solve(t(basis) %*% basis) %*% t(basis)
  
  fields = vapply(1:nfield, function(x) basis %*% (proj %*% as.vector(fields[,,x])),
                  FUN.VALUE = matrix(0, 30, 30))
  
  ## demean
  # fields = vapply(1:nfield, function(x) fields[,,x] - mean(fields[,,x]), FUN.VALUE = matrix(0, 30, 30))
  
  ## split
  ind = sample(1:nfield, nfield/2, replace = F)
  prior = fields[,,ind]
  posterior = fields[,,-ind]
  toc()
  
  ### Pointwise Comparisons
  tic("ED")
  prior_flat = flatten(prior)
  post_flat = flatten(posterior)

  kdist[s] = ks.dist(prior_flat, post_flat)
  kdept[s] = ks.depth(prior_flat, post_flat)
  kinte[s] = ks.int(prior_flat, post_flat)
  kelas[s] = elas.depth(prior_flat, post_flat)
  kmedi[s] = ks.med(post_flat, prior_flat)
  
  toc()
  cat("\n")
}

mean(kdist < 0.05)
mean(kdept < 0.05)
mean(kinte < 0.05)
mean(kelas < 0.05)
mean(kmedi < 0.05)

plot(kdist)
plot(kdept)
plot(kinte)
plot(kelas)
plot(kmedi)


