rm(list = ls())

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

regions = 64
sims = 100
batches = 10

data.dir = "/Users/trevh/research/assimilation-cfr/simdata/run1/"
list.files(data.dir)
files = list.files(data.dir)

dep = array(0, dim=c(regions, sims, batches))
bonf = array(0, dim=c(regions, sims, batches))

d = 1
b = 1
for(f in files) {
  if(grepl("Depth", f)) {
    dep[,,d] = readRDS(paste0(data.dir, f))
    d = d+1
  }
  if(grepl("Bonferroni", f)) {
    bonf[,,b] = readRDS(paste0(data.dir, f))
    b = b+1
  }
  
  if(grepl("post_mu", f)) {
    post_mu = readRDS(paste0(data.dir, f))
  }
}

# reformat 
d_power = matrix(0, batches, regions)
b_power = matrix(0, batches, regions)
for(k in 1:10) {
  d_power[k,] = apply(dep[,,k], 1, mean)
  b_power[k,] = apply(bonf[,,k], 1, mean)
}

# sort by L2
regions_mu = matsplitter(post_mu, 5, 5)
l2 = sapply(1:regions, function(x) sum(regions_mu[,,x]^2))

# plot
plot(d_power[1,], type = "l", col = "red")
lines(b_power[1,], col = "blue")
for (k in 2:10) {
  lines(d_power[k,], col = "red")
  lines(b_power[k,], col = "blue")
}



l2 = log(l2 + 1)

plot(l2, d_power[1,], col = "red")
points(l2, b_power[1,], col = "blue")
for (k in 2:10) {
  points(l2, d_power[k,], col = "red")
  points(l2, b_power[k,], col = "blue")
}


# boxplot(d_power, use.cols = T)

