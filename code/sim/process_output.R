rm(list = ls())

regions = 64
sims = 100
batches = 10

data.dir = "/Users/trevh/research/assimilation-cfr/simdata/simdata/"
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
regions_l2 = sapply(1:regions, function(x) sum(regions_mu[,,x]^2))
regions_order = order(regions_l2)

# plot
plot(d_power[1,], type = "l", col = "red")
lines(b_power[1,], col = "blue")
for (k in 2:10) {
  lines(d_power[k,], col = "red")
  lines(b_power[k,], col = "blue")
}


