regions = 64
sims = 100
batches = 11

data.dir = "/Users/trevh/research/assimilation-cfr/outdata/pointwise/"
list.files(data.dir)
files = list.files(data.dir)

dep = array(0, dim=c(regions, sims, batches))
bonf = array(0, dim=c(regions, sims, batches))

d = 1
b = 1
for(f in files) {
  if(grepl("depth_skew_cr_high", f)) {
    dep[,,d] = readRDS(paste0(data.dir, f))
    d = d+1
  }
  if(grepl("bf_skew_cr_high", f)) {
    bonf[,,b] = readRDS(paste0(data.dir, f))
    b = b+1
  }
}

# reformat 
d_power = matrix(0, batches, regions)
b_power = matrix(0, batches, regions)
power_order = c(10, 1, 2, 3, 4, 5, 6, 7 ,8, 9, 11)
for(k in power_order) {
  d_power[which(power_order == k),] = apply(dep[,,k], 1, mean)
  b_power[which(power_order == k),] = apply(bonf[,,k], 1, mean)
}


# plot
plot(d_power[,1], type = "l", col = "red")
lines(b_power[,1], col = "blue")
for (k in 2:regions) {
  lines(d_power[,k], col = "red")
  lines(b_power[,k], col = "blue")
}


# high power group
high_d = d_power[, which(d_power[6,] < 1 & d_power[6,] > 0.8)]
high_b = b_power[, which(b_power[6,] < 1 & b_power[6,] > 0.8)]

# medium power group
medium_d = d_power[, which(d_power[6,] < 0.8 & d_power[6,] > 0.4)]
medium_b = b_power[, which(b_power[6,] < 0.8 & b_power[6,] > 0.4)]

# lower power group
low_d = d_power[, which(d_power[6,] < 0.4 & d_power[6,] > 0.2)]
low_b = b_power[, which(b_power[6,] < 0.4 & b_power[6,] > 0.2)]

# lines for each grouph
plot(apply(high_d, 1, mean), type = "l", col = "red")
lines(apply(high_b, 1, mean), col = "blue")

lines(apply(medium_d, 1, mean), col = "red")
lines(apply(medium_b, 1, mean), col = "blue")

lines(apply(low_d, 1, mean), col = "red")
lines(apply(low_b, 1, mean), col = "blue")

