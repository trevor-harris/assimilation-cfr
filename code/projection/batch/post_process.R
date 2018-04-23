rm(list = ls())

data.dir = "/Users/trevh/research/assimilation-cfr/outdata/"
files = list.files(data.dir)

outdata = vector("list", length(files))
for (f in 1:length(files)) {
  load(paste0(data.dir, files[f]))
  outdata[[f]] = diffs
}
names(outdata) = gsub(".RData", "", files)

files

# hard coded COVARIANCE stuff
# over each region
scales = c(0.01, 0.05, 0.1, 0.3, 0.5, 0.7, 1.3, 1.5, 1.7, 2, 3, 4, 5)
cov_powers = matrix(0, length(scales), 9)
for(i in scales) {
  j = which(scales == i)
  cov_powers[j,] = sapply(1:9, function(x) sum(outdata[[paste0("cov_scale_", i)]][x,])/100)
}
# plots
plot(scales, cov_powers[,1], type = "l", main = "Covariance Power")
for(i in 2:13) {
  lines(scales, cov_powers[,i])
}


# hard coded MEAN stuff
# parabolic mean over each region
para_mu = matrix(0, 10, 9)
for(i in 1:10) {
  para_mu[i,] = sapply(1:9, function(x) sum(outdata[[paste0("mu_parab_", i)]][x,])/100)
}
# plots
plot(para_mu[,1], type = "l", main = "Parabolic Mean Power")
for(i in 2:9) {
  lines(para_mu[,i])
}


# Partial mean over each region
partial_mu = matrix(0, 13, 9)
for(i in 1:10) {
  partial_mu[i,] = sapply(1:9, function(x) sum(outdata[[paste0("mu_partial_", i)]][x,])/100)
}
# plots
plot(partial_mu[,1], type = "l", main = "Partial Mean Power")
for(i in 2:13) {
  lines(partial_mu[,i])
}

