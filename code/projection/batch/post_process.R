rm(list = ls())

data.dir = "/Users/trevh/research/assimilation-cfr/outdata/"
files = list.files(data.dir)

outdata = vector("list", length(files))
for (f in 1:length(files)) {
  load(paste0(data.dir, files[f]))
  outdata[[f]] = diffs
}
names(outdata) = gsub(".RData", "", substring(files, 6))



# hard coded COVARIANCE stuff
# over each region
scales = c(1, 3, 5, 7, 13, 15, 17, 20, 30, 40, 50)
cov_powers = matrix(0, length(scales), 9)
for(i in scales) {
  j = which(scales == i)
  cov_powers[j,] = sapply(1:9, function(x) sum(outdata[[paste0("cov_scale_", i)]][x,])/100)
}
# plots
plot(scales/10, cov_powers[,1], type = "l", main = "Covariance Power")
for(i in 2:9) {
  lines(scales/10, cov_powers[,i])
}


# hard coded MEAN stuff
# over each region
mu_powers = matrix(0, 10, 9)
for(i in 1:10) {
  mu_powers[i,] = sapply(1:9, function(x) sum(outdata[[paste0("mu_parab_", i)]][x,])/100)
}
# plots
plot(mu_powers[,1], type = "l", main = "Parabolic Mean Power")
for(i in 2:9) {
  lines(mu_powers[,i])
}

# specific one
sapply(1:9, function(x) sum(outdata[["mu_parab_5"]][x,])/100)

