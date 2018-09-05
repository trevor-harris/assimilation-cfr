rm(list = ls())

load_power_data = function(data.dir) {
  files = list.files(data.dir)
  
  outdata = vector("list", length(files))

  for (f in 1:length(files)) {
    outdata[[f]] = readRDS(paste0(data.dir, files[f]))
  }
  names(outdata) = gsub(".rds", "", files)
  return(outdata)
}

get_prefix = function(power_data) {
  # inner part removed numbers, outer removed decimals
  unique(gsub('\\..*', '', gsub('[0-9]+', '', names(power_data))))
}

get_pts = function(power_data, prefix) {
  is.pre = grepl(prefix, names(power_data))
  pre_power = names(power_data[is.pre])
  pre_pts = rep(0, length(pre_power))
  for(i in 1:length(pre_power)) {
    pre_pts[i] = as.numeric(unlist(regmatches(pre_power[i],gregexpr("(?>-)*[[:digit:]]+\\.*[[:digit:]]*",pre_power[i], perl=TRUE))))
  }
  return(sort(pre_pts))
}

viz_power = function(ds, prefix, pts, title, xlab) {
  powers = matrix(0, length(pts), 9)
  
  for(i in pts) {
    j = which(pts == i)
    powers[j,] = sapply(1:9, function(x) sum(ds[[paste0(prefix, i)]][x,])/100)
  }
  # plots
  plot(pts, powers[,1], type = "l", main = title, xlab = xlab)
  for(i in 2:9) {
    lines(pts, powers[,i])
  }
  lines(pts, apply(powers, 1, mean), col = "red")
}

viz_power_comp = function(ds, prefix1, prefix2, pts1, pts2, title, xlab, legend.vars) {
  powers1 = matrix(0, length(pts1), 64)
  for(i in pts1) {
    j = which(pts1 == i)
    powers1[j,] = sapply(1:64, function(x) sum(ds[[paste0(prefix1, i)]][x,])/100)
  }
  
  powers2 = matrix(0, length(pts2), 64)
  for(i in pts2) {
    j = which(pts2 == i)
    powers2[j,] = sapply(1:64, function(x) sum(ds[[paste0(prefix2, i)]][x,])/100)
  }
  
  # plots
  plot(pts1, powers1[,1], type = "l", col = "red", main = title, xlab = xlab)
  lines(pts2, powers2[,1], col = "blue")
  for(i in 2:64) {
    lines(pts1, powers1[,i], col = "red")
    lines(pts2, powers2[,i], col = "blue")
  }
  lines(pts1, apply(powers1, 1, mean), col = "darkred", lwd = 3)
  lines(pts2, apply(powers2, 1, mean), col = "darkblue", lwd = 3)
  # legend("right", legend = legend.vars, col=c("red", "blue"), lty=c(1, 1))
}


viz_power_means = function(ds, prefix1, prefix2, pts1, pts2, title, xlab, legend.vars) {
  powers1 = matrix(0, length(pts1), 9)
  for(i in pts1) {
    j = which(pts1 == i)
    powers1[j,] = sapply(1:64, function(x) sum(ds[[paste0(prefix1, i)]][x,])/100)
  }
  
  powers2 = matrix(0, length(pts2), 9)
  for(i in pts2) {
    j = which(pts2 == i)
    powers2[j,] = sapply(1:64, function(x) sum(ds[[paste0(prefix2, i)]][x,])/100)
  }
  
  # plots
  plot(pts1, apply(powers1, 1, mean), type = "l", col = "red", main = title, xlab = xlab)
  lines(pts2, apply(powers2, 1, mean), col = "blue")
  legend("right", legend = legend.vars, col=c("red", "blue"), lty=c(1, 1))
}

data.dir = "/Users/trevh/research/assimilation-cfr/outdata/pointwise/"
list.files(data.dir)

# extract the data components
power_data = load_power_data(data.dir)
prefix = get_prefix(power_data)
pts = lapply(prefix, function(x) get_pts(power_data, x))
names(pts) = prefix

# plot em (9 regions)
viz_power_comp(power_data, prefix[7], prefix[1], pts[[7]], pts[[1]],
               "Depth vs Bonferroni Full", "shift", c("Depth", "Bonferroni"))

viz_power_comp(power_data, prefix[9], prefix[3], pts[[9]], pts[[3]],
               "Depth vs Bonferroni Partial", "shift", c("Depth", "Pointwise"))

viz_power_comp(power_data, prefix[11], prefix[5], pts[[11]], pts[[5]],
               "Depth vs Bonferroni Skew Partial", "shift", c("Depth", "Bonferroni"))

# plot em (64 regions)
viz_power_comp(power_data, prefix[8], prefix[2], pts[[8]], pts[[2]],
               "Depth vs Bonferroni Full", "shift", c("Depth", "Bonferroni"))

viz_power_comp(power_data, prefix[10], prefix[4], pts[[10]], pts[[4]],
               "Depth vs Bonferroni Partial", "shift", c("Depth", "Pointwise"))

viz_power_comp(power_data, prefix[12], prefix[6], pts[[12]], pts[[6]],
               "Depth vs Bonferroni Skew Partial", "shift", c("Depth", "Bonferroni"))





viz_power_means(power_data, prefix[4], prefix[1], pts[[4]], pts[[1]],
               "Depth vs Bonferroni Full", "shift", c("Depth", "Bonferroni"))

viz_power_means(power_data, prefix[5], prefix[2], pts[[5]], pts[[2]],
               "Depth vs Bonferroni Partial", "shift", c("Depth", "Pointwise"))

viz_power_means(power_data, prefix[6], prefix[3], pts[[6]], pts[[3]],
               "Depth vs Bonferroni Skew Partial", "shift", c("Depth", "Bonferroni"))




sum(power_data[["pw_cr_0"]]) / 900
sum(power_data[["bf_cr_0"]]) / 900
sum(power_data[["depth_cr_0"]]) / 900

plot(apply(power_data[["pw_cr_0.4"]], 1, mean), type="l", col = "green", ylim=c(0, 1))
lines(apply(power_data[["bf_cr_0.4"]], 1, mean), col = "blue")
lines(apply(power_data[["depth_cr_0.4"]], 1, mean), col = "red")
 # legend("right", legend = c("Pointwise", "Bonferroni", "Depth"), col=c("green", "red", "blue"), lty=c(1, 1))



# graph em
viz_power(power_data, prefix[2], pts[[2]], 
          "Covariance Power", xlab = "Scale")

viz_power(power_data, prefix[3], pts[[3]], 
          "Constant Mean Power", xlab = "Shift")

viz_power(power_data, prefix[4], pts[[4]], 
          "Constant Mean Power (Eigen)", xlab = "Shift")

viz_power_comp(power_data, prefix[3], prefix[4], pts[[3]], pts[[4]],
               "Random vs Eigen", "shift", c("Random", "Eigen"))

viz_power(power_data, prefix[5], pts[[5]], 
          "Parabolic Mean Power", xlab = "Shift")

viz_power(power_data, prefix[6], pts[[6]], 
          "Partial Mean Power", xlab = "Shift")
