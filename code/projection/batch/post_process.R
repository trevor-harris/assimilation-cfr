rm(list = ls())

data.dir = "/Users/trevh/research/assimilation-cfr/outdata/"

load_power_data = function(data.dir) {
  files = list.files(data.dir)
  
  outdata = vector("list", length(files))
  for (f in 1:length(files)) {
    load(paste0(data.dir, files[f]))
    outdata[[f]] = diffs
  }
  names(outdata) = gsub(".RData", "", files)
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


# extract the data components
power_data = load_power_data(data.dir)
prefix = get_prefix(power_data)
pts = sapply(prefix, function(x) get_pts(power_data, x))

# graph em
viz_power(power_data, prefix[1], pts[[1]], 
          "Covariance Power", xlab = "Scale")

viz_power(power_data, prefix[2], pts[[2]], 
          "Constant Mean Power", xlab = "Shift")

viz_power(power_data, prefix[3], pts[[3]], 
          "Parabolic Mean Power", xlab = "Shift")

viz_power(power_data, prefix[4], pts[[4]], 
          "Partial Mean Power", xlab = "Shift")
