rm(list = ls())
gc()

library(extdepth)
library(ncdf4)
library(dplyr)
library(reshape2)
library(ggplot2)
library(OpenImageR)

# prepare data
prep_prior = function(nc.prior) {
  
  n.lon = nc.prior$dim$lon$len
  n.lat = nc.prior$dim$lat$len
  n.ens = nc.prior$dim$time2$len
  
  # extract data from the ncdf4 objects
  prior = ncvar_get(nc.prior, attributes(nc.prior$var)$names[1], start = c(1, 1, 1), count = c(-1, -1, -1))
  
  # transpose for intuitive (to me) layout
  prior = aperm(prior, c(2, 1, 3))
  
  # remove lat means
  # prior = vapply(1:n.ens, function(x) prior[,,x] - rowMeans(prior[,,x]), FUN.VALUE = matrix(0, nrow = n.lat, ncol = n.lon))
  
  # normalize
  lats = as.vector(nc.prior$dim$lat$vals)
  latmat = matrix(rep(lats, n.lon), n.lat, n.lon)
  latmat = sqrt(abs(cos(latmat*pi/180)))
  
  prior = vapply(1:n.ens, function(x) prior[,,x]*latmat, FUN.VALUE = matrix(0, nrow = n.lat, ncol = n.lon))
  
  return(prior)
}
prep_post = function(nc.post, t) {
  
  n.lon = nc.post$dim$lon$len
  n.lat = nc.post$dim$lat$len
  n.ens = nc.post$dim$sub_ens$len
  
  # extract data from the ncdf4 objects
  ens = ncvar_get(nc.post, attributes(nc.post$var)$names[1], start = c(1, 1, t, 1), count = c(-1, -1, 1, -1))
  
  # transpose for intuitive (to me) layout
  ens = aperm(ens, c(2, 1, 3))
  
  # remove lat means
  # ens = vapply(1:n.ens, function(x) ens[,,x] - rowMeans(ens[,,x]), FUN.VALUE = matrix(0, nrow = n.lat, ncol = n.lon))
  
  # normalize
  lats = as.vector(nc.post$dim$lat$vals)
  latmat = matrix(rep(lats, n.lon), n.lat, n.lon)
  latmat = sqrt(abs(cos(latmat*pi/180)))
  
  ens = vapply(1:n.ens, function(x) ens[,,x]*latmat, FUN.VALUE = matrix(0, nrow = n.lat, ncol = n.lon))
  
  return(ens)
}
flatten = function(mat) {
  matrix(mat, prod(dim(mat)[1:2]), dim(mat)[3])
}

# depth CDF
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
meandepth = function(gmat, fmat) {
  apply(gmat, 2, function(x) mean(depth(x, fmat)))
}
ks.mean = function(f, g) {
  fed = meandepth(f, f)
  ged = meandepth(g, f)
  
  f.surv = rev(c(0, sort(fed)))
  g.cdf = sapply(1:length(f.surv), function(x) mean(ged > f.surv[x]))
  f.cdf = sapply(1:length(f.surv), function(x) mean(fed > f.surv[x]))
  
  ks = max(abs(f.cdf - g.cdf))
  rate = sqrt((ncol(g)*ncol(f)) / (ncol(g) + ncol(f)))
  ks_pval(rate*ks)
}
ks_pval = function(t, n = 20) {
  2*(sum(sapply(1:n, function(x) (-1)^(x-1) * exp(-2*(x^2)*t^2))))
}

# plots
remove_cr = function(cr, gmat, downsamp=1) {
  lower = cr$lower
  upper = cr$upper
  out = rowMeans((gmat - lower)*(lower > gmat) + (gmat - upper)*(upper < gmat))

  matrix(out, 96/downsamp, 144/downsamp)
}
remove_cr2 = function(cr, gmat, downsamp=1) {
  gmat = post
  lower = cr$lower
  upper = cr$upper
  out = (gmat - lower)*(lower > gmat) + (gmat - upper)*(upper < gmat)
  out[out == 0] = NA
  out = rowMeans(out, na.rm = TRUE)
  out[is.nan(out)] = 0
  matrix(out, 96/downsamp, 144/downsamp)
}
field_plot <- function(field, nc, main = "", downsamp = 1, zlim = c(-max(abs(field)), max(abs(field)))) {
  
  lats = as.vector(nc$dim$lat$vals)[seq(1, 96, by=downsamp)]
  lons = as.vector(nc$dim$lon$vals)[seq(1, 144, by=downsamp)]
  dimnames(field) = list(lats, ifelse(lons >= 180, lons - 360, lons))
  
  field.gg = melt(field)
  colnames(field.gg) = c("lat", "lon", "value")
  
  world = map_data("world")
  world = world[world$long <= 178, ]
  
  zlim = c(-max(abs(field)), max(abs(field)))
  
  ggplot() +
    geom_raster(data = field.gg, aes(x=lon, y=lat, fill=value), interpolate = TRUE) +
    geom_polygon(data = world, aes(x=long, y=lat, group=group), fill = NA, color="black") +
    coord_cartesian() +
    scale_fill_gradient2(midpoint=0, low="blue", mid="white", high="red") +
    # scale_fill_gradient2(midpoint=0, low="#4393c3", mid="white", high="#d6604d") +
    # scale_fill_gradient2(midpoint=0, low="#2166ac", mid="white", high="#b2182b") +
    theme_void() +
    ggtitle(main) +
    theme(plot.title = element_text(hjust = 0.5))
}
save_dir = "../Presentation/results/"

# pvals = read.csv("../research/assimilation-cfr/cfr/adjusted_pvals")$x
pvals = read.csv("adjusted_pvals")$x

signif = as.numeric(pvals < 0.05)
time = 1:998
signif.model = glm(signif ~ time, family="binomial")

sig.logit = data.frame(time = time, sig = signif)
ggplot(sig.logit, aes(x = time, y = sig)) +
  geom_point(alpha = 0.2) +
  geom_smooth(method = "glm", 
              method.args = list(family = "binomial"), 
              se = TRUE) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5)) +
  ylab("Significant") +
  ggtitle("Significance over time")
ggsave(paste0(save_dir, "sig_time.png"), width = 5, height = 3.2)


sig.exp = data.frame(time = time, pexp = -log10(pvals))
ggplot(sig.exp, aes(x = time, y = pexp)) +
  geom_point(alpha = 0.2) +
  # geom_smooth() +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5)) +
  ylab("-log(pvalue)") +
  ggtitle("-log(p-values) over time")
ggsave(paste0(save_dir, "pval_time.png"), width = 5, height = 3.2)

# inspect some of the times with the bigggest differences
# nc.post = nc_open('../research/climate_data/tas_ens_da_hydro_r.1000-2000_d.16-Feb-2018.nc')
# nc.prior = nc_open('../research/climate_data/tas_prior_da_hydro_r.1000-2000_d.16-Feb-2018.nc')

nc.post = nc_open('../data/tas_ens_da_hydro_r.1000-2000_d.16-Feb-2018.nc')
nc.prior = nc_open('../data/tas_prior_da_hydro_r.1000-2000_d.16-Feb-2018.nc')

downsamp = 1

prior = prep_prior(nc.prior)
prior = sapply(1:dim(prior)[3], function(x) down_sample_image(prior[,,x], downsamp))

prior.depths = meandepth(prior, prior)

prior.surv = rev(c(0, sort(prior.depths)))
prior.cdf = sapply(1:length(prior.surv), function(x) mean(prior.depths > prior.surv[x]))

##### CDF plots #####
times = c(1, 300, 600, 998)
for(t in times) {
  post = prep_post(nc.post, t)
  post = sapply(1:dim(post)[3], function(x) down_sample_image(post[,,x], downsamp))
  
  post.depths = meandepth(post, prior)
  post.cdf = sapply(1:length(prior.surv), function(x) mean(post.depths > prior.surv[x]))
  
  cdfgg = data.frame(prior = prior.cdf, posterior = post.cdf)
  cdfgg = melt(cdfgg)
  cdfgg["depth"] = rep(seq(0, 1, length.out = 999), 2)
  
  pval.t = formatC(pvals[t], format = "e", digits = 1)
  
  plt.t = ggplot(cdfgg) +
    geom_line(
      aes(
        x = depth,
        y = value,
        color = variable),
      size = 0.9) +
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5)) +
    ggtitle(paste0("Year ", t, " (p = ", pval.t, ")")) +
    xlab("alpha") + 
    ylab(paste0("Fraction of Images in C(X)"))
  print(plt.t)
  ggsave(paste0(save_dir, "year", t, "cdf.png"), width = 5, height = 3.2)
  
}


##### DIVERGENCE plots #####
# CDF of prior to get central regions
prior.ranks = rank(prior.depths) / length(prior.depths)
cr = central_region(prior, prior.ranks, 0.05)

# times = as.integer(seq(1, 998, length.out = 5))
times = c(1, 300, 600, 998)
for(t in times) {
  
  post = prep_post(nc.post, t)
  post = sapply(1:dim(post)[3], function(x) down_sample_image(post[,,x], downsamp))
  post.depths = meandepth(post, prior)
  
  pval.t = formatC(pvals[t], format = "e", digits = 1)
  print(field_plot(remove_cr(cr, post, downsamp),
                   nc.prior,
                   downsamp = downsamp,
                   main = paste0("Year ", t, " (p = ", pval.t, ")")))
  ggsave(paste0(save_dir, "year", t, "div.png"), width = 5, height = 3.2)
}


##### AVERAGE ALL #####
times = 1:998
downsamp = 1

# resetup prior
prior = prep_prior(nc.prior)
prior = sapply(1:dim(prior)[3], function(x) down_sample_image(prior[,,x], downsamp))

prior.depths = meandepth(prior, prior)
prior.ranks = rank(prior.depths) / length(prior.depths)
cr = central_region(prior, prior.ranks, 0.05)

# temp = matrix(0, (96/downsamp)*(144/downsamp), length(times))
# for(t in 1:length(times)) {
#   post = prep_post(nc.post, times[t])
#   post = sapply(1:dim(post)[3], function(x) down_sample_image(post[,,x], downsamp))
#   temp[,t] = as.vector(remove_cr(cr, post, downsamp))
#   
#   cat("Time ", times[t], "\n")
# }
# 
# temp = matrix(rowMeans(temp), (96/downsamp), (144/downsamp))
# field_plot(temp, nc.post, main = "Average temperature difference", downsamp = downsamp)


#### AVERAGE ANCIENT
times = (1:200)[as.logical(signif[1:200])]
temp = matrix(0, (96/downsamp)*(144/downsamp), length(times))
for(t in 1:length(times)) {
  post = prep_post(nc.post, times[t])
  post = sapply(1:dim(post)[3], function(x) down_sample_image(post[,,x], downsamp))
  temp[,t] = as.vector(remove_cr(cr, post, downsamp))
  
  cat("Time ", times[t], "\n")
}

temp = matrix(rowMeans(temp), (96/downsamp), (144/downsamp))
field_plot(temp, nc.post, main = "Average Exceedence (First 200 Years)", downsamp = downsamp)
ggsave(paste0(save_dir, "ancient.png"), width = 5, height = 3.2)


#### AVERAGE MODERN
times = (798:998)[as.logical(signif[798:998])]
temp = matrix(0, (96/downsamp)*(144/downsamp), length(times))
for(t in 1:length(times)) {
  post = prep_post(nc.post, times[t])
  post = sapply(1:dim(post)[3], function(x) down_sample_image(post[,,x], downsamp))
  temp[,t] = as.vector(remove_cr(cr, post, downsamp))
  
  cat("Time ", times[t], "\n")
}

temp = matrix(rowMeans(temp), (96/downsamp), (144/downsamp))
field_plot(temp, nc.post, main = "Average Exceedence (Last 200 Years)", downsamp = downsamp)
ggsave(paste0(save_dir, "modern.png"), width = 5, height = 3.2)


#### AVERAGE SIGNIFICANT
times = (1:998)[as.logical(signif)]
temp = matrix(0, (96/downsamp)*(144/downsamp), length(times))
for(t in 1:length(times)) {
  post = prep_post(nc.post, times[t])
  post = sapply(1:dim(post)[3], function(x) down_sample_image(post[,,x], downsamp))
  temp[,t] = as.vector(remove_cr(cr, post, downsamp))
  
  cat("Time ", times[t], "\n")
}

temp = matrix(rowMeans(temp), (96/downsamp), (144/downsamp))
field_plot(temp, nc.post, main = "Average Exceedence (All Years)", downsamp = downsamp)
ggsave(paste0(save_dir, "average.png"), width = 5, height = 3.2)


#### AVERAGE SIGNIFICANT TS
times = (1:998)
temp = matrix(0, (96/downsamp)*(144/downsamp), length(times))
for(t in 1:length(times)) {
  post = prep_post(nc.post, times[t])
  post = sapply(1:dim(post)[3], function(x) down_sample_image(post[,,x], downsamp))
  temp[,t] = as.vector(remove_cr(cr, post, downsamp))
  
  cat("Time ", times[t], "\n")
}

temp.ts = colMeans(temp)
plot(temp.ts)


##### ACTUAL CLIMATE ##### 
field_plot2 <- function(field, nc, main = "", downsamp = 1, zlim = c(-max(abs(field)), max(abs(field)))) {
  
  lats = as.vector(nc$dim$lat$vals)[seq(1, 96, by=downsamp)]
  lons = as.vector(nc$dim$lon$vals)[seq(1, 144, by=downsamp)]
  dimnames(field) = list(lats, ifelse(lons >= 180, lons - 360, lons))
  
  field.gg = melt(field)
  colnames(field.gg) = c("lat", "lon", "value")
  
  world = map_data("world")
  world = world[world$long <= 178, ]
  
  zlim = c(-max(abs(field)), max(abs(field)))
  
  ggplot() +
    geom_raster(data = field.gg, aes(x=lon, y=lat, fill=value), interpolate = TRUE) +
    geom_polygon(data = world, aes(x=long, y=lat, group=group), fill = NA, color="black") +
    coord_cartesian() +
    # scale_fill_gradient2(midpoint=0, low="#4393c3", mid="white", high="#d6604d") +
    # scale_fill_gradient2(midpoint=0, low="#2166ac", mid="white", high="#b2182b") +
    scale_fill_distiller(palette = "RdBu") +
    theme_void() +
    ggtitle(main) +
    theme(legend.position="none") + 
    theme(plot.title = element_text(hjust = 0.5))
}

post = prep_post(nc.post, times[1])
field_plot2(post[,,1], nc.post, main = "Posterior ensemble member (Year 1)", downsamp = downsamp)
ggsave(paste0(save_dir, "posterior.png"), width = 5, height = 3.2)

field_plot2(matrix(prior[,1], 96, 144), nc.post, main = "Prior ensemble member", downsamp = downsamp)
ggsave(paste0(save_dir, "prior.png"), width = 5, height = 3.2)
