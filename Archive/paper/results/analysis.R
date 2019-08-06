rm(list = ls())
gc()

years = 851:1848
ens = 100

library(extdepth)
library(ncdf4)
library(dplyr)
library(reshape2)
library(ggplot2)
library(OpenImageR)
library(tictoc)
library(future)
library(future.apply)

source("research/assimilation-cfr/code/depth_tests.R")
source("research/assimilation-cfr/code/depths.R")
source("research/assimilation-cfr/code/simulation.R")

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
read_era = function(dir, file) {
  cbind(readRDS(paste0(dir, file)), era = as.numeric(strsplit(file, "\\D+")[[1]][-1]))
}
field_plot <- function(field, nc, main = "", downsamp = 1, zlim = c(-max(abs(field)), max(abs(field)))) {
  
  lats = as.vector(nc$dim$lat$vals)[seq(1, 96, by=downsamp)]
  lons = as.vector(nc$dim$lon$vals)[seq(1, 144, by=downsamp)]
  dimnames(field) = list(lats, ifelse(lons >= 180, lons - 360, lons))
  
  field.gg = melt(field)
  colnames(field.gg) = c("lat", "lon", "Temp")
  
  world = map_data("world")
  world = world[world$long <= 178, ]
  
  zlim = c(-max(abs(field)), max(abs(field)))
  
  ggplot() +
    geom_raster(data = field.gg, aes(x=lon, y=lat, fill=Temp), interpolate = TRUE) +
    geom_polygon(data = world, aes(x=long, y=lat, group=group), fill = NA, color="black") +
    coord_cartesian() +
    scale_fill_gradient2(midpoint=0, low="blue", mid="white", high="red") +
    # scale_fill_gradient2(midpoint=0, low="#4393c3", mid="white", high="#d6604d") +
    # scale_fill_gradient2(midpoint=0, low="#2166ac", mid="white", high="#b2182b") +
    theme_void() +
    ggtitle(main) +
    theme(plot.title = element_text(hjust = 0.5))
}

# regionalizations
matsplitter = function(M, r, c) {
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
matcombiner = function(arr, r) {
  c = dim(arr)[3] / r
  mlist = lapply(seq(dim(arr)[3]), function(x) arr[,,x])
  do.call(rbind, lapply(0:(c-1), function(x) do.call(cbind, mlist[(1:r) + x*c])))
}
rmap = function(val, r, c) {
  vapply(val, function(x) matrix(x, r, c), matrix(0, r, c))
}

# equal area regionalization
rescale = function(vec, l, u) {
  r = ((vec - min(vec)) / (max(vec) - min(vec))) * (u - l) + l
  as.integer(r)
}
area_splitter = function(arr, regions = 16, rx = sqrt(regions)) {
  mrow = dim(arr)[1]
  mcol = dim(arr)[2]
  mens = dim(arr)[3]
  
  rx = rx
  ry = regions / rx
  
  indx = 1:mrow - (1+mrow)/2
  indx = -acos(indx / (max(indx)))
  indx = indx[seq(1, mrow, length.out = rx+1)]
  indx = rescale(indx, 1, mrow)
  
  indy = as.integer(seq(1, mcol, length.out = ry+1))
  
  out = lapply(1:rx, function(x) {
          lapply(1:ry, function(y) arr[(indx[x] + as.integer(x>1)):indx[x+1], 
                                       (indy[y] + as.integer(y>1)):indy[y+1],])
        })
  unlist(out, recursive = FALSE)
}


save_dir = "research/assimilation-cfr/paper/results/"

# import raw size data
# dir = "../temp/power/independent/"
dir = "../research/assimilation-cfr/paper/results/results/"
files = list.files(dir)

temperature = read_era(dir, files[1])
for(f in 3:length(files)) {
  temperature = rbind(temperature, read_era(dir, files[f]))
}
temperature = rbind(temperature, read_era(dir, files[2]))
temperature[["time"]] = years

temperature$stat = temperature$stat / (sqrt((ens*ens)/(ens + ens)))

ggplot(temperature, aes(time, stat)) +
  geom_point(alpha = 0.7) +
  geom_smooth(se = F, color = "red") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5)) +
  ylab("K") +
  xlab("Year")
  # ggtitle("K(F, G) over time")
ggsave("../research/assimilation-cfr/paper/results/effect_over_time.png", width = 5, height = 3.2)

temperature$pval = p.adjust(temperature$pval)

ggplot(temperature, aes(time, pval)) +
  geom_point(alpha = 0.7) +
  # geom_smooth(se = F, color = "red") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5)) +
  ylab("p-values")+
  xlab("Year")
  # ggtitle("P-values of K(F, G) over time")
ggsave("../research/assimilation-cfr/paper/results/pval_over_time.png", width = 5, height = 3.2)




#### Regaionlized significant differences
nc.post = nc_open('research/assimilation-cfr/data/tas_ens_da_hydro_r.1000-2000_d.16-Feb-2018.nc')
nc.prior = nc_open('research/assimilation-cfr/data/tas_prior_da_hydro_r.1000-2000_d.16-Feb-2018.nc')
prior_ind = read.csv("research/assimilation-cfr/data/prior_ens.txt", header = F)$V1


# 16 region setting
times = 998
lats = 96
lons = 144
ens = 100
regions = 16
rs = sqrt(regions)
kfield = matrix(0, regions, times)
pfield = matrix(0, regions, times)

plan(multiprocess)

# import all of prior
prior = prep_prior(nc.prior)
prior = prior[,,prior_ind]
prior.split = area_splitter(prior, regions)

post = vapply(1:times, function(t) prep_post(nc.post, t), array(0, dim=c(lats, lons, ens)))


for(t in 1:times) {
  tic(paste0("Time: ", t))
  
  post.t = post[,,,t]
  post.split = area_splitter(post.t)
  
  ktest = future_sapply(1:regions, function(x) kolm(flatten(prior.split[[x]]), flatten(post.split[[x]])))
  kfield[,t] = ktest[1,]
  pfield[,t] = ktest[2,]
  
  toc()
}

# reassemble into series
kseries = melt(t(kfield)) %>% 
  mutate(Region = as.factor(Var2),
         value = value / sqrt(ens*ens / (ens + ens)))

pseries = melt(t(pfield)) %>%
  mutate(Region = as.factor(Var2),
         value = p.adjust(value, method = "BY"))

ggplot(kseries, aes(Var1, value, color = Region)) +
  geom_point(size = 0.2) +
  geom_smooth(color = "black") +
  theme_classic() + 
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.position="none") +
  scale_x_continuous(breaks = c(1, 500, 998)) +
  scale_y_continuous(breaks = c(0, 0.5, 1)) +
  xlab("Time") +
  ylab("K") +
  ggtitle("K over time by region") +
  facet_wrap(vars(Region), 4, 4)
ggsave(paste0(save_dir, "k_region.png"), width = 5, height = 3.2)

ggplot(pseries, aes(Var1, value, color = Region)) +
  geom_point(size = 0.2) +
  theme_classic() + 
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.position="none") +
  scale_x_continuous(breaks = c(1, 500, 998)) +
  scale_y_continuous(breaks = c(0, 0.5, 1)) +
  xlab("Time") +
  ylab("p-value") +
  ggtitle("P-values over time by region") +
  facet_wrap(vars(Region), 4, 4)
ggsave(paste0(save_dir, "pval_region.png"), width = 5, height = 3.2)


# # 64 region setting
# times = 998
# lats = 96
# lons = 144
# ens = 100
# regions = 64
# rs = sqrt(regions)
# kfield = matrix(0, regions, times)
# pfield = matrix(0, regions, times)
# 
# # import all of prior
# prior.split = area_splitter(prior, regions)
# 
# 
# for(t in 1:times) {
#   tic(paste0("Time: ", t))
#   
#   post.t = post[,,,t]
#   post.split = area_splitter(post.t, regions)
#   
#   ktest = future_sapply(1:regions, function(x) kolm(flatten(prior.split[[x]]), flatten(post.split[[x]])))
#   kfield[,t] = ktest[1,]
#   pfield[,t] = ktest[2,]
#   
#   toc()
# }
# 
# # reassemble into series
# kseries = melt(t(kfield)) %>% 
#   mutate(Region = as.factor(Var2))
# 
# pseries = melt(t(pfield)) %>%
#   mutate(Region = as.factor(Var2),
#          value = p.adjust(value, method = "BY"))
# 
# ggplot(kseries, aes(Var1, value, color = Region)) +
#   geom_point(size = 0.2) +
#   geom_smooth(color = "black") +
#   theme_classic() + 
#   theme(plot.title = element_text(hjust = 0.5)) +
#   theme(legend.position="none") +
#   scale_x_continuous(breaks = c(1, 500, 998)) +
#   scale_y_continuous(breaks = c(0, 0.5, 1)) +
#   xlab("Time") +
#   ylab("K") +
#   ggtitle("K over time by region") +
#   facet_wrap(vars(Region), 8, 8)
# # ggsave(paste0(save_dir, "k_region.png"), width = 5, height = 3.2)
# 
# ggplot(pseries, aes(Var1, value, color = Region)) +
#   geom_point(size = 0.2) +
#   theme_classic() + 
#   theme(plot.title = element_text(hjust = 0.5)) +
#   theme(legend.position="none") +
#   scale_x_continuous(breaks = c(1, 500, 998)) +
#   scale_y_continuous(breaks = c(0, 0.5, 1)) +
#   xlab("Time") +
#   ylab("p-value") +
#   ggtitle("P-values over time by region") +
#   facet_wrap(vars(Region), 8, 8)
# # ggsave(paste0(save_dir, "pval_region.png"), width = 5, height = 3.2)
# 


# reassemble into maps
kmaps = vapply(1:times, function(t) {
  matcombiner(rmap(kfield[,t], lats/4, lons/4), 4)
}, matrix(0, lats, lons))

pmaps = vapply(1:times, function(t) {
  matcombiner(rmap(pfield[,t], lats/4, lons/4), 4)
}, matrix(0, lats, lons))

field_plot(kmaps[,,3], nc.post)
field_plot(pmaps[,,3], nc.post)

