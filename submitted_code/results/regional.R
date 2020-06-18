rm(list = ls())
gc()




########### READ ME #############

# you must change the working directory to be the submitted_code folder
# none of this will work otherwise
# mine is left here as an example

########## Example
# setwd("/Users/trevh/research/assimilation-cfr/submitted_code/")

#################################





library(ncdf4)
library(tictoc)
library(ggplot2)
library(dplyr)
library(reshape2)
library(latex2exp)
library(devtools)

devtools::install_github('trevor-harris/kstat')
library(kstat)

# code for importing and processing the ensemble data
source("util/import_data.R")


###### DATA PREP

# reg_ind = numbers used to identify regions in the mask file
reg_ind = c(c(10, 20, 30, 40, 50), c(1000, 2000, 3000, 4000, 5000, 6000, 7000))
reg_names = c(c("Arctic Ocean", "Indian Ocean", "Pacific Ocean", "Atlantic Ocean", "Southern Ocean"), 
              c("Antarctica", "South America", "North America", "Africa", "Europe", "Asia", "Australia"))

# open connection to NCDF4 data files
nc.post = nc_open('data/tas_ens_da_hydro_r.1000-2000_d.16-Feb-2018.nc')
nc.prior = nc_open('data/tas_prior_da_hydro_r.1000-2000_d.16-Feb-2018.nc')

# this file contains a list of ensemble members used to subset the background (prior)
prior_ind = read.csv("data/prior_ens.txt", header = F)$V1

# import all of prior then subset to ensemble members contained in prior_ind
prior = prep_prior(nc.prior)
prior = flatten(prior[,,prior_ind])

# Mask file identifies regions
mask = read.csv("data/mask.csv", stringsAsFactors = F)
mask = as.matrix(mask)
mask = apply(mask, 2, rev)


###### SIGNIFICANCE TESTING
# data dimensions
years = 851:1848 
times = 998
lats = 96
lons = 144
ens = 100
reg = length(reg_ind)
kfield = matrix(0, times, reg)
pfield = matrix(0, times, reg)

prior_reg = vector("list", length(reg_ind))
for(r in 1:length(reg_ind)) {
  prior_reg[[r]] =  matrix(prior[mask == reg_ind[r]], ncol = 100)
}

# run test for each year in reconstruction
for(i in 1:times) {
  tic(paste0("t = ", years[i]))
  
  post.t = flatten(prep_post(nc.post, i))
  post_reg = vector("list", length(reg_ind))
  for(r in 1:length(reg_ind)) {
    post_reg[[r]] =  matrix(post.t[mask == reg_ind[r]], ncol = 100)
  }
  
  for(r in 1:length(reg_ind)) {
    post_reg.r = post_reg[r][[1]]
    prior_reg.r = prior_reg[r][[1]]
    
    ktest = kstat(prior_reg.r, post_reg.r)
    kfield[i,r] = ktest[1]
    pfield[i,r] = ktest[2]
  }
  toc()
}

saveRDS(kfield , file = "results/kfield.RDS")
saveRDS(pfield , file = "results/pfield.RDS")

pval2 = p.adjust(pfield, "BY")
pval2 = matrix(pval2, nrow(pfield), ncol(pfield))

# unadjusted and FDR corrected pvals. FDR correction assumes all tests unique but this isnt true.
# this FDR correction was done originally and is left it as comparison with the correct one (pval3)
pval0 = k.t[,2]
pval1 = p.adjust(pval0, "BY")

# determine which ensembles are the same (out to 5 decimals) due to assimilating the exact same proxies
ex = matrix(0, 13824, times)
for(i in 1:times) {
  post.t = flatten(prep_post(nc.post, i))
  post.t = post.t - rowMeans(post.t)
  ex[,i] = round(post.t[, 10], 5)
}

# go through each par of saved ensembles and check if they are the same
equiv = matrix(F, times, times)
for(i in 1:(times-1)) {
  for(j in (i+1):times) {
    equiv[i, j] = identical(ex[,i], ex[,j])
  }
}

# orig = 1st appearance of a nonunique ensemble
# repeats = > 1st appearance
orig = which(apply(equiv, 1, sum) > 0)
repeats = which(apply(equiv, 2, sum) > 0)

pval0 = readRDS("results/pfield.RDS")

# Adjust only the years which are unique (drops 19 repeated years)
pval2 = p.adjust(pval0[-repeats,], "BY")
pval2 = matrix(pval2, length(pval2)/length(reg_ind), length(reg_ind))
pval3 = matrix(0, times, length(reg_ind))

# add in p-values for repeated years (use original instance p-value)
pval3[-repeats,] = pval2
pval3[repeats,] = pval2[orig]

# rename for consistency
pfield = pval3


# reassemble into series
kseries = melt(kfield) %>%
  mutate(Region = as.factor(Var2),
         Year = rep(years, 12),
         value = value)
levels(kseries$Region) = reg_names

pseries = melt(pfield) %>%
  mutate(Region = as.factor(Var2),
         Year = rep(years, 12),
         value = p.adjust(value, method = "BY"))
levels(pseries$Region) = reg_names


###### PLOT RESULTS

# plot KD values
ggplot(kseries, aes(Year, value)) +
  geom_point(size = 0.2) +
  geom_smooth(color = "red") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.position="none") +
  scale_x_continuous(breaks = c(1000, 1400, 1800)) +
  scale_y_continuous(breaks = c(0.1, 0.5, 1)) +
  scale_color_identity(guide = 'legend') +
  xlab("Year") +
  ylab("KD") +
  facet_wrap(vars(Region), 3, 4, scales = "free_x")
# ggsave("results/kd_region.png", width = 6, height = 3.8)


# plot p-values (grey out non siginificant values)
nonsigpseries = pseries[pseries$value > 0.05,]

ggplot() +
  geom_point(data = pseries, aes(Year, value), size = 0.2, color = "black") +
  geom_point(data = nonsigpseries, aes(Year, value), size = 0.2, color = "grey") +
  theme_classic() + 
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.position="none") +
  scale_x_continuous(breaks = c(1000, 1400, 1800)) +
  scale_y_continuous(breaks = c(0, 0.5, 1)) +
  scale_color_identity(guide = 'legend') +
  xlab("Year") +
  ylab("p-value") +
  facet_wrap(vars(Region), 3, 4, scales = "free_x")
# ggsave("results/pval_region.png", width = 6, height = 3.8)
