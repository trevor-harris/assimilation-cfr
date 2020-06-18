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

# open connection to NCDF4 data files
nc.post = nc_open('data/tas_ens_da_hydro_r.1000-2000_d.16-Feb-2018.nc')
nc.prior = nc_open('data/tas_prior_da_hydro_r.1000-2000_d.16-Feb-2018.nc')

# this file contains a list of ensemble members used to subset the background (prior)
prior_ind = read.csv("data/prior_ens.txt", header = F)$V1

# import all of prior then subset to ensemble members contained in prior_ind
prior = prep_prior(nc.prior)
prior = flatten(prior[,,prior_ind])


###### SIGNIFICANCE TESTING
times = 998
years = 851:1848
k.t = matrix(0, times, 2)

# run test for each year in reconstruction
for(i in 1:times) {
  tic(paste0("t = ", years[i]))
  post.t = flatten(prep_post(nc.post, i))
  k.t[i,] = kstat(prior, post.t)
  toc()
}

# unadjusted and FDR corrected pvals. FDR correction assumes all tests unique but this isnt true.
# this FDR correction was done originally and is left it as comparison with the correct one (pval3)
pval0 = k.t[,2]
pval1 = p.adjust(pval0, "BY")

# determine which ensembles are the same (out to 5 decimals) due to assimilating the exact same proxies
times = 998
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

pval2 = p.adjust(pval0[-repeats], "BY")
pval3 = rep(0, times)
pval3[-repeats] = pval2
pval3[repeats] = pval2[orig]

pval1 = readRDS("results/global_pval_full_adjusted.RDS")
pval3 = readRDS("results/global_pval_correct_adjusted.RDS")
plot(pval1 - pval3)

max(pval3 - pval1)
min(pval3 - pval1)

saveRDS(pval0, "results/global_pval_unadjusted.RDS");
saveRDS(pval1, "results/global_pval_full_adjusted.RDS");
saveRDS(pval3, "results/global_pval_correct_adjusted.RDS");


###### PLOT RESULTS
kstats = data.frame(stat = k.t[,1], pval = k.t[,2], Year = years)

ggplot(kstats, aes(x=Year, y=stat)) +
  geom_point(size = 2) +
  geom_smooth(color = "red", size = 2, se = F) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5), 
        text = element_text(size=24)) +
  ylab("KD")
# ggsave("results/kd_global.png", width = 8, heigh = 6)


ggplot(kstats, aes(x=Year, y=pval)) +
  geom_point(size = 2) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5), 
        text = element_text(size=24)) +
  ylab("p-values")
# ggsave("results/pval_global.png", width = 8, heigh = 6)