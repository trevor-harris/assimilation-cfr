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

# FAD test
library(refund)

# code for running the QI and FAD tests
source("util/other_methods.R")

# code for importing and processing the ensemble data
source("util/import_data.R")

# prior
nc.prior = nc_open('data/tas_prior_da_hydro_r.1000-2000_d.16-Feb-2018.nc')
prior_ind = read.csv("data/prior_ens.txt", header = F)$V1

prior = prep_prior(nc.prior)
prior = prior[,,prior_ind]
prior = flatten(prior)

# post
nc.post = nc_open('data/tas_ens_da_hydro_r.1000-2000_d.16-Feb-2018.nc')

times = 998
years = 851:1848

# compute the FAD test p-value for each year in the reconstruction
fad = matrix(0, times, 2)
for(i in 1:times) {
  tic(paste0("Year ", i))
  
  post = prep_post(nc.post, i)
  post = flatten(post)
  
  fad[i, ] = fadtest(prior, post)[2]
  toc()
}

# plot
fadvals = data.frame(time = years, value = fad[,1])
ggplot(data = fadvals, aes(years, value)) +
  geom_point() +
  geom_smooth(color = "red", size = 2, se = F) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5), 
        text = element_text(size=20)) +
  ylab("FAD") +
  xlab("Year")

