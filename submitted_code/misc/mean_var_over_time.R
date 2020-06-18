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
library(future)
library(future.apply)

# FAD test
library(refund)

# code for importing and processing the ensemble data
source("util/import_data.R")

# prior
nc.prior = nc_open('data/tas_prior_da_hydro_r.1000-2000_d.16-Feb-2018.nc')
prior_ind = read.csv("data/prior_ens.txt", header = F)$V1

prior = prep_prior(nc.prior)
prior = prior[,,prior_ind]

# post
nc.post = nc_open('data/tas_ens_da_hydro_r.1000-2000_d.16-Feb-2018.nc')

# data dimensions
times = 998
years = 851:1848
ens = 100

# compute the pointwise, within ensemble, mean and standard deviation
# of the background (prior) and the analysis (post)
# Compute the mean squared difference between means
# compute the mean squared ratio of the standard deviations
param = matrix(0, 998, 2)
for(y in 1:10) {
  tic(paste0("Year ", y))
  post = prep_post(nc.post, y)

  mu1 = apply(prior, c(1, 2), mean)
  mu2 = apply(post, c(1, 2), mean)
  
  sd1 = apply(prior, c(1, 2), sd)
  sd2 = apply(post, c(1, 2), sd)

  param[y,1] = mean((mu1 - mu2)^2)
  param[y,2] = mean((sd1 / sd2)^2)

  toc()
}

plot(param[,1])
plot(param[,2])

# plot mean MSEs
mu = data.frame(time = years, value = param[,1])
ggplot(data = mu, aes(years, value)) +
  geom_point() +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5), 
        text = element_text(size=20)) +
  ylab("Average difference squared") +
  xlab("Year") +
  ggtitle("Average difference between background and analysis means")
# ggsave("../research/assimilation-cfr/paper/misc/means.png", width = 8, heigh = 6)

# plot sd mean square ratios
sig = data.frame(time = years, value = param[,2])
ggplot(data = sig, aes(years, value)) +
  geom_point() +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5), 
        text = element_text(size=20)) +
  ylab("Average ratio squared") +
  xlab("Year") +
  ggtitle("Average ratio of background over analysis standard deviations")
# ggsave("../research/assimilation-cfr/paper/misc/sds.png", width = 8, heigh = 6)



