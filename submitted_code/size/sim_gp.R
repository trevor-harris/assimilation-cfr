rm(list = ls())
gc()



########### READ ME #############

# you must change the working directory to be the submitted_code folder
# none of this will work otherwise
# mine is left here as an example

########## Example
# setwd("/Users/trevh/research/assimilation-cfr/submitted_code/")

#################################




library(tictoc)
library(future.apply)

# FAD test
library(refund)

# Band test
devtools::install_version('roahd', version = "1.4.1", repos = "http://cran.us.r-project.org")
library(roahd)

# KD test
devtools::install_github('trevor-harris/kstat')
library(kstat)

# code for simulating guassian processes, t-processes, and plotting functions
source("util/simulation.R")

# code for running the QI and FAD tests
source("util/other_methods.R")

#### SIZE
plan(multicore)

i=1
sims = 2000

# do not actually run this entire loop. The inner part was run in parallel on a cluster.
for (n1 in c(50,100,200,300)) {
  for (n2 in c(50,100,200,300)) {
    for (rng in c(0.2,0.3,0.4,0.5)) {
      for (nu in c(0.5,1.0,1.5)) {
        
        set.seed(072393 + i)
        
        vals = future_sapply(1:sims, function(x) {
          
          # simulate 2d gaussian process data
          f = gp2d(fields = n1, range = rng, nu = nu)
          g = gp2d(fields = n2, range = rng, nu = nu)
          
          # flattent into 1d vectors
          f = flatten(f)
          g = flatten(g)
          
          # compute each competing statistic on the generated data and save the p-value
          c(kstat(f, g)[2], quality(f, g)[2], fadtest(f, g)[2], bandtest(f, g)[2])
          
        })
        
        size = data.frame(pval = c(vals[1,], vals[2,], vals[3,], vals[4,]),
                          method = rep(c("Kolm", "Qual", "FAD", "MBD"), each=sims),
                          rng = rng,
                          nu = nu,
                          n1 = n1,
                          n2 = n2)
        
        saveRDS(size, file = paste0("size/out_gp/","size_sim_",i,".RDS"))
        
        i = i + 1
      }
    }
  }
}

