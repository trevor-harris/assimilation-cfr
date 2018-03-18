source('code/setup.R')

##### CONNECT TO DATA #####
# read ensembles and prior ncdf4 objects
nc.post = nc_open('/Users/trevh/Research/assimilation-cfr/data/tas_ens_da_hydro_r.1000-2000_d.16-Feb-2018.nc')
nc.prior = nc_open('/Users/trevh/Research/assimilation-cfr/data/tas_prior_da_hydro_r.1000-2000_d.16-Feb-2018.nc')


##### CONNECT TO DATA #####
# read ensembles and prior ncdf4 objects
nc.post = nc_open('/Users/Trevor/research/assimilation-cfr/data/tas_ens_da_hydro_r.1000-2000_d.16-Feb-2018.nc')
nc.prior = nc_open('/Users/Trevor/research/assimilation-cfr/data/tas_prior_da_hydro_r.1000-2000_d.16-Feb-2018.nc')


# commonly use 17x30 and 10x15
lat.basis = 28
lon.basis = 32

basis = create_basis(lat.basis, lon.basis, nc.prior)
proj = solve(t(basis) %*% basis) %*% t(basis)

# prep prior
prior.sub = read.csv("data/prior_ens.txt", header = F)
prior.sub = as.vector(prior.sub[,1])

prior.ens = prep_prior(nc.prior)
prior.ens = prior.ens[,,prior.sub]

# compute the central regions and find the ED of each prior ensemble member in the list
# also save the CR difference fields
ens_num = 1:100
eds = 1:100
for(e in ens_num) {
    
    # slice prior to member e
    prior = prior.ens[,,e]
    
    # prep posterior at member e
    post = prep_post_time(nc.post, e)
    
    # this time with only the last half of the ensemble
    # post = post[,,1:501]
    
    ##### FIT BASIS AND FIND ED #####
    prior.coef = proj %*% as.vector(prior)
    post.coef = sapply(1:dim(post)[3], function(x) proj %*% as.vector(post[,,x]))
    
    eds[e] = edepth(prior.coef, post.coef)
}

plot(eds, main = "Depth of each prior ensemble member", ylab = "Extremal Depth")

round(eds[c(1, 25, 50, 75, 100)], 3)
write.csv(eds, "edepths_28x32")

eds.df = data.frame(ed = eds, ind = 1:100)
ggplot(eds.df, aes(x = ind, y = ed)) +
  geom_point(aes(color = eds), show.legend = FALSE) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black")) +
  theme(plot.title = element_text(hjust = 0.5)) + 
  labs(x = "Ensemble No.", 
       y = "Extremal Depth", 
       title = "Depth of Prior Ensemble Members")

ggsave(paste0("paper/figures/prior_depths_tail", ".png"), width = 5, height = 3.2)




