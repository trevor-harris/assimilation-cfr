##### HORIZONTAL #####

# define data and import fuctions
source('code/setup.R')



#### OVERALL MAPS #####
# this is run on the main computer. Do not run on laptop



#### ED vs Coeff number #####
ens.rng = 1:100
eds = 1:100
for(e in ens.rng) {
  
  # slice prior to member e
  prior = prior.ens[,,e]
  
  # prep posterior at member e
  post = prep_post_time(nc.post, e)
  
  # this time with only the last half of the ensemble
  # post = post[,,1:501]
  
  # FIT BASIS AND FIND ED
  prior.coef = proj %*% as.vector(prior)
  post.coef = sapply(1:dim(post)[3], function(x) proj %*% as.vector(post[,,x]))
  
  eds[e] = edepth(prior.coef, post.coef)
}

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

ggsave(paste0("paper/figures/prior_depths", ".png"), width = 5, height = 3.2)


#### Difference Fields #####
ens.rng = c(10, 25, 45, 76)
cr_diff_fields = array(0, dim = c(dim(prior.ens)[1], dim(prior.ens)[2], length(ens.rng)))
for(e in ens.rng) {
  
  # slice prior to member e
  prior = prior.ens[,,e]
  
  # prep posterior at member e
  post = prep_post_time(nc.post, e)
  
  # FIT BASIS AND FIND ED
  prior.coef = proj %*% as.vector(prior)
  post.coef = sapply(1:dim(post)[3], function(x) proj %*% as.vector(post[,,x]))
  
  # calculate extremal depth
  ed = edepth_set(post.coef)
  central = central_region(post.coef, ed, 0.5)
  
  # Use the x% central regions to measure outlyingness
  alpha.dev = rep(0, length(prior.coef))
  for (a in 1:length(alpha.dev)) {
    if (central[[1]][a] > prior.coef[a]) {
      alpha.dev[a] = prior.coef[a] - central[[1]][a]
    } 
    if (prior.coef[a] > central[[2]][a]) {
      alpha.dev[a] = prior.coef[a] - central[[2]][a]
    }
  }
  prior.vec = basis %*% alpha.dev
  prior.diff = matrix(prior.vec , nrow(prior), ncol(prior))
  
  cr_diff_fields[,,which(ens.rng == e)] = prior.diff
  # eds_sm[which(ens.rng == e)] = edepth(prior.coef, post.coef)
}

# plot the difference fields
for(e in 1:length(ens.rng)) {
  print(field_plot(cr_diff_fields[,,e], nc.prior, main = paste0("Ensemble ", ens.rng[e]), zlim = c(-2, 2)))
  ggsave(paste0("paper/figures/ens_diff_", ens.rng[e], ".png"), width = 5, height = 3.2)
}

# corresponding ED values
round(eds[ens.rng], 3)


#### Coef Examples ####
ens.rng = c(10, 76)
# cr_diff_fields = array(0, dim = c(dim(prior.ens)[1], dim(prior.ens)[2], length(ens.rng)))
for(e in ens.rng) {
  
  # slice prior to member e
  prior = prior.ens[,,e]
  
  # prep posterior at member e
  post = prep_post_time(nc.post, e)
  
  # FIT BASIS AND FIND ED
  prior.coef = proj %*% as.vector(prior)
  post.coef = sapply(1:dim(post)[3], function(x) proj %*% as.vector(post[,,x]))
  
  ed = edepth_set(post.coef)
  central = central_region(post.coef, ed, 0.05)
  median = post.coef[,ed == 1]
  
  # create the prior/posterior comparison graph data
  coef.df = data.frame(ind = 1:1040, post.coef)
  coef.df = coef.df[coef.df$ind <= 15,]
  coef.df = melt(coef.df, id.vars = "ind")
  
  prior.df = data.frame(ind = 1:1040, prior=prior.coef)
  prior.df = prior.df[prior.df$ind <= 15,]
  
  # create the prior/CR regions comparison graph data
  cr.df = data.frame(ind = 1:1040, lower=central[[1]], upper=central[[2]])
  cr.df = cr.df[cr.df$ind <= 15,]
  cr.df = melt(cr.df, id.vars = "ind")
  
  med.df = data.frame(ind = 1:1040, median=median)
  med.df = med.df[med.df$ind <= 15,]
  
  all.plot = ggplot() +
    geom_line(data = coef.df, aes(x = ind, y = value, group = variable)) +
    geom_line(data = prior.df, aes(x = ind, y = prior, col = "Prior"),
              show.legend = FALSE) +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank(), 
          axis.line = element_line(colour = "black")) +
    theme(plot.title = element_text(hjust = 0.5)) +
    labs(x = "Coefficient",
         y = "Value",
         title = paste0(" Ens.", e, ": Prior vs Posterior"))
  print(all.plot)
  ggsave(paste0("paper/figures/prior_post_coef_", e, ".png"), width = 5, height = 3.2)
  
  cr.plot = ggplot() +
    geom_line(data = cr.df, aes(x = ind, y = value, group = variable)) +
    geom_line(data = med.df, aes(x = ind, y = median, color = "#F8766D"),
              show.legend = FALSE) +
    geom_line(data = prior.df, aes(x = ind, y = prior, color = "#00BFC4"),
              show.legend = FALSE) +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank(), 
          axis.line = element_line(colour = "black")) +
    theme(plot.title = element_text(hjust = 0.5)) +
    labs(x = "Coefficient",
         y = "Value",
         title = paste0(" Ens.", e, ": Prior vs Central Regions"),
         fill = "sdfa")
  print(cr.plot)
  ggsave(paste0("paper/figures/prior_cr_coef_", e, ".png"), width = 5, height = 3.2)
}
