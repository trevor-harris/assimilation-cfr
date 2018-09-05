source('code/setup.R')

##### VIZ THE INDIVIDUAL FIELDS #####
post_ens = prep_post_ens(nc.post, 100)
prior_ens = prior.ens[,,prior.sub]

post_df = melt(post_ens[,,1])
ggplot(data = post_df, aes(x = Var2, y = Var1)) +
  geom_raster(aes(fill = value), interpolate = T) +
  scale_fill_distiller(palette="Spectral") +
  labs(title = "Posterior Field") +
  theme_void() +
  theme(plot.title = element_text(hjust = 0.5))

# ggsave(paste0("plots/postfield.png"), width = 5, height = 3.2)

prior_df = melt(prior_ens[,,1])
ggplot(data = prior_df, aes(x = Var2, y = Var1)) +
  geom_raster(aes(fill = value), interpolate = T) +
  scale_fill_distiller(palette="Spectral") +
  labs(title = "Prior Field") +
  theme_void() +
  theme(plot.title = element_text(hjust = 0.5))

# ggsave(paste0("plots/priorfield.png"), width = 5, height = 3.2)



# Add in the grid lines



# t = 100 is ALL different
t = 900
prior_ens = prior.ens[,,prior.sub]
post_ens = prep_post_ens(nc.post, t)

regions = 216
lil_priors = vapply(1:100, function(x) matsplitter(prior_ens[,,x], 8, 8),
                    FUN.VALUE = array(0, dim = c(8, 8, regions)))


lil_posts  = vapply(1:100, function(x) matsplitter(post_ens[,,x], 8, 8),
                    FUN.VALUE = array(0, dim = c(8, 8, regions)))


# find the observed kst field
kol.field = kst.field(lil_priors, lil_posts, 100)

# find the permutation distribution
perm.fields = kst.permute(lil_priors, lil_posts, 100, 100)

perm.ed = edepth_set(perm.fields, depth_function = "rank")
perm.cr = central_region(perm.fields, perm.ed)

ymin = min(c(kol.field, perm.cr[[2]]))
ymax = max(c(kol.field, perm.cr[[2]]))

plot(kol.field, type = "l", col = "red", ylim = c(ymin, ymax), main = paste0("t = ", t))
lines(perm.cr[[2]], lwd = 2)
