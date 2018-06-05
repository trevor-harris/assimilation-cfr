rm(list = ls())
source("code/ks_field_functions.R")
source("code/sim_functions.R")

library(extdepth)
library(plotly)

# marginal number of points in the field (field is pts x pts)
pts = 40

# number of regions to subdivide the fields into
regions = 64

# standard flat prior mean
prior_mu = matrix(0, pts, pts)
post_mu = matrix(0, pts, pts)

prior_mu = as.vector(prior_mu)
post_mu = as.vector(post_mu)

prior = sim_gp(mu = prior_mu, l = 0.05, pts = pts)
post = sim_gp(mu = post_mu, l = 0.05, pts = pts)

# plot_ly(showscale = F) %>%
#   add_surface(z = ~prior[,,10])

# split em all
prior.split = vapply(1:100, function(x) matsplitter(prior[,,x], 5, 5),
                        FUN.VALUE = array(0, dim = c(5, 5, regions)))

post.split = vapply(1:100, function(x) matsplitter(post[,,x], 5, 5),
                       FUN.VALUE = array(0, dim = c(5, 5, regions)))

# ks = kst.permute2(prior.split, post.split, 100, 100)
# 
# plot(ks[["obs"]], type = "l", col = "red", ylim = c(0, 0.4))
# for(i in 1:100) {
#   lines(ks[["perm"]][,i])
# }

kst = kst.field(prior.split, post.split, 100)
ksp = kst.permute(prior.split, post.split, 500, 100)

microbenchmark(kst.permute(prior.split, post.split, 50, 100), times = 3)
microbenchmark(kst.permute(prior.split, post.split, 50, 2), times = 3)

# # Depth central regions
# perm.ed = edepth_set(ksp, depth_function = "rank")
# perm.cr = central_region(ksp, perm.ed)
# 
# kst.ed = edepth(kst, ksp, depth_function = "rank")
# 
# plot(kst, type = "l", col = "red", ylim = c(0, 0.4))
# for(i in 1:500) {
#   lines(ksp[,i])
# }
# lines(kst, col = "red")
# lines(perm.cr[[1]], col = "blue", lwd = 2)
# lines(perm.cr[[2]], col = "blue", lwd = 2)



# Depth central regions
perm.ed = edepth_set(ksp, depth_function = "rank")
perm.cr = central_region(ksp, perm.ed)
edepth(kst, ksp, depth_function = "rank")

# All functions
kst_df = data.frame(kst = kst, ind = 1:length(kst))
ksp_df = melt(ksp)
ksp_df[["ED"]] = rep(perm.ed, each=regions)

upp_df = data.frame(upper = perm.cr[[2]], ind = 1:length(kst))

ggplot() +
  geom_line(data = ksp_df, aes(x = Var1, y = value, group = ED, color = ED), alpha=0.5) +
  geom_line(data = kst_df, aes(x = ind, y = kst), color = "red") +
  geom_line(data = upp_df, aes(x = ind, y = upper), size = 1.1) +
  theme_classic() +
  theme(panel.background = element_rect(fill = "lightslategray")) +
  scale_color_distiller(palette = "YlGnBu") +
  xlab("Region") +
  ylab("Smoothed Permutation Distribution") +
  labs(color = "Ext Depth") +
  ggtitle("Permutation Distribution")



# Manual loess
colsmooth = function(fmat, degree) {
  fmat = as.matrix(fmat)
  fmat = data.frame(fmat, ind = 1:nrow(fmat))
  fmat = lapply(head(colnames(fmat), -1), 
                function(y) loess(as.formula(paste0(y, " ~ ind")), data=fmat, span=degree)$fitted)
  return(as.matrix(do.call(cbind, fmat)))
}

# Depth central regions
kst_df = data.frame(kst = kst, ind = 1:length(kst))
kst_df[["loess"]] = colsmooth(kst, 0.5)

ksp_sm = colsmooth(ksp, 0.5)
ksp_df = melt(ksp_sm)

perm.ed = edepth_set(ksp_sm, depth_function = "rank")
perm.cr = central_region(ksp_sm, perm.ed)
edepth(kst_df[["loess"]], ksp_sm, depth_function = "rank")

ksp_df[["ED"]] = rep(perm.ed, each=regions)


upp_df = data.frame(upper = perm.cr[[2]], ind = 1:length(kst))

ggplot() +
  geom_line(data = ksp_df, aes(x = Var1, y = value, group = ED, color = ED), alpha=0.5) +
  geom_line(data = kst_df, aes(x = ind, y = loess), color = "red") +
  geom_line(data = upp_df, aes(x = ind, y = upper), size = 1.1) +
  theme_classic() +
  theme(panel.background = element_rect(fill = "lightslategray")) +
  scale_color_distiller(palette = "YlGnBu") +
  xlab("Region") +
  ylab("Smoothed Permutation Distribution") +
  labs(color = "Ext Depth") +
  ggtitle("Permutation Distribution")


# ggplot() +
#   stat_smooth(data = ksp_df, aes(x = Var1, y = value, group = ED, color = ED), geom='line', alpha=0.5, se=FALSE) +
#   stat_smooth(data = kst_df, aes(x = ind, y = kst), geom='line', se=FALSE) +
#   theme_classic() +
#   theme(panel.background = element_rect(fill = "lightslategray")) +
#   scale_color_distiller(palette = "YlGnBu") +
#   xlab("Region") +
#   ylab("Smoothed Permutation Distribution") +
#   labs(color = "Ext Depth") +
#   ggtitle("Permutation Distribution")

# mean_ksp = rowMeans(ksp)
# plot(mean_ksp, type = "l", ylim = c(0, 0.4))
# mean_ksp[6:10]
# 
# var_ksp = apply(ksp, 1, sd)
# lines(var_ksp)
# lines(kst, col = "red")
# lines(perm.cr[[1]], col = "blue", lwd = 2)
# lines(perm.cr[[2]], col = "blue", lwd = 2)






