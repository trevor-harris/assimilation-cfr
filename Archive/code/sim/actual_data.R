regions = 216
# times = as.integer(seq(1, 999, length.out = 20))
times = 1

KS = matrix(0, length(times), regions)
KU = matrix(0, length(times), regions)

prior_ens = prior.ens[,,prior.sub]
for (t in times) {
  prior_ens = prior.ens[,,prior.sub]
  post_ens = prep_post_ens(nc.post, t)
  
  lil_priors = vapply(1:100, function(x) matsplitter(prior_ens[,,x], 8, 8),
                      FUN.VALUE = array(0, dim = c(8, 8, regions)))
  
  
  lil_posts  = vapply(1:100, function(x) matsplitter(post_ens[,,x], 8, 8),
                      FUN.VALUE = array(0, dim = c(8, 8, regions)))
  
  
  # find the observed kst field
  kol.field = kst.field(lil_priors, lil_posts, 100)
  
  # find the permutation distribution
  perm.fields = kst.permute(lil_priors, lil_posts, 100, 100)
  
  # perm.ed = edepth_set(perm.fields, depth_function = "rank")
  # perm.cr = central_region(perm.fields, perm.ed)
  # 
  # KS[which(times==t), ] = kol.field
  # KU[which(times==t), ] = perm.cr[[2]]
  
  perm.upper = sapply(1:length(kol.field), function(r) quantile(perm.fields[r,], 0.95))
  
  KS[which(times==t), ] = kol.field
  KU[which(times==t), ] = perm.upper
}

# save.image("actual.RData")

# 
KS_df = melt(t(KS))
ggplot(data = KS_df, aes(x = Var1, y = value, group = Var2)) +
  geom_smooth(aes(color = Var2), level = 0) +
  theme_classic()




