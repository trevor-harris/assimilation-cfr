library(plotly)

pts = 40

remove_IQR = function(x) {
  thresh = min(abs(quantile(x, 0.25)), abs(quantile(x, 0.75)))
  ifelse(abs(x) < thresh, 0, ifelse(x < 0, x+thresh, x-thresh))
}

set.seed(102391)
gp_mu = matrix(sim_gp(fields = 1, mu = 0, l = 0.2, pts = 40), 40, 40) 
gp_mu = gp_mu - mean(gp_mu)
gp_mu_thre = remove_IQR(gp_mu)

# plot_ly(showscale = F) %>%
#   add_surface(z = ~gp_mu_thre)

saveRDS(matrix(post_mu, pts, pts), file = paste0("/Users/trevh/research/assimilation-cfr/simdata/post_mu.rds"))
# image(gp_mu_thre)


