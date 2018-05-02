library(plotly)

pts = 40

remove_IQR = function(x) {
  thresh = min(abs(quantile(x, 0.25)), abs(quantile(x, 0.75)))
  ifelse(abs(x) < thresh, 0, ifelse(x < 0, x+thresh, x-thresh))
}

set.seed(4)
gp_mu1 = matrix(sim_gp(fields = 1, mu = 0, l = 2, pts = 40), 40, 40) 
gp_mu2 = matrix(sim_gp(fields = 1, mu = 0, l = 0.1, pts = 40), 40, 40)
gp_mu3 = matrix(sim_gp(fields = 1, mu = 0, l = 0.05, pts = 40), 40, 40) 
gp_mu = (gp_mu1 + gp_mu2 + gp_mu3)
gp_mu = gp_mu - mean(gp_mu)
gp_mu_thre = sign(gp_mu) * sqrt(abs(remove_IQR(gp_mu))) / 2

plot_ly(showscale = F) %>%
  add_surface(z = ~gp_mu_thre)

saveRDS(matrix(gp_mu_thre, pts, pts), file = paste0("/Users/trevh/research/assimilation-cfr/simdata/post_mu.rds"))
image(gp_mu_thre)


