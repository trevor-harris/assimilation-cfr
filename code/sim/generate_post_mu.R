library(plotly)

pts = 40

remove_IQR = function(x, range = 0.5) {
  lower = 0.5 - range/2
  upper = 0.5 + range/2
  thresh = max(abs(quantile(gp_mu, lower)), abs(quantile(gp_mu, upper)))
  ifelse(abs(x) < thresh, 0, ifelse(x < 0, x+thresh, x-thresh))
}

set.seed(102393)
gp_mu = matrix(sim_gp(fields = 1, mu = 0, l = 1, pts = 40, ker = "exp"), 40, 40)
gp_mu = gp_mu - mean(gp_mu)
gp_mu_thre = remove_IQR(gp_mu, 0.1) / 2

summary(as.vector(gp_mu_thre))

plot_ly(showscale = F) %>%
  add_surface(z = ~gp_mu_thre)


saveRDS(matrix(gp_mu_thre, pts, pts), file = paste0("/Users/trevh/research/assimilation-cfr/simdata/run3/post_mu.rds"))
# image(gp_mu_thre)



