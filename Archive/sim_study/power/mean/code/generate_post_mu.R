library(ggplot2)
library(reshape2)
source('sim_study/shared_code/gaussian_process.R')


set.seed(420)
post_mu = sim_gp(1, mu = 0, l = 10, 40)

mu_df = melt(post_mu*2)
ggplot(data = mu_df, aes(x = Var2, y = Var1)) +
  geom_raster(aes(fill = value), interpolate = T) +
  scale_fill_distiller(palette="Spectral") +
  labs(title = "Simulated Posterior Mean Field") +
  theme_void() +
  theme(plot.title = element_text(hjust = 0.5))

ggsave(paste0("sim_study/power/mean/plots/post_mu.png"), width = 5, height = 3.2)

saveRDS(matrix(post_mu, 40, 40), file = "sim_study/power/mean/data/post_mu.rds")

