library(ggplot2)
library(reshape2)
source('code/sim_functions.R')

set.seed(420)
post_mu = sim_gp(1, mu = 0, l = 1, 40)

mu_df = melt(post_mu)
ggplot(data = mu_df, aes(x = Var2, y = Var1)) +
  geom_raster(aes(fill = value), interpolate = T) +
  scale_fill_distiller(palette="Spectral") +
  labs(title = "Simulated Posterior Mean Field") +
  theme_void() +
  theme(plot.title = element_text(hjust = 0.5))

ggsave(paste0("sim_study/power/mean/plots/post_mu.png"), width = 5, height = 3.2)

saveRDS(matrix(post_mu, 40, 40), file = paste0("/Users/trevh/research/assimilation-cfr/sim_study/power/mean/data/post_mu.rds"))



