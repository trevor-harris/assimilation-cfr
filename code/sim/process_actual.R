rm(list = ls())
gc()

library(extdepth)
source("/Users/trevh/research/assimilation-cfr/code/ks_field_functions.R")
source("/Users/trevh/research/assimilation-cfr/code/prep_functions.R")
source("/Users/trevh/research/assimilation-cfr/code/sim_functions.R")

data.dir = "/Users/trevh/research/assimilation-cfr/cfrdata/run3/"
files = list.files(data.dir)

regions = 216
batches = 100
times = 998

stats = matrix(0, times, regions)
upper = matrix(0, times, regions)

for(batch in 1:batches) {
  stats.dat = readRDS(paste0(data.dir, "values", batch, ".rds"))
  upper.dat = readRDS(paste0(data.dir, "upper", batch, ".rds"))
  
  
  ind = 1:nrow(stats.dat) + 10*(batch - 1)
  
  for(i in ind) {
    stats[i, ] = stats.dat[i - 10*(batch - 1), ]
    upper[i, ] = upper.dat[i - 10*(batch - 1), ]
  }
}


# Individual functions against there upper bounds
for (t in c(1, 250, 750, 998)) {
  ksku_df = melt(t(rbind(stats[t,], upper[t,])))
  ksku_df[["Var2"]] = as.factor(ksku_df[["Var2"]])
  levels(ksku_df[["Var2"]]) <- list("Kolmogorov"="1", "Upper Bound"="2")
  
  plt = ggplot(data = ksku_df, aes(x = Var1, y = value, group = Var2)) +
    # stat_smooth(geom='line', se = TRUE, level = 0.95, aes(color = Var2)) +
    geom_smooth(aes(color = Var2), alpha = 0.3) +
    geom_line(aes(color = Var2)) +
    theme_classic() +
    xlab("Region") +
    ylab("Smoothed Kolmogorov-Smirnov") +
    labs(color = "") +
    ggtitle(paste0("Kolmogorov-Smirnov function vs Upper Bound at t = ", t))
  
  print(plt)
}


# All functions
KS_df = melt(t(stats))
ggplot(data = KS_df, aes(x = Var1, y = value, group = Var2)) +
  stat_smooth(geom='line', alpha=0.5, se=FALSE, aes(color = Var2)) +
  theme_classic() +
  theme(panel.background = element_rect(fill = "lightslategray")) +
  scale_color_distiller(palette = "YlGnBu") +
  xlab("Region") +
  ylab("Smoothed Kolmogorov-Smirnov") +
  labs(color = "Time") +
  ggtitle("Mean Kolmogorov-Smirnov function over Time")




