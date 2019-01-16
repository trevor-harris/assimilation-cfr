rm(list = ls()); gc()

library(reshape2)
library(ggplot2)

size = readRDS("research/assimilation-cfr/paper/size/independent1d/sim1d.RDS")
size[["Statistic"]] = as.factor(size[["Stat"]])
levels(size[["Statistic"]]) = c("K", "Q")

ggplot(size, aes(x = Functions, y = Size, color = Statistic)) + 
  geom_boxplot() +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5)) + 
  geom_hline(yintercept = 0.05) +
  xlab("Number of functions (N)")
  # ggtitle("Size of K(F, G) v.s. Q(F, G)")

ggsave("research/assimilation-cfr/paper/size/size1d.png", width = 5, height = 3.2)
