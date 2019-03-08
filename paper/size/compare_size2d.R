rm(list = ls()); gc()

library(reshape2)
library(ggplot2)

setwd("../research/assimilation-cfr/paper/size/")

ksize = readRDS("ksize2d.RDS")
qsize = readRDS("qsize2d.RDS")

ksize[["Statistic"]] = "K"
qsize[["Statistic"]] = "Q"

size = rbind(ksize, qsize)

size[["Statistic"]] = as.factor(size[["Statistic"]])
size[["Functions"]] = as.factor(size[["nfun"]])

ggplot(size, aes(x = Functions, y = size, fill = Statistic)) +
  geom_boxplot() +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5)) + 
  geom_hline(yintercept = 0.05) +
  xlab("Number of functions (N)") +
  facet_wrap(. ~ range)
ggsave("size2d.png", width = 7)

