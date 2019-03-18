rm(list = ls()); gc()

library(reshape2)
library(ggplot2)

setwd("../research/assimilation-cfr/paper/size/out2")
setwd("../out2/")

ksize = rbind(readRDS("ksize2d_1.RDS"), readRDS("ksize2d_2.RDS"))
qsize = rbind(readRDS("qsize2d_1.RDS"), readRDS("qsize2d_2.RDS"))

ksize[["Statistic"]] = "K"
qsize[["Statistic"]] = "Q"

size = rbind(ksize, qsize)

size[["Statistic"]] = as.factor(size[["Statistic"]])
size[["Functions"]] = as.factor(size[["nfun"]])

size = size[size$range < 25, ]

ggplot(size, aes(x = Functions, y = size, fill = Statistic)) +
  geom_boxplot() +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5)) + 
  geom_hline(yintercept = 0.05) +
  xlab("Number of functions (n)") +
  ylab("Rejection proportion") +
  facet_wrap(. ~ range, nrow = 1)
ggsave("../size2d2.png", width = 9, height = 3)
