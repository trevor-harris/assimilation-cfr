rm(list = ls()); gc()

library(reshape2)
library(ggplot2)
library(reshape2)
library(dplyr)

setwd("../../size/out/")

ksize = readRDS("big_ksize.RDS")
qsize = readRDS("big_qsize.RDS")

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
ggsave("../size1d2.png",  width = 9, height = 3)

