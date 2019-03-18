rm(list = ls()); gc()

library(reshape2)
library(ggplot2)
library(dplyr)

setwd("../../conv/out/")

conv = readRDS("conv.RDS")
conv = conv[conv$range < 25, ]
conv2 = melt(conv[,-1], id.vars = c("nfun", "range"), variable.name = "Critical")

conv[["Range"]] = as.factor(conv[["range"]])
conv[["nfun"]] = as.factor(conv[["nfun"]])

conv[["mse"]] = conv[["cdf_diff"]]
conv[["rmse"]] = sqrt(conv[["cdf_diff"]])

# summary(conv2$value)
ggplot(conv, aes(x = nfun, y = rmse, fill = Range)) +
  geom_boxplot() +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5)) + 
  geom_hline(yintercept = 0) +
  xlab("Number of functions (n)") +
  ylab("RMSE") +
  facet_wrap(. ~ Range, nrow = 1)
  # ylim(c(0, 0.03))
ggsave("../conv1.png", width = 9, height = 3)

conv2[["Range"]] = as.factor(conv2[["range"]])
conv2[["nfun"]] = as.factor(conv2[["nfun"]])
conv2[["Critical2"]] = factor(conv2[["Critical"]], labels =  c("0.90", "0.95", "0.99"))

ggplot(conv2, aes(x = nfun, y = value, fill = Critical2)) +
  geom_boxplot() +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5)) + 
  geom_hline(yintercept = 0) +
  xlab("Number of functions (n)") +
  ylab("Difference in critical values") +
  ylim(c(-0.018, 0.004)) +
  facet_wrap(. ~ Range, nrow = 1) +
  guides(fill=guide_legend("Level"))
ggsave("../conv2.png", width = 9, height = 3)

