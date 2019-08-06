# rm(list = ls()); gc()

library(ggplot2)
library(dplyr)
library(reshape2)

setwd("C:/Users/trevorh2/research/")
setwd("../research/assimilation-cfr/paper/power/outstoch2/")

dir = getwd()
files = list.files(dir)
power_data = readRDS(paste0(dir, "/", files[1]))

for(f in 2:length(files)) {
  power_data = rbind(power_data, readRDS(paste0(dir, "/", files[f])))
}


ggplot(power_data, aes(x=V2, y=V1, color=Stat)) +
  geom_line() +
  geom_hline(yintercept = 0.05) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5)) +
  xlab("Overall Parameter Change") +
  ylab("Rejection Proportion") +
  facet_wrap(. ~ param, nrow = 1, scales = "free_x")

# ggsave("../stoch.png", width = 9, heigh = 3)
