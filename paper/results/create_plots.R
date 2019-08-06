rm(list = ls()); gc()

library(ggplot2)
library(dplyr)
library(reshape2)
library(latex2exp)


dir = "../research/proxy/results/data/"
files = list.files(dir)
cl_data = readRDS(paste0(dir, "/", files[1]))
for(f in 3:length(files)) {
  cl_data = rbind(cl_data, readRDS(paste0(dir, "/", files[f])))
}
cl_data = rbind(cl_data, readRDS(paste0(dir, "/", files[2])))
cl_data[["Year"]] = 850:(850+997)
cl_data[["stat"]] = cl_data[["stat"]] * sqrt(100 + 100) / 100

ggplot(cl_data, aes(x=Year, y=stat)) +
  geom_point(size = 2) +
  geom_smooth(color = "red", size = 2, se = F) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5), 
        text = element_text(size=24)) +
  ylab("KD") 

ggsave("../research/proxy/results/effect_over_time.png", width = 8, heigh = 6)


ggplot(cl_data, aes(x=Year, y=pval)) +
  geom_point(size = 2) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5), 
        text = element_text(size=24)) +
  ylab("p-values") 

ggsave("../research/proxy/results/pval_over_time.png", width = 8, heigh = 6)
