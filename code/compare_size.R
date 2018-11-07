rm(list = ls()); gc()

library(ggplot2)
library(dplyr)

getwd()

# import raw size data
dir = "../research/assimilation-cfr/paper/size/output/"
files = list.files(dir)
size_data = readRDS(paste0(dir, files[1]))
for(f in 2:length(files)) {
  size_data = rbind(size_data, readRDS(paste0(dir, files[f])))
}

size = size_data %>%
  select(pval, functions, range, method) %>%
  group_by(method, functions, range,) %>%
  summarize(size = mean(pval < 0.05))


ggplot(size, aes(x=method, y=size, color=as.factor(functions))) +
  geom_boxplot() +
  geom_hline(yintercept = 0.05) +
  theme_classic()

ggplot(size, aes(x=method, y=size, color=as.factor(range))) +
  geom_boxplot() +
  geom_hline(yintercept = 0.05) +
  theme_classic()

ggplot(size, aes(x=as.factor(range), y=size, color=as.factor(method))) +
  geom_boxplot() +
  geom_hline(yintercept = 0.05) +
  theme_classic()

