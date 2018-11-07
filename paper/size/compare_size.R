rm(list = ls()); gc()

library(ggplot2)
library(dplyr)
library(reshape2)

# import raw size data
dir = "../temp/size/independent/"
dir = "../temp/size/output/"
files = list.files(dir)
size_data = readRDS(paste0(dir, files[1]))
for(f in 2:length(files)) {
  size_data = rbind(size_data, readRDS(paste0(dir, files[f])))
}

# Compare mean shifts
size = size_data %>%
  select(pval, method, functions, range) %>%
  group_by(method, functions, range) %>%
  summarize(size = mean(pval < 0.05)) %>%
  ungroup() %>%
  mutate(method = recode(method, 
                         K_WA = "K Statistic",
                         Q_XD = "Quality Index"),
         functions = as.factor(functions),
         range = as.factor(range)
  )

ggplot(size, aes(x=method, y=size, color=functions)) +
  geom_boxplot() +
  geom_hline(yintercept = 0.05) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5)) +
  xlab("Method") +
  ylab("Size") +
  ggtitle("Type 1 error under F = G")

ggplot(size, aes(x=method, y=size, color=range)) +
  geom_boxplot() +
  geom_hline(yintercept = 0.05) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5)) +
  xlab("Method") +
  ylab("Size") +
  ggtitle("Type 1 error under F = G")
