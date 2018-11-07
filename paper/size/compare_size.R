rm(list = ls()); gc()

library(ggplot2)
library(dplyr)
library(reshape2)

getwd()

# import raw size data
dir = "../research/assimilation-cfr/paper/size/output_rev/"
files = list.files(dir)
size_data = readRDS(paste0(dir, files[1]))
for(f in 2:length(files)) {
  size_data = rbind(size_data, readRDS(paste0(dir, files[f])))
}

size = size_data %>%
  select(pval, functions, range, method) %>%
  group_by(method, functions, range) %>%
  summarize(size = mean(pval < 0.05)) %>%
  ungroup() %>%
  mutate(method = recode(method, 
                         K_WA = "K Statistic",
                          Q_XD = "Quality Index")
  )


ggplot(size, aes(x=method, y=size, color=method)) +
  geom_boxplot() +
  geom_hline(yintercept = 0.05) +
  theme_classic() +
  xlab("Method") +
  ylab("Size")

ggplot(size, aes(x=method, y=size, color=as.factor(functions))) +
  geom_boxplot() +
  geom_hline(yintercept = 0.05) +
  theme_classic() +
  xlab("Method") +
  ylab("Size")

ggplot(size, aes(x=method, y=size, color=as.factor(range))) +
  geom_boxplot() +
  geom_hline(yintercept = 0.05) +
  theme_classic()

ggplot(size, aes(x=as.factor(range), y=size, color=as.factor(method))) +
  geom_boxplot() +
  geom_hline(yintercept = 0.05) +
  theme_classic()


size = size_data %>%
  filter(method == "K_WA", functions != 20) %>%
  select(range, functions, pval)

# little helper
set.seed(1023)
booty = function(size, rng, funcs) {
  pvals = size %>% 
    filter(range == rng, functions == funcs) %>% 
    select(pval) %>% 
    unlist() %>% 
    as.numeric()
  boots = sapply(1:2000, function(x) mean(sample(pvals, replace = T) < 0.05))
  data.frame(range = rng, functions = funcs, reps = boots)
}

boot = data.frame(range = numeric(0), functions = numeric(0), reps = numeric(0))
for(r in unique(size$range)) {
  for(f in unique(size$functions)) {
    boot = rbind(boot, booty(size, r, f))
  }
}
boot[["functions"]] = as.factor(boot[["functions"]])
boot[["range"]] = as.factor(boot[["range"]])

ggplot(boot, aes(range, reps, color=functions)) +
  geom_boxplot() + 
  geom_hline(yintercept = 0.05) +
  theme_classic()

ggplot(boot, aes(functions, reps, color=range)) +
  geom_boxplot() + 
  geom_hline(yintercept = 0.05) +
  theme_classic()
