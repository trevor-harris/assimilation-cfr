rm(list = ls()); gc()

library(ggplot2)
library(dplyr)
library(reshape2)
library(xtable)

# set to the top level folder
setwd("/Users/trevorh2/research/assimilation-cfr/submit/")

# import raw size data
dir = "size/out_gp/"

files = list.files(dir)
size_data = readRDS(paste0(dir, files[1]))
for(f in 2:length(files)) {
  size_data = rbind(size_data, readRDS(paste0(dir, files[f])))
}

# calculate size and summarize data
size = size_data %>%
  dplyr::select(pval, method, n1, n2, rng, nu) %>%
  group_by(method, n1, n2, rng, nu) %>%
  summarize(size = mean(pval < 0.05)) %>%
  ungroup() %>%
  mutate(stat = as.factor(method),
         n1 = as.factor(n1),
         n2 = as.factor(n2),
         range = as.factor(rng),
         nu = as.factor(nu),
         stat = recode(method, Kolm = "KD", Qual = "QI")
  ) %>% 
  dplyr::select(
    stat, n1, n2, range, nu, size
  ) %>%
  filter(
    stat %in% c("KD", "QI")
  )



### full table
size.tab = size %>%
  mutate(
    size = sprintf("%.2f", round(size, 2)),
    size = ifelse(stat == "QI",
                  paste0("(", size, ")"),
                  size)
  ) %>%
  dcast(
    n1 + n2 + stat ~ nu + range, value.var = "size"
  )
xtable(size.tab, booktabs = T)


