rm(list = ls())
gc()



########### READ ME #############

# you must change the working directory to be the submitted_code folder
# none of this will work otherwise
# mine is left here as an example

########## Example
# setwd("/Users/trevh/research/assimilation-cfr/submitted_code/")

#################################




library(ggplot2)
library(dplyr)
library(reshape2)
library(xtable)


# import raw size data
dir = "size/out_gp/"

files = list.files(dir)
size_data = readRDS(paste0(dir, files[1]))
for(f in 2:length(files)) {
  size_data = rbind(size_data, readRDS(paste0(dir, files[f])))
}

# calculate size and summarize data by sample sizes, range, and smoothness
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



### print the full table
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

# latex formatted (in the paper)
xtable(size.tab, booktabs = T)


