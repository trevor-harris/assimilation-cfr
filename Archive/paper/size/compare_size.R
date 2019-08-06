rm(list = ls()); gc()

library(ggplot2)
library(dplyr)
library(reshape2)
library(xtable)

# import raw size data
dir = "../research/assimilation-cfr/paper/size/out/"
files = list.files(dir)
size_data = readRDS(paste0(dir, files[1]))
for(f in 2:length(files)) {
  size_data = rbind(size_data, readRDS(paste0(dir, files[f])))
}

size = size_data %>%
  select(pval, method, n1, n2, rng, nu) %>%
  group_by(method, n1, n2, rng, nu) %>%
  summarize(size = mean(pval < 0.05)) %>%
  ungroup() %>%
  mutate(stat = recode(method, 
                            K = "K",
                            Q = "QI"),
         n1 = as.factor(n1),
         n2 = as.factor(n2),
         range = as.factor(rng),
         nu = as.factor(nu)
  ) %>% 
  select(
    stat, n1, n2, range, nu, size
  )

### long
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


### wide
size.tab = size %>%
  mutate(
    size = as.character(round(size, 2)),
    size = ifelse(stat == "QI",
                  paste0("(", size, ")"),
                  size)
  ) %>%
  dcast(
    n1 + n2 + nu + range ~ stat, value.var = "size"
  ) %>%
  mutate(
    size = paste0(K, " ", QI)
  ) %>%
  select(
    n1, n2, nu, range, size
  ) %>%
  dcast(
    nu + n1 + n2 ~ range, value.var = "size"
  )
xtable(size.tab, booktabs = T)


### multi
size.tab = size %>%
  mutate(
    size = sprintf("%.2f" ,  round(size, 2)),
    size = ifelse(stat == "QI",
                  paste0("(", size, ")"),
                  size)
  ) %>%
  dcast(
    n1 + n2 + nu + range ~ stat, value.var = "size"
  ) %>%
  mutate(
    size = paste0(K, " ", QI)
  ) %>%
  select(
    n1, n2, nu, range, size
  ) %>%
  dcast(
    n1 + n2 + nu ~ range, value.var = "size"
  )
xtable(size.tab %>% 
         filter(nu == "0.5") %>%
         select(-nu), booktabs = T)
xtable(size.tab %>% 
         filter(nu == "1") %>%
         select(-nu), booktabs = T)
xtable(size.tab %>% 
         filter(nu == "1.5") %>%
         select(-nu), booktabs = T)

size.tab = size %>%
  filter(
    range == 0.4,
    nu == 1
  ) %>%
  select(
    stat, n1, n2, size
  ) %>%
  dcast(
   n1 + stat ~ n2, value.var = "size"
  )
xtable(size.tab, booktabs = T)


size.tab = size %>%
  filter(n1 == n2) %>%
  dcast(
     range + stat ~ nu + n1, value.var = "size"
  )


size.tab = size %>%
  filter(n1 == n2) %>%
  dcast(
    stat + nu + range ~ n1, value.var = "size"
  )
xtable(size.tab, booktabs = T)
>>>>>>> e2d04c7a935436b7d6ee9a61b920dad510b5ad5a

