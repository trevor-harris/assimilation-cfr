rm(list = ls()); gc()

library(ggplot2)
library(dplyr)
library(reshape2)

# import raw size data
dir = "../temp/power/independent/"
dir = "../research/assimilation-cfr/paper/power/independent/"

files = list.files(dir)
power_data = readRDS(paste0(dir, files[1]))
for(f in 2:length(files)) {
  power_data = rbind(power_data, readRDS(paste0(dir, files[f])))
}

# Compare mean shifts
power_mu = power_data %>%
  filter(feature == "mean") %>%
  select(pval, method, functions, mu1, mu2) %>%
  group_by(method, functions, mu2) %>%
  summarize(power = mean(pval < 0.05)) %>%
  ungroup() %>%
  mutate(method = recode(method, 
                         K = "K Statistic",
                         Q = "Quality Index")
  )
ggplot(power_mu, aes(x=mu2, y=power, color=method)) +
  geom_line() +
  geom_hline(yintercept = 0.05) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5)) +
  xlab("Posterior mean") +
  ylab("Power") +
  ggtitle("Power against mean change")

# Compare Scale changes
power_sd = power_data %>%
  filter(feature == "sd") %>%
  select(pval, method, functions, sd1, sd2) %>%
  group_by(method, functions, sd2) %>%
  summarize(power = mean(pval < 0.05)) %>%
  ungroup() %>%
  mutate(method = recode(method, 
                         K = "K Statistic",
                         Q = "Quality Index")
  )
ggplot(power_sd, aes(x=sd2, y=power, color=method)) +
  geom_line() +
  geom_hline(yintercept = 0.05) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5)) +
  xlab("Posterior scale") +
  ylab("Power") +
  ggtitle("Power against scale change")


# Compare correlation changes
power_corr = power_data %>%
  filter(feature == "corr") %>%
  select(pval, method, functions, r1, r2) %>%
  group_by(method, functions, r2) %>%
  summarize(power = mean(pval < 0.05)) %>%
  ungroup() %>%
  mutate(method = recode(method, 
                         K = "K Statistic",
                         Q = "Quality Index")
  )
ggplot(power_corr, aes(x=r2, y=power, color=method)) +
  geom_line() +
  geom_hline(yintercept = 0.05) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5)) +
  xlab("Posterior correlation") +
  ylab("Power") +
  ggtitle("Power against correlation change")
