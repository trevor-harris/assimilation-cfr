rm(list = ls()); gc()

library(ggplot2)
library(dplyr)
library(reshape2)


dir = "../research/assimilation-cfr/paper/power/const/"
files = list.files(dir)
power_data = readRDS(paste0(dir, "/", files[1]))
for(f in 2:length(files)) {
  power_data = rbind(power_data, readRDS(paste0(dir, "/", files[f])))
}

# Compare mean shifts
pow_mu = power_data %>%
  filter(feature == "mean") %>%
  select(pval, method, mu1, mu2) %>%
  group_by(method, mu2) %>%
  summarize(power = mean(pval < 0.05)) %>%
  ungroup() %>%
  mutate(Statistic = recode(method, 
                            K = "K",
                            Q = "QI"),
         parameter = "Mean",
         vary = mu2
  ) %>%
  select(
    Statistic, parameter, vary, power
  )

pow_sd = power_data %>%
  filter(feature == "sd") %>%
  select(pval, method, sd1, sd2) %>%
  group_by(method, sd2) %>%
  summarize(power = mean(pval < 0.05)) %>%
  ungroup() %>%
  mutate(Statistic = recode(method, 
                            K = "K",
                            Q = "QI"),
         parameter = "Std. Dev.",
         vary = sd2
  ) %>%
  select(
    Statistic, parameter, vary, power
  )


pow_cor = power_data %>%
  filter(feature == "corr") %>%
  select(pval, method, r1, r2) %>%
  group_by(method, r2) %>%
  summarize(power = mean(pval < 0.05)) %>%
  ungroup() %>%
  mutate(Statistic = recode(method, 
                            K = "K",
                            Q = "QI"),
         parameter = "Range",
         vary = r2
  ) %>%
  select(
    Statistic, parameter, vary, power
  )


pow_sm = power_data %>%
  filter(feature == "smooth") %>%
  select(pval, method, nu1, nu2) %>%
  group_by(method, nu2) %>%
  summarize(power = mean(pval < 0.05)) %>%
  ungroup() %>%
  mutate(Statistic = recode(method, 
                            K = "K",
                            Q = "QI"),
         parameter = "Smoothness",
         vary = nu2
  ) %>%
  select(
    Statistic, parameter, vary, power
  )

power = rbind(pow_mu, pow_sd, pow_cor, pow_sm)

ggplot(power, aes(x=vary, y=power, color=Statistic)) +
  geom_line() +
  geom_hline(yintercept = 0.05) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5)) +
  xlab("Parameter Change") +
  ylab("Rejection Proportion") +
  facet_wrap(. ~ parameter, nrow = 1, scales = "free_x")

ggsave("../research/assimilation-cfr/paper/power/power_const.png", width = 9, heigh = 3)
