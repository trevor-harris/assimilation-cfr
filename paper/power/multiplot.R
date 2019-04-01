rm(list = ls()); gc()

library(ggplot2)
library(dplyr)
library(reshape2)

setwd("C:/Users/trevorh2/research/")
setwd("../research/assimilation-cfr/paper/power/het/new")

dir = getwd()
files = list.files(dir)
power_data = readRDS(paste0(dir, "/", files[1]))
for(f in 2:length(files)) {
  power_data = rbind(power_data, readRDS(paste0(dir, "/", files[f])))
}

power_data[["Statistic"]] = factor(power_data[["Stat"]], labels = c("K", "QI"))


ggplot(power_data, aes(x=V2, y=V1, color=Statistic)) +
  geom_line() +
  geom_hline(yintercept = 0.05) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5)) +
  xlab("Parameter Change") +
  ylab("Rejection Proportion") +
  facet_wrap(. ~ param, nrow = 1, scales = "free_x")

ggsave("../multi_het.png", width = 9, heigh = 3)


##### SMALLER R1

# Compare mean shifts
power_mu = power_data %>%
  filter(feature == "mean") %>%
  select(pval, method, functions, mu1, mu2, r1) %>%
  group_by(method, functions, mu2, r1) %>%
  summarize(power = mean(pval < 0.05)) %>%
  ungroup() %>%
  mutate(Statistic = recode(method, 
                            K = "K",
                            Q = "QI"),
         parameter = "Mean",
         vary = mu2
  ) %>% 
  filter(
    r1 == 10
  ) %>%
  select(
    method, power, Statistic, parameter, vary
  )


# Compare Scale changes
power_sd = power_data %>%
  filter(feature == "sd") %>%
  select(pval, method, functions, sd1, sd2, r1) %>%
  group_by(method, functions, sd2, r1) %>%
  summarize(power = mean(pval < 0.05)) %>%
  ungroup() %>%
  mutate(Statistic = recode(method, 
                            K = "K",
                            Q = "QI"),
         parameter = "Standard Deviation",
         vary = sd2
  ) %>% 
  filter(
    r1 == 10
  )  %>%
  select(
    method, power, Statistic, parameter, vary
  )


# Compare correlation changes
power_corr = power_data %>%
  filter(feature == "corr") %>%
  select(pval, method, functions, r1, r2) %>%
  group_by(method, functions, r1, r2) %>%
  summarize(power = mean(pval < 0.05)) %>%
  ungroup() %>%
  mutate(Statistic = recode(method, 
                            K = "K",
                            Q = "QI"),
         parameter = "Range",
         vary = r2 - 10
         
  ) %>% 
  filter(
    r1 == 10
  ) %>%
  select(
    method, power, Statistic, parameter, vary
  )


power2d = rbind(power_mu, power_sd, power_corr)

ggplot(power2d, aes(x=vary, y=power, color=Statistic)) +
  geom_line() +
  geom_hline(yintercept = 0.05) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5)) +
  xlab("Parameter Change") +
  ylab("Rejection Proportion") +
  facet_wrap(. ~ parameter, nrow = 1, scales = "free_x")

ggsave("../multi2d2.png", width = 9, heigh = 3)

############# 1D data

setwd("../../power/out1d/")

dir = getwd()
files = list.files(dir)
power_data = readRDS(paste0(dir, "/", files[1]))
for(f in 2:length(files)) {
  power_data = rbind(power_data, readRDS(paste0(dir, "/", files[f])))
}

# Compare mean shifts
power_mu = power_data %>%
  filter(feature == "mean") %>%
  select(pval, method, functions, mu1, mu2, r1) %>%
  group_by(method, functions, mu2, r1) %>%
  summarize(power = mean(pval < 0.05)) %>%
  ungroup() %>%
  mutate(Statistic = recode(method, 
                            K = "K",
                            Q = "Q"),
         parameter = "Mean",
         vary = mu2
  ) %>% 
  filter(
    r1 == 20
  ) %>%
  select(
    method, power, Statistic, parameter, vary
  )


# Compare Scale changes
power_sd = power_data %>%
  filter(feature == "sd") %>%
  select(pval, method, functions, sd1, sd2, r1) %>%
  group_by(method, functions, sd2, r1) %>%
  summarize(power = mean(pval < 0.05)) %>%
  ungroup() %>%
  mutate(Statistic = recode(method, 
                            K = "K",
                            Q = "Q"),
         parameter = "Standard Deviation",
         vary = sd2
  ) %>% 
  filter(
    r1 == 20
  )  %>%
  select(
    method, power, Statistic, parameter, vary
  )


# Compare correlation changes
power_corr = power_data %>%
  filter(feature == "corr") %>%
  select(pval, method, functions, r1, r2) %>%
  group_by(method, functions, r1, r2) %>%
  summarize(power = mean(pval < 0.05)) %>%
  ungroup() %>%
  mutate(Statistic = recode(method, 
                            K = "K",
                            Q = "Q"),
         parameter = "Range",
         vary = r2 - 20
         
  ) %>% 
  filter(
    r1 == 20
  ) %>%
  select(
    method, power, Statistic, parameter, vary
  )


power1d = rbind(power_mu, power_sd, power_corr)

ggplot(power1d, aes(x=vary, y=power, color=Statistic)) +
  geom_line() +
  geom_hline(yintercept = 0.05) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5)) +
  xlab("Parameter Change") +
  ylab("Power") +
  facet_wrap(. ~ parameter, nrow = 1, scales = "free_x")

ggsave("../multi1d.png", width = 9, heigh = 3)

