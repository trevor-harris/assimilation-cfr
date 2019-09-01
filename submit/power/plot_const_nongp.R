rm(list = ls()); gc()

library(ggplot2)
library(dplyr)
library(reshape2)
library(latex2exp)


dir = "power/out_const_nongp/"
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
  summarize(power = mean(pval < 0.05, na.rm = T)) %>%
  ungroup() %>%
  mutate(Statistic = as.factor(method),
         Statistic = recode(method, 
                            Kolm = "KD",
                            Qual = "QI"),
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
  summarize(power = mean(pval < 0.05, na.rm = T)) %>%
  ungroup() %>%
  mutate(Statistic = as.factor(method),
         Statistic = recode(method, 
                            Kolm = "KD",
                            Qual = "QI"),
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
  summarize(power = mean(pval < 0.05, na.rm = T)) %>%
  ungroup() %>%
  mutate(Statistic = as.factor(method),
         Statistic = recode(method, 
                            Kolm = "KD",
                            Qual = "QI"),
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
  summarize(power = mean(pval < 0.05, na.rm = T)) %>%
  ungroup() %>%
  mutate(Statistic = as.factor(method),
         Statistic = recode(method, 
                            Kolm = "KD",
                            Qual = "QI"),
         parameter = "Smoothness",
         vary = nu2
  ) %>%
  select(
    Statistic, parameter, vary, power
  )

power = rbind(pow_mu, pow_sd, pow_cor, pow_sm)
power$parameter = as.factor(power$parameter)

power = power[power$Statistic %in% c("KD", "QI"),]
power$Statistic = factor(power$Statistic, labels = c("KD", "QI"))

levels(power$parameter) = c(TeX("$\\mu"), TeX("$r"), TeX("$\\nu"), TeX("$\\sigma"))

ggplot(power, aes(x=vary, y=power, color=Statistic)) +
  geom_line(size = 1) +
  geom_hline(yintercept = 0.05) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5), 
        text = element_text(size=18)) +
  xlab("Parameter Change") +
  ylab("Power") +
  facet_wrap(. ~ parameter, nrow = 1, scales = "free_x",
             labeller = label_parsed)

ggsave("power/power_const_nongp.png", width = 10, heigh = 3)
