rm(list = ls()); gc()

library(ggplot2)
library(dplyr)
library(reshape2)

setwd("/Users/trevh/research")

dir = "/Users/trevh/research/assimilation-cfr/paper/ngphetpower/data/"
files = list.files(dir)
power_data = readRDS(paste0(dir, "/", files[1]))
for(f in 2:length(files)) {
  power_data = rbind(power_data, readRDS(paste0(dir, "/", files[f])))
}

# Compare mean shifts
pow_mu = power_data %>%
  filter(feature == "mean") %>%
  dplyr::select(pval, method, amp) %>%
  group_by(method, amp) %>%
  summarize(power = mean(pval < 0.05, na.rm = T)) %>%
  ungroup() %>%
  mutate(Statistic = recode(method, 
                            Kolm = "KD",
                            Qual = "QI"),
         parameter = "Mean",
         vary = amp
  ) %>%
  dplyr::select(
    Statistic, parameter, vary, power
  )

pow_sd = power_data %>%
  filter(feature == "sd") %>%
  dplyr::select(pval, method, amp) %>%
  group_by(method, amp) %>%
  summarize(power = mean(pval < 0.05, na.rm = T)) %>%
  ungroup() %>%
  mutate(Statistic = recode(method, 
                            Kolm = "KD",
                            Qual = "QI"),
         parameter = "Std. Dev.",
         vary = amp
  ) %>%
  dplyr::select(
    Statistic, parameter, vary, power
  )

power = rbind(pow_mu, pow_sd)
power$parameter = as.factor(power$parameter)
levels(power$parameter) = c(TeX("$\\mu"), TeX("$\\sigma"))

ggplot(power, aes(x=vary, y=power, color=Statistic)) +
  geom_line(size = 1) +
  geom_hline(yintercept = 0.05) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5), 
        text = element_text(size=18)) +
  xlab(TeX("$\\kappa")) +
  ylab("Power") +
  facet_wrap(. ~ parameter, nrow = 1, scales = "free_x",
             labeller = label_parsed)

ggsave("assimilation-cfr/paper/ngphetpower/ngppower_het.png", width = 9, heigh = 3)

