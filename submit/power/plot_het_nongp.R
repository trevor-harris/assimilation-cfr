rm(list = ls()); gc()

library(ggplot2)
library(dplyr)
library(reshape2)


dir = "power/out_het_nongp"
files = list.files(dir)
power_data = readRDS(paste0(dir, "/", files[1]))
for(f in 2:length(files)) {
  power_data = rbind(power_data, readRDS(paste0(dir, "/", files[f])))
}

# Compare mean shifts
pow_mu = power_data %>%
  filter(feature == "mean") %>%
  select(pval, method, amp) %>%
  group_by(method, amp) %>%
  summarize(power = mean(pval < 0.05)) %>%
  ungroup() %>%
  mutate(Statistic = recode(method, 
                            Kolm = "KD",
                            Qual = "QI"),
         parameter = "Mean",
         vary = amp
  ) %>%
  select(
    Statistic, parameter, vary, power
  )

pow_sd = power_data %>%
  filter(feature == "sd") %>%
  select(pval, method, amp) %>%
  group_by(method, amp) %>%
  summarize(power = mean(pval < 0.05)) %>%
  ungroup() %>%
  mutate(Statistic = recode(method, 
                            Kolm = "KD",
                            Qual = "QI"),
         parameter = "Std. Dev.",
         vary = amp
  ) %>%
  select(
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

ggsave("power/power_het_nongp.png", width = 10, heigh = 3)
