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
library(latex2exp)



dir = "power/out_het_gp"
files = list.files(dir)
power_data = readRDS(paste0(dir, "/", files[1]))
for(f in 2:length(files)) {
  power_data = rbind(power_data, readRDS(paste0(dir, "/", files[f])))
}


#  convert mean shift power comparison into ggplot2 compliant format
pow_mu = power_data %>%
  filter(feature == "mean") %>%
  dplyr::select(pval, method, amp) %>%
  group_by(method, amp) %>%
  summarize(power = mean(pval < 0.05)) %>%
  ungroup() %>%
  mutate(Statistic = recode(method, 
                            K = "KD",
                            Q = "QI"),
         parameter = "Mean",
         vary = amp
  ) %>%
  dplyr::select(
    Statistic, parameter, vary, power
  )


#  convert standard deviation power comparison into ggplot2 compliant format
pow_sd = power_data %>%
  filter(feature == "sd") %>%
  dplyr::select(pval, method, amp) %>%
  group_by(method, amp) %>%
  summarize(power = mean(pval < 0.05)) %>%
  ungroup() %>%
  mutate(Statistic = recode(method, 
                            K = "KD",
                            Q = "QI"),
         parameter = "Std. Dev.",
         vary = amp
  ) %>%
  dplyr::select(
    Statistic, parameter, vary, power
  )

# combine back together and plot side by side
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

# ggsave("power/power_het_gp.png", width = 10, heigh = 3)
