rm(list = ls()); gc()

library(ggplot2)
library(dplyr)
library(reshape2)

setwd("C:/Users/trevorh2/research/")
setwd("../research/assimilation-cfr/paper/power/out2d/")

# setwd("../../power/out2d/")

# import raw size data
# dir = "../temp/power/independent/"
# dir = "research/assimilation-cfr/paper/power/independent1d/"

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
                            Q = "QI")
  )
ggplot(power_mu, aes(x=mu2, y=power, color=Statistic)) +
  geom_line() +
  geom_hline(yintercept = 0.05) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5)) +
  xlab("Mean shift") +
  ylab("Power") +
  facet_wrap(. ~ r1, nrow = 1)
# ggtitle("Power against mean change")
ggsave("../mu2d.png", width = 9)

# Compare Scale changes
power_sd = power_data %>%
  filter(feature == "sd") %>%
  select(pval, method, functions, sd1, sd2, r1) %>%
  group_by(method, functions, sd2, r1) %>%
  summarize(power = mean(pval < 0.05)) %>%
  ungroup() %>%
  mutate(Statistic = recode(method, 
                            K = "K",
                            Q = "QI")
  )
ggplot(power_sd, aes(x=sd2, y=power, color=Statistic)) +
  geom_line() +
  geom_hline(yintercept = 0.05) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5)) +
  xlab("Std. dev. multiple") +
  ylab("Power") +
  facet_wrap(. ~ r1, nrow = 1)
# ggtitle("Power against scale change")
ggsave("../sd2d.png", width = 9)


# Compare correlation changes
power_corr = power_data %>%
  filter(feature == "corr") %>%
  select(pval, method, functions, r1, r2) %>%
  group_by(method, functions, r1, r2) %>%
  summarize(power = mean(pval < 0.05)) %>%
  ungroup() %>%
  mutate(Statistic = recode(method, 
                            K = "K",
                            Q = "QI")
  )
ggplot(power_corr, aes(x=r2, y=power, color=Statistic)) +
  geom_line() +
  geom_hline(yintercept = 0.05) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5)) +
  xlab("Range Shift") +
  ylab("Power") +
  facet_wrap(. ~ r1, nrow = 1)
# ggtitle("Power against correlation change")
ggsave("../cr2d.png", width = 9, height = 3)

