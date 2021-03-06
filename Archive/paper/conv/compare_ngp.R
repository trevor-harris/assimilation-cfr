rm(list = ls()); gc()

library(ggplot2)
library(dplyr)
library(reshape2)
library(latex2exp)

# import raw size data
dir = "../research/matern/nongp/"
files = list.files(dir)
conv_data = readRDS(paste0(dir, files[1]))
for(f in 2:length(files)) {
  conv_data = rbind(conv_data, readRDS(paste0(dir, files[f])))
}

conv = conv_data %>% na.omit() %>%
  mutate(range = as.factor(range),
         n = as.factor(n1),
         nu = as.factor(nu))


conv = conv %>%
  select(-one_of(c("n1", "n2"))) %>%
  mutate(r = as.factor(range),
         nu = as.factor(nu)
  )

nu_names = c(`0.5` = TeX("$\\nu = 0.5"), `1` = TeX("$\\nu = 1.0"), `1.5` = TeX("$\\nu = 1.5"))

ggplot(conv, aes(x = n, y = sqrt(cdf_diff), fill = r)) +
  geom_boxplot() +
  theme_bw() +
  theme(axis.title = element_text(size=16),
        strip.text = element_text(size=14),
        legend.title = element_text(size=16),
        legend.text = element_text(size=14)) +
  theme(plot.title = element_text(hjust = 0.5)) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  xlab("Number of functions (n)") +
  ylab(TeX("$L^2$ Distance")) +
  facet_grid(nu ~ r, 
             switch = "y", 
             labeller = label_both) +
  theme(strip.text.y = element_text(angle = 180)) +
  guides(fill=guide_legend("Range"))
ggsave("../research/proxy/paper/conv/ngp_l2.png", width = 10, height = 5.3)

crit = conv %>% 
  select(n, r, nu, cval_90, cval_95, cval_99) %>%
  melt(id.vars = c("n", "r", "nu"), variable.name = "critical") %>%
  mutate(critical = recode(critical,
                           "cval_90" = "0.90",
                           "cval_95" = "0.95",
                           "cval_99" = "0.99"))

ggplot(crit, aes(x = n, y = value, fill = critical)) +
  geom_boxplot() +
  theme_bw() +
  theme(axis.title = element_text(size=16),
        strip.text = element_text(size=14),
        legend.title = element_text(size=16),
        legend.text = element_text(size=14)) +
  theme(plot.title = element_text(hjust = 0.5)) + 
  geom_hline(yintercept = 0) +
  xlab("Number of functions (n)") +
  ylab("Critical value difference") +
  ylim(c(-0.06, 0.06)) +
  facet_grid(nu ~ r, 
             switch = "y", 
             labeller = label_both) +
  theme(strip.text.y = element_text(angle = 180)) +
  guides(fill=guide_legend("Level"))
ggsave("../research/proxy/paper/conv/ngp_cv.png", width = 10, height = 5.3)
