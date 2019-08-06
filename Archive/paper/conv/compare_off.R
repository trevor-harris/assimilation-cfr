rm(list = ls()); gc()

library(ggplot2)
library(dplyr)
library(reshape2)

# import raw size data
dir = "../research/assimilation-cfr/paper/conv/off/"
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
  mutate(range = recode(range,
                        "0.2" = "r = 0.2",
                        "0.3" = "r = 0.3",
                        "0.4" = "r = 0.4",
                        "0.5" = "r = 0.5"
  )
  )

ggplot(conv, aes(x = n, y = sqrt(cdf_diff), fill = range)) +
  geom_boxplot() +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5)) + 
  geom_hline(yintercept = 0) +
  xlab("Number of functions (n)") +
  ylab(bquote(''*L^2*'  Distance')) +
  facet_wrap(nu ~ range, nrow = 3) +
  guides(fill=guide_legend("Range"))
ggsave("../research/assimilation-cfr/paper/conv/off_l2.png")


crit = conv %>% 
  select( -cdf_diff) %>%
  melt(id.vars = c("n", "range", "nu"), variable.name = "critical") %>%
  mutate(critical = recode(critical,
                           "cval_90" = "0.90",
                           "cval_95" = "0.95",
                           "cval_99" = "0.99"))

ggplot(crit, aes(x = n, y = value, fill = critical)) +
  geom_boxplot() +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5)) + 
  geom_hline(yintercept = 0) +
  xlab("Number of functions (n)") +
  ylab("Kolmogorov - Permutation critical values") +
  ylim(c(-0.04, 0.04)) +
  facet_wrap(nu ~ range, nrow = 3) +
  guides(fill=guide_legend("Level"))
ggsave("../research/assimilation-cfr/paper/conv/off_cv.png")



# smaller plots 
conv.r = conv %>% 
  filter(nu == 0.5)

ggplot(conv.r, aes(x = n, y = sqrt(cdf_diff), fill = range)) +
  geom_boxplot() +
  theme_classic() +
  ggtitle("Increasing range") +
  geom_hline(yintercept = 0) +
  xlab("Replicates (n)") +
  ylab(bquote(''*L^2*'  Distance')) +
  facet_wrap(. ~ range, nrow = 1)


conv.nu = conv %>% 
  filter(range == "r = 0.2")

ggplot(conv.nu, aes(x = n, y = sqrt(cdf_diff), fill = nu)) +
  geom_boxplot() +
  theme_classic() +
  ggtitle("Increasing smoothness") +
  geom_hline(yintercept = 0) +
  xlab("Replicates (n)") +
  ylab(bquote(''*L^2*'  Distance')) +
  facet_wrap(. ~ nu, nrow = 1)





