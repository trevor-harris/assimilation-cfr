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


# import raw size data
dir = "conv/out_gp/"
files = list.files(dir)
conv_data = readRDS(paste0(dir, files[1]))
for(f in 2:length(files)) {
  conv_data = rbind(conv_data, readRDS(paste0(dir, files[f])))
}

# convert to ggplot2 compliant format and subset to only equal sample sizes
conv = conv_data %>% na.omit() %>%
  mutate(range = as.factor(range),
         n = as.factor(n1),
         nu = as.factor(nu))

conv = conv %>%
  dplyr::select(-one_of(c("n1", "n2"))) %>%
  mutate(r = as.factor(range),
         nu = as.factor(nu)
  )

# use greek characters for labelling
levels(conv$nu) = c(TeX("$\\nu = 0.5"), TeX("$\\nu = 1.0"), TeX("$\\nu = 1.5"))
levels(conv$r) = c(TeX("$r = 0.2"), TeX("$r = 0.3"), TeX("$r = 0.4"), TeX("$r = 0.5"))

# plot the L2 distance metric
ggplot(conv, aes(x = n, y = sqrt(cdf_diff), fill = range)) +
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
             labeller = label_parsed) +
  guides(fill=guide_legend("Range"))
# ggsave("conv/gp_l2.png", width = 9.1, height = 5)

# pull out only the critical value results
crit = conv %>% 
  dplyr::select(n, r, nu, cval_90, cval_95, cval_99) %>%
  melt(id.vars = c("n", "r", "nu"), variable.name = "critical") %>%
  mutate(critical = recode(critical,
                           "cval_90" = "0.90",
                           "cval_95" = "0.95",
                           "cval_99" = "0.99"))

# plot critical value differences
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
  ylim(c(-0.04, 0.04)) +
  facet_grid(nu ~ r, 
             switch = "y",
             labeller = label_parsed) +
  guides(fill=guide_legend("Level"))
# ggsave("conv/gp_cv.png", width = 9.1, height = 5)
