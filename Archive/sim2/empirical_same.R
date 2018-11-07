rm(list = ls()); gc()

library(ggplot2)
library(dplyr)

getwd()

source("research/assimilation-cfr/sim/reference.R")

# second run (more reliable since higher iterations)
dir = "research/assimilation-cfr/asymptotic/out/"
files = list.files(dir)

sims = data.frame()
for(f in files) {
  load(paste0(dir, f))
  sims = rbind(sims, out)
}

t = seq(0, 2, length.out = 5000)
perm.cdf = sapply(t, function(x) mean(sims$perms <= x))
theo.cdf = sapply(t, function(x) (ks_cdf(x)))

cdf.gg = data.frame(Method = rep(c("Permutation", "Kolmogorov"), each=length(t)),
                    Value = c(t, t),
                    CDF = c(perm.cdf, theo.cdf))

ggplot(cdf.gg, aes(Value, CDF, color=Method)) +
  geom_line() +
  geom_vline(xintercept = t[min(which(perm.cdf > 0.95))], color="#00BFC4") +
  geom_vline(xintercept = t[min(which(theo.cdf > 0.95))], color="#F8766D") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5)) +
  ggtitle("Kolmogorov Distribution v.s. Permutation Distribution") +
  xlab("K Value") +
  ylab("Probability")
ggsave(paste0("research/assimilation-cfr/paper/misc/", "distribution.png"), width = 5, height = 3.2)


t[min(which(perm.cdf > 0.95))] - t[min(which(theo.cdf > 0.95))]
