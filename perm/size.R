rm(list = ls()); gc()

library(ggplot2)
library(dplyr)

getwd()

# second run (more reliable since higher iterations)
dir = "research/assimilation-cfr/perm/out/"
files = list.files(dir)
sims = matrix(0, length(files), 9)
for(f in files) {
  load(paste0(dir, f))
  sims[which(files == f),] = as.numeric(out[1,])
}

sims = as.data.frame(sims)
names(sims) = c("Seed", "i", "nfun", "mean", "sd", "range", "nperm", "empirical", "theoretical")

empiric = sapply(1:100, function(x) mean(sample(sims$empirical, replace = T)))
theory = sapply(1:100, function(x) mean(sample(sims$theoretical, replace = T)))

size = data.frame(method = rep(c("Permutation", "Theoretical"), each = 100),
                  Size = c(empiric, theory))

ggplot(size, aes(x = method, y = Size, color = method)) +
  geom_boxplot() +
  geom_abline(intercept = 0.05, slope = 0) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5)) +
  ggtitle("Size")
