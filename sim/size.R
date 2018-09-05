rm(list = ls()); gc()

library(ggplot2)
library(dplyr)

# third run (same as second with new seed)
dir = "sim/out/"
files = list.files(dir)
sims3 = matrix(0, length(files), 7)
for(f in 1:length(files)) {
  load(paste0(dir, "sim", f, ".RData"))
  sims3[f,] = out
}

sims3 = as.data.frame(sims3)
names(sims3) = c("Seed", "Sims", "PriorFields", "PostFields", "Dimension", "Range", "Size")


# second run (more reliable since higher iterations)
dir = "sim/seed1023_sims1000/"
files = list.files(dir)
sims2 = matrix(0, length(files), 7)
for(f in 1:length(files)) {
  load(paste0(dir, "sim", f, ".RData"))
  sims2[f,] = out
}

sims2 = as.data.frame(sims2)
names(sims2) = c("Seed", "Sims", "PriorFields", "PostFields", "Dimension", "Range", "Size")

# first run
dir = "sim/seed1_sims500/"
files = list.files(dir)
sims1 = matrix(0, length(files), 5)
for(f in 1:length(files)) {
  load(paste0(dir, "sim", f, ".RData"))
  sims1[f,] = out
}
sims1 = as.data.frame(sims1)
names(sims1) = c("PriorFields", "PostFields", "Dimension", "Range", "Size")
sims1[["Seed"]] = 1
sims1[["Sims"]] = 500
sims1 = sims2[c(6, 7, 1, 2, 3, 4, 5)]

sims = rbind(sims1, sims2, sims3)

# write.table(sims, "research/assimilation-cfr/sim/summary", sep="\t")

sims = sims %>%
  mutate(Range = as.factor(Range),
         Dimension = as.factor(Dimension),
         PriorFields = as.factor(PriorFields),
         PostFields = as.factor(PostFields))

levels(sims$Dimension) = c("10 x 10", "20 x 20", "30 x 30")

# Dim vs Range
ggplot(data = sims, aes(x = Dimension, y = Size, color = Range)) +
  # geom_point() +
  geom_boxplot() +
  geom_abline(intercept = 0.05, slope = 0) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5)) +
  ggtitle("Size v.s. Dimension and Range") 

# Dim vs Prior
ggplot(data = sims, aes(x = Dimension, y = Size, color = PriorFields)) +
  # geom_point() +
  geom_boxplot() +
  geom_abline(intercept = 0.05, slope = 0) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5)) +
  ggtitle("Size v.s. Dimension and no. Prior Fields")

# Prior vs Range
ggplot(data = sims, aes(x = PriorFields, y = Size, color = Range)) +
  # geom_point() +
  geom_boxplot() +
  geom_abline(intercept = 0.05, slope = 0) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5)) +
  ggtitle("Size v.s. no. Prior Fields and Range") 

# Prior vs Post
ggplot(data = sims, aes(x = PriorFields, y = Size, color = PostFields)) +
  # geom_point() +
  geom_boxplot() +
  geom_abline(intercept = 0.05, slope = 0) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5)) +
  ggtitle("Size v.s. no. Prior Fields and no. Posterior Fields") 


