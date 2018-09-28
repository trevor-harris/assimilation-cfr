rm(list = ls()); gc()

library(ggplot2)
library(dplyr)

getwd()

load("../sim/asym/seed777_sims100.RData")
sims4 = sims

# second run (more reliable since higher iterations)
dir = "../sim/asym/seed177_sims1000/"
files = list.files(dir)
sims3 = matrix(0, length(files), 7)
for(f in 1:length(files)) {
  load(paste0(dir, "sim", f, ".RData"))
  sims3[f,] = out
}

sims3 = as.data.frame(sims3)
names(sims3) = c("Seed", "Sims", "PriorFields", "PostFields", "Dimension", "Range", "Size")


# second run (more reliable since higher iterations)
dir = "../sim/asym/seed1023_sims1000/"
files = list.files(dir)
sims2 = matrix(0, length(files), 7)
for(f in 1:length(files)) {
  load(paste0(dir, "sim", f, ".RData"))
  sims2[f,] = out
}

sims2 = as.data.frame(sims2)
names(sims2) = c("Seed", "Sims", "PriorFields", "PostFields", "Dimension", "Range", "Size")

# first run
dir = "../sim/asym/seed1_sims500/"
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

sims = rbind(sims1, sims2, sims3, sims4)

save_dir = "../Presentation/results/"

sims = sims %>%
  mutate(Range = as.factor(Range),
         Dimension = as.factor(Dimension),
         PriorFields = as.factor(PriorFields),
         PostFields = as.factor(PostFields))

names(sims) = c("Smoothness", "Size", "Seed", "Sims", "XFields", "YFields", "Dimension")
levels(sims$Dimension) = c("10 x 10", "10 x 10", "20 x 20", "20 x 20", "30 x 30", "30 x 30")
sims$XFields = factor(sims$XFields, levels = c("100", "500", "1000"), ordered = TRUE)
sims$Smoothness = factor(sims$Smoothness, levels = c("1", "5", "10"), ordered = TRUE)

# Dim vs Range
ggplot(data = sims, aes(x = Dimension, y = Size, color = Smoothness)) +
  # geom_point() +
  geom_boxplot() +
  geom_abline(intercept = 0.05, slope = 0) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5)) +
  ggtitle("Size v.s. Dimension and Range") +
ggsave(paste0(save_dir, "size1.png"), width = 5, height = 3.2)

# Dim vs Prior
ggplot(data = sims, aes(x = Dimension, y = Size, color = XFields)) +
  # geom_point() +
  geom_boxplot() +
  geom_abline(intercept = 0.05, slope = 0) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5)) +
  ggtitle("Size v.s. Dimension and #X Fields")
ggsave(paste0(save_dir, "size2.png"), width = 5, height = 3.2)

# Prior vs Range
ggplot(data = sims, aes(x = XFields, y = Size, color = Smoothness)) +
  # geom_point() +
  geom_boxplot() +
  geom_abline(intercept = 0.05, slope = 0) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5)) +
  ggtitle("Size v.s.  #X Fields and Range") 
ggsave(paste0(save_dir, "size3.png"), width = 5, height = 3.2)

# Prior vs Range
ggplot(data = sims[sims["XFields"] != "100",], aes(x = XFields, y = Size, color = Smoothness)) +
  # geom_point() +
  geom_boxplot() +
  geom_abline(intercept = 0.05, slope = 0) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5)) +
  ggtitle("Size v.s. #X Fields and Range")
ggsave(paste0(save_dir, "size3b.png"), width = 5, height = 3.2)

# Prior vs Post
ggplot(data = sims, aes(x = XFields, y = Size, color = YFields)) +
  # geom_point() +
  geom_boxplot() +
  geom_abline(intercept = 0.05, slope = 0) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5)) +
  ggtitle("Size v.s. #X Fields and #Y Fields")
ggsave(paste0(save_dir, "size4.png"), width = 5, height = 3.2)


# third run (same as second with new seed)
dir = "../sim/out/"
files = list.files(dir)
sims3 = matrix(0, length(files), 7)
for(f in 1:length(files)) {
  load(paste0(dir, "sim", f, ".RData"))
  sims3[f,] = out
}

sims = as.data.frame(sims3)
names(sims) = c("Seed", "Sims", "PriorFields", "PostFields", "Dimension", "Range", "Size")

save(sims, file = "../sim/asym/seed777_sims100.RData")
