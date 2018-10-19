rm(list = ls()); gc()

library(ggplot2)
library(dplyr)

getwd()

# second run (more reliable since higher iterations)
dir = "../paper/size/out/"
files = list.files(dir)
sims = matrix(0, length(files), 7)
for(f in files) {
  load(paste0(dir, f))
  sims[which(files == f),] = out
}

sims = as.data.frame(sims)
names(sims) = c("Seed", "Sims", "PriorFields", "PostFields", "Dimension", "Range", "Size")

save_dir = "../paper/size/"

sims = sims %>%
  mutate(Range = as.factor(Range),
         Dimension = as.factor(Dimension),
         PriorFields = as.factor(PriorFields),
         PostFields = as.factor(PostFields))

names(sims) = c("Seed", "Sims", "XFields", "YFields","Dimension", "Smoothness", "Size")
levels(sims$Dimension) = c("10 x 10","20 x 20", "30 x 30")
sims$XFields = factor(sims$XFields, levels = c("100", "500", "1000"), ordered = TRUE)
sims$Smoothness = factor(sims$Smoothness, levels = c("1", "5", "10"), ordered = TRUE)

# Dim vs Range
ggplot(data = sims, aes(x = Dimension, y = Size, color = Smoothness)) +
  # geom_point() +
  geom_boxplot() +
  geom_abline(intercept = 0.05, slope = 0) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5)) +
  ggtitle("Size v.s. Dimension and Range")
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
  ggtitle("Size by Number of fields in X and Smoothness") +
  xlab("Number of fields in X")
ggsave(paste0(save_dir, "size3.png"), width = 5, height = 3.2)

# Prior vs Range
ggplot(data = sims[sims["XFields"] != "100",], aes(x = XFields, y = Size, color = Smoothness)) +
  # geom_point() +
  geom_boxplot() +
  geom_abline(intercept = 0.05, slope = 0) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5)) +
  ggtitle("Size by Number of fields in X and Smoothness") +
  xlab("Number of fields in X")
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

