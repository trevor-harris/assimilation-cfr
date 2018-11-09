rm(list = ls()); gc()

library(ggplot2)
library(dplyr)
library(reshape2)

read_era = function(dir, file) {
  cbind(readRDS(paste0(dir, file)), era = as.numeric(strsplit(file, "\\D+")[[1]][-1]))
}

# import raw size data
dir = "../temp/power/independent/"
dir = "research/assimilation-cfr/paper/results/results/"
files = list.files(dir)

temperature = read_era(dir, files[1])
for(f in 3:length(files)) {
  temperature = rbind(temperature, read_era(dir, files[f]))
}
temperature = rbind(temperature, read_era(dir, files[2]))
temperature[["time"]] = 1:998

temperature$stat = temperature$stat / (sqrt((100*100)/(100 + 100)))

ggplot(temperature, aes(time, stat)) +
  geom_point(alpha = 0.7) +
  geom_smooth(se = F, color = "red") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5)) +
  ylab("Effect Size") +
  xlab("Time") +
  ggtitle("K(F, G) over time")
ggsave("research/assimilation-cfr/paper/results/effect_over_time.png", width = 5, height = 3.2)

