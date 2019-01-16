rm(list = ls()); gc()

library(ggplot2)
library(dplyr)
library(reshape2)

# import raw size data
dir = "research/assimilation-cfr/paper/size/independent/"
files = list.files(dir)
size_data = readRDS(paste0(dir, files[1]))
for(f in 2:length(files)) {
  size_data = rbind(size_data, readRDS(paste0(dir, files[f])))
}

# size = size_data %>%
#   select(pval, method, functions, range) %>%
#   group_by(method, functions, range) %>%
#   summarize(size = mean(pval < 0.05)) %>%
#   ungroup() %>%
#   mutate(Statistic = recode(method, 
#                          K_WA = "K",
#                          Q_XD = "Q"),
#          Functions = as.factor(functions),
#          Range = as.factor(range)
#   )
# 
# ggplot(size, aes(x=Statistic, y=size, color=Functions)) +
#   geom_boxplot() +
#   geom_hline(yintercept = 0.05) +
#   theme_classic() +
#   theme(plot.title = element_text(hjust = 0.5)) +
#   xlab("Number of Functions (N)") +
#   ylab("Size")
#   # ggtitle("Size v.s. Number of Functions")
# ggsave(paste0("../assimilation-cfr/paper/size/", "size_functions.png"), width = 5, height = 3.2)
# 
# ggplot(size, aes(x=method, y=size, color=Range)) +
#   geom_boxplot() +
#   geom_hline(yintercept = 0.05) +
#   theme_classic() +
#   theme(plot.title = element_text(hjust = 0.5)) +
#   xlab("Method") +
#   ylab("Size") +
#   ggtitle("Size v.s. Range")
# ggsave(paste0("../assimilation-cfr/paper/size/", "size_range.png"), width = 5, height = 3.2)
# 
# 
# # boxplots over the seeds by range and function
# size = size_data %>%
#   select(pval, method, functions, range, seed) %>%
#   group_by(method, functions, range, seed) %>%
#   summarize(size = mean(pval < 0.05)) %>%
#   ungroup() %>%
#   mutate(Method = recode(method, 
#                          K_WA = "K Statistic",
#                          Q_XD = "Quality Index"),
#          Functions = as.factor(functions),
#          Range = as.factor(range)
#   )
# 
# ggplot(size, aes(x=Range, y=size, color=Method)) +
#   geom_boxplot() +
#   geom_hline(yintercept = 0.05) +
#   theme_classic() +
#   theme(plot.title = element_text(hjust = 0.5)) +
#   theme(strip.placement = "outside") +
#   facet_wrap(~Functions, strip.position = "bottom") +
#   xlab("Range by Number of Functions") +
#   ylab("Size") +
#   ggtitle("Size")


# boxplots over the seeds by function
size = size_data %>%
  select(pval, method, functions, seed) %>%
  group_by(method, functions, seed) %>%
  summarize(size = mean(pval < 0.05)) %>%
  ungroup() %>%
  mutate(Statistic = recode(method, 
                         K = "K",
                         Q = "Q"),
         Functions = as.factor(functions)
  )

ggplot(size, aes(x=Functions, y=size, color=Statistic)) +
  geom_boxplot() +
  geom_hline(yintercept = 0.05) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5)) +
  xlab("Number of functions (N)") +
  ylab("Size")
  # ggtitle("Size of K(F, G) v.s. Q(F, G)")
ggsave("research/assimilation-cfr/paper/size/size.png", width = 5, height = 3.2)



