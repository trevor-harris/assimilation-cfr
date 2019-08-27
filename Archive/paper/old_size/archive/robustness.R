rm(list = ls()); gc()
library(tictoc)

# get args
seed = 042696
pts = 20
sims = 500

source("research/assimilation-cfr/code/depth_tests.R")
source("research/assimilation-cfr/code/depths.R")
source("research/assimilation-cfr/code/simulation.R")


#### SIZE
set.seed(seed)

n = 25
trials = 1:20
npost = c(10, 25, 50, 100, 200, 300)

rows = length(trials)*length(npost)*sims

kvals = data.frame(pval = numeric(rows),
                   trial = numeric(rows),
                   n2 = numeric(rows))

qvals = data.frame(pval = numeric(rows),
                   trial = numeric(rows),
                   n2 = numeric(rows))

i=1
for(t in trials) {
  for(n2 in 1:length(npost)) {
    for(s in 1:sims) {
      tic("Total")
      cat("Simulation ", s, "\n")
      
      f = gp1d(fields = n, pts = pts, l = 30)
      g = gp1d(fields = npost[n2]+1, pts = pts, l = 30)
      
      kvals[i,] = c(kolm(f, g)[2], t, npost[n2])
      qvals[i,] = c(quality(f, g)[2], t, npost[n2])

      i = i + 1

      toc()
    }
  }
}

kvals[["method"]] = "K"
qvals[["method"]] = "Q"
vals = rbind(kvals, qvals)

size = vals %>%
  select(pval, method, n2, trial) %>%
  group_by(method, n2, trial) %>%
  summarize(size = mean(pval < 0.05)) %>%
  ungroup() %>%
  mutate(Method = recode(method, 
                         K = "K Statistic",
                         Q = "Quality Index"),
         Functions = as.factor(n2)
  )

ggplot(size, aes(x=Functions, y=size, color=Method)) +
  geom_boxplot() +
  geom_hline(yintercept = 0.05) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5)) +
  xlab("Functions in Y (X has 25)") +
  ylab("Size") +
  ggtitle("Size")

saveRDS(vals, file = "research/assimilation-cfr/paper/size/function_robustness.RDS")
