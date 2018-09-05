rm(list = ls())

matsplitter<-function(M, r, c) {
  # splits 1 matrix into c MxR matricies
  # I have no idea how this works
  
  rg <- (row(M)-1)%/%r+1
  cg <- (col(M)-1)%/%c+1
  rci <- (rg-1)*max(cg) + cg
  N <- prod(dim(M))/r/c
  cv <- unlist(lapply(1:N, function(x) M[rci==x]))
  dim(cv)<-c(r,c,N)
  cv
}

regions = 64
sims = 100
batches = 25
run = 1

data.dir = paste0("/Users/trevh/research/assimilation-cfr/simdata/run", run, "/")
list.files(data.dir)
files = list.files(data.dir)

dep = array(0, dim=c(regions, sims, batches))
bonf = array(0, dim=c(regions, sims, batches))
pw = array(0, dim=c(regions, sims, batches))

d = 1
b = 1
p = 1
for(f in files) {
  if(grepl("Depth", f)) {
    dep[,,d] = readRDS(paste0(data.dir, f))
    d = d+1
  }
  if(grepl("Bonferroni", f)) {
    bonf[,,b] = readRDS(paste0(data.dir, f))
    b = b+1
  }
  if(grepl("Pointwise", f)) {
    pw[,,p] = readRDS(paste0(data.dir, f))
    p = p+1
  }
  if(grepl("post_mu", f)) {
    post_mu = readRDS(paste0(data.dir, f))
  }
}

# reformat 
d_power = matrix(0, batches, regions)
b_power = matrix(0, batches, regions)
p_power = matrix(0, batches, regions)
for(k in 1:batches) {
  d_power[k,] = apply(dep[,,k], 1, mean)
  b_power[k,] = apply(bonf[,,k], 1, mean)
  p_power[k,] = apply(pw[,,k], 1, mean)
}

# sort by L2
regions_mu = matsplitter(post_mu, 5, 5)
l2 = sapply(1:regions, function(x) sum(regions_mu[,,x]^2))

# ggplot2 version
library(reshape2)
library(ggplot2)

d_df = melt(d_power)
d_df[["Type"]] = "Depth"
d_df[["L2"]] = rep(l2, each=batches)
names(d_df) = c("Batch", "Region", "Power", "Type", "L2")

b_df = melt(b_power)
b_df[["Type"]] = "Bonferroni"
b_df[["L2"]] = rep(l2, each=batches)
names(b_df) = c("Batch", "Region", "Power", "Type", "L2")

p_df = melt(p_power)
p_df[["Type"]] = "Pointwise"
p_df[["L2"]] = rep(l2, each=batches)
names(p_df) = c("Batch", "Region", "Power", "Type", "L2")

powers = rbind(d_df, b_df, p_df)

mu_df = melt(post_mu)
ggplot(data = mu_df, aes(x = Var2, y = Var1)) +
  geom_raster(aes(fill = value), interpolate = T) +
  scale_fill_distiller(palette="Spectral") +
  labs(title = paste0("Posterior Mean - Test ", run)) +
  theme_void() +
  theme(plot.title = element_text(hjust = 0.5))

ggsave(paste0("plots/run", run, "mean", ".png"), width = 5, height = 3.2)

ggplot(data = powers, aes(x = L2, y = Power, color = Type)) +
  geom_point(alpha = 0.3) +
  geom_smooth(method = "loess") +
  ylim(0,1) +
  labs(title = paste0("Power Comparison - Test ", run)) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5))

ggsave(paste0("plots/run", run, "power", ".png"), width = 5, height = 3.2)
  

