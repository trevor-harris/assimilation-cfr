rm(list = ls())
gc()




########### READ ME #############

# you must change the working directory to be the submitted_code folder
# none of this will work otherwise
# mine is left here as an example

########## Example
# setwd("/Users/trevh/research/assimilation-cfr/submitted_code/")

#################################





library(R.matlab)
library(ggplot2)


prox = readMat("maps/proxydata_aprmar_lmr_v0.2.0_pages2k_v2.mat")


### count the number of non-missing values each year to get the number of proxies each year
tmap = function(ind) {
  tmap = prox$lmr2k.data[851:1848,ind]
  tmap = 1-apply(tmap, 1, is.nan)
  colSums(tmap)
}

## stratify by archive (proxy) type
arch = lapply(prox$archive, function(x) {
  as.character(unlist(x))
})

# big 3
tree_rings = which(arch == "Tree Rings")
corals = which(arch == "Corals and Sclerosponges")
ice_cores = which(arch == "Ice Cores")
major = which(arch %in% c("Tree Rings", "Corals and Sclerosponges", "Ice Cores"))

# subset by archive type
tmap_tr = tmap(tree_rings)
tmap_cr = tmap(corals)
tmap_ic = tmap(ice_cores)
tmap_ot = tmap(-major)
tmap_gg = data.frame("tr" = tmap_tr, "cr" = tmap_cr, "ic" = tmap_ic, "ot" = tmap_ot, "Year" = 851:1848)


# plot with area plots
ggplot(tmap_gg) +
  geom_area(aes(Year, tr,
                fill = "Tree Rings"),
            color = 'black') +
  geom_area(aes(Year, ic,
                fill = "Ice Cores"),
            color = 'black') +
  geom_area(aes(Year, ot,
                fill = "Other"),
            color = 'black') +
  geom_area(aes(Year, cr,
                fill = "Corals"),
            color = 'black') +
  scale_fill_manual("", 
                    breaks = c("Tree Rings", "Ice Cores", "Corals", "Other"),
                    values = c("coral", "skyblue", "grey30", "forestgreen")) +
  theme_classic() +
  theme(legend.position="bottom",
        text = element_text(size = 16)) +
  ylab("Number of sites")
ggsave("maps/temporal_distribution.png", width = 6, height = 4)


# plot with area plots (log scale)
ggplot(tmap_gg) +
  geom_area(aes(Year, log1p(tr),
                fill = "Tree Rings")) +
  geom_area(aes(Year, log1p(ic),
                fill = "Ice Cores")) +
  geom_area(aes(Year, log1p(ot),
                fill = "Other")) +
  geom_area(aes(Year, log1p(cr),
                fill = "Corals")) +
  scale_fill_manual("", 
                    breaks = c("Tree Rings", "Ice Cores", "Corals", "Other"),
                    values = c("coral", "skyblue", "grey30", "forestgreen")) +
  theme_classic() +
  theme(legend.position="bottom",
        text = element_text(size = 16)) +
  # ylab("log counts") +
  scale_y_continuous("Number of sites", c(0, 2, 4, 6, 8),
                     labels = as.integer(exp(c(0, 2, 4, 6, 8))))
ggsave("maps/log_temporal_distribution.png", width = 6, height = 4)