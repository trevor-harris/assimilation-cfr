rm(list = ls())
gc()




########### READ ME #############

# you must change the working directory to be the submitted_code folder
# none of this will work otherwise
# mine is left here as an example

########## Example
# setwd("/Users/trevh/research/assimilation-cfr/submitted_code/")

#################################





library(ncdf4)
library(tictoc)
library(ggplot2)
library(dplyr)
library(reshape2)
library(R.matlab)
library(rgdal)
library(mapproj)
library(raster)
library(gridExtra)
library(cowplot)

devtools::install_github('trevor-harris/kstat')
library(kstat)


# helper function for multiplot
get_legend<-function(myggplot){
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

# import world map overlay
# this throws a warning. Ignore it.
data("wrld_simpl", package = "maptools")                                                                   
wm <- crop(wrld_simpl, extent(-180, 180, 45, 90))

# read in proxy data
prox = readMat("maps/proxydata_aprmar_lmr_v0.2.0_pages2k_v2.mat")

## stratify by archive (proxy) type
arch = lapply(prox$archive, function(x) {
  as.character(unlist(x))
})

### add proxy locations (before 1600)
smap = prox$lmr2k.data[1500,]
smap = 1-sapply(smap, is.nan)
smap = which(smap == 1)

# extract lat and lon and reparameterize lon to match ggplot
lat = prox$p.lat[smap]
lon = prox$p.lon[smap]
lon = ifelse(lon >= 180, lon - 360, lon)
prox_loc = data.frame(lat = lat, lon = lon)

# extract labels for each of the three major proxy types
arch = arch[smap]
tree_rings = which(arch == "Tree Rings")
corals = which(arch == "Corals and Sclerosponges")
ice_cores = which(arch == "Ice Cores")
major = which(arch %in% c("Tree Rings", "Corals and Sclerosponges", "Ice Cores"))

prox_loc[["archive"]] = "Other"
prox_loc[tree_rings, 'archive'] = "Tree Rings"
prox_loc[corals, 'archive'] = "Corals and \nSclerosponges"
prox_loc[ice_cores, 'archive'] = "Ice Cores"

# subset to the arctic and count number of tree rings
prox_loc = prox_loc[!is.na(prox_loc$lat),]
prox_loc = prox_loc[lat > 50,]
nrow(prox_loc[prox_loc$archive == "Tree Rings",])


# Defines the x axes required
x_lines <- seq(-120,180, by = 60)


# Pre 1600 plot (lp)
lp <- ggplot() +
  geom_polygon(data = wm, aes(x=long, y=lat, group=group), 
               fill = "grey90", colour = "black", alpha = 0.8) +
  geom_point(data = prox_loc, aes(x=lon, y=lat, fill=archive),
             color = "black", shape = 25, size = 4) +
  coord_map("ortho", orientation = c(90, 0, 0)) +
  scale_y_continuous(breaks = seq(45, 90, by = 10), labels = NULL) +
  scale_x_continuous(breaks = NULL) +
  xlab("") + 
  ylab("") +
  theme_void() +
  
  # Adds labels
  # geom_text(aes(x = 180, y = seq(55, 85, by = 10), hjust = -0.2, label = paste0(seq(55, 85, by = 10), "°N"))) +
  geom_text(aes(x = x_lines, y = 36, label = c("120°W", "60°W", "0°", "60°E", "120°E", "180°W"))) +
  
  # Adds axes
  geom_hline(aes(yintercept = 45), size = 0.8)  +
  geom_segment(aes(y = 45, yend = 90, x = x_lines, xend = x_lines), linetype = "dashed") +
  
  theme(panel.background = element_blank(),
        panel.grid.major = element_line(size = 0.25, linetype = 'dashed',
                                        colour = "black"),
        axis.ticks=element_blank()) +
  scale_fill_manual("Proxy Type", 
                    breaks = c("Tree Rings", "Ice Cores", "Other"),
                    values = c("#00BFC4", "grey30", "#7CAE00"))


### add proxy locations (after 1600)
smap = prox$lmr2k.data[1700,]
smap = 1-sapply(smap, is.nan)
smap = which(smap == 1)

## stratify by archive (proxy) type
arch = lapply(prox$archive, function(x) {
  as.character(unlist(x))
})

# extract lat and lon and reparameterize lon to match ggplot
lat = prox$p.lat[smap]
lon = prox$p.lon[smap]
lon = ifelse(lon >= 180, lon - 360, lon)
prox_loc = data.frame(lat = lat, lon = lon)

# extract labels for each of the three major proxy types
arch = arch[smap]
tree_rings = which(arch == "Tree Rings")
corals = which(arch == "Corals and Sclerosponges")
ice_cores = which(arch == "Ice Cores")
major = which(arch %in% c("Tree Rings", "Corals and Sclerosponges", "Ice Cores"))

prox_loc[["archive"]] = "Other"
prox_loc[tree_rings, 'archive'] = "Tree Rings"
prox_loc[corals, 'archive'] = "Corals and \nSclerosponges"
prox_loc[ice_cores, 'archive'] = "Ice Cores"

# subset to the arctic and count number of tree rings
prox_loc = prox_loc[!is.na(prox_loc$lat),]
prox_loc = prox_loc[lat > 50,]
nrow(prox_loc[prox_loc$archive == "Tree Rings",])


# Post 1600 plot (rp)
rp <- ggplot() +
  geom_polygon(data = wm, aes(x=long, y=lat, group=group), 
               fill = "grey90", colour = "black", alpha = 0.8) +
  geom_point(data = prox_loc, aes(x=lon, y=lat, fill=archive),
             color = "black", shape = 25, size = 4) +
  coord_map("ortho", orientation = c(90, 0, 0)) +
  scale_y_continuous(breaks = seq(45, 90, by = 45), labels = NULL) +
  scale_x_continuous(breaks = NULL) +
  xlab("") + 
  ylab("") +
  theme_void() +
  geom_text(aes(x = x_lines, y = 36, label = c("120°W", "60°W", "0°", "60°E", "120°E", "180°W"))) +
  
  # Adds axes
  geom_hline(aes(yintercept = 45), size = 0.8)  +
  geom_segment(aes(y = 45, yend = 90, x = x_lines, xend = x_lines), linetype = "dashed") +
  
  theme(panel.background = element_blank(),
        panel.grid.major = element_line(size = 0.25, linetype = 'dashed',
                                        colour = "black"),
        axis.ticks=element_blank()) +
  theme(text = element_text(size = 16)) +
  scale_fill_manual("Proxy Type",
                    breaks = c("Tree Rings", "Ice Cores", "Other"),
                    values = c("#00BFC4", "grey30", "#7CAE00")) + 
  theme(legend.position="none")

# Save the legend
legend <- get_legend(lp)
lp <- lp + theme(legend.position="none")

# Arrange plot
arctic_tree_rings = grid.arrange(lp, rp, legend, ncol=3, widths=c(2.4, 2.4, 0.8))

save_plot("maps/arctic_tree_rings.png", arctic_tree_rings, base_aspect_ratio = 2.9)
