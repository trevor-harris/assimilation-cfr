nlat.centers = 20
nlon.centers = 20
nc = nc.prior

n.center = nlat.centers * nlon.centers

lats = as.vector(nc$dim$lat$vals)
lons = as.vector(nc$dim$lon$vals)
n.lon = nc$dim$lon$len
n.lat = nc$dim$lat$len

# lats.cen = seq(min(lats), max(lats), length.out = nlat.centers)
lons.cen = seq(min(lons), max(lons), length.out = nlon.centers)
lons.cen = lons.cen -180

# cosine normalize to add more basis towards tropics and away from poles
lats.pi = cos((lats / max(lats)) * (pi/2))
lats.cen = seq(min(-lats.pi), max(lats.pi), length.out = nlat.centers)
lats.cen = (acos(lats.cen) / (pi/2) * max(lats)) - 90
lats.cen = -lats.cen

grid = expand.grid(lons.cen, lats.cen)
colnames(grid) = c("Longitude", "Latitude")
ggplot(grid) +
  borders("world", colour = "grey30") +
  geom_point(aes(x = Longitude, y = Latitude, color = "red"), show.legend = F) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        plot.title = element_text(hjust = 0.5)) +
  labs(x = "", y = "")
ggsave("paper/figures/basis_grid.png", width = 5, height = 3.2)

