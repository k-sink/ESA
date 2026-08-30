##############################
### MAPS ###
##############################
# Katharine Sink

# maps for spatial distribution of watersheds, partitions, and groups
# additional maps created in QGIS

library(tidyverse)
library(sf)
library(terra)
library(rnaturalearth)
library(rnaturalearthdata)
library(maps)

###################
# LOAD DATA
###################
# add group numbers to watershed list
ESA_par_groups <- read_csv("kmeans_r_results.csv") %>% 
  dplyr::select(SITENO, Partition, dec_lat_va, dec_long_va) %>% 
  left_join(partition_groups %>% dplyr::select(Partition, ESArstd_group, ESAtails_group), by = "Partition")

# coordinates as sf object
sites_sf <- st_as_sf(ESA_par_groups, coords = c("dec_long_va", "dec_lat_va"), crs = 4326)

# add all factor columns 
sites_sf$Partition <- factor(sites_sf$Partition, levels = 1:10)

levels_5 <- c("1", "2", "3", "4", "5")
sites_sf$ESA  <- factor(sites_sf$ESArstd_group,  levels = levels_5)
sites_sf$dESA <- factor(sites_sf$ESAtails_group, levels = levels_5)

#############################
# DEFINE EXTENT, PROJECTION
#############################
# create bounding box for watershed extent
lon_min <- -105; lon_max <- -50
lat_min <- 20; lat_max <- 60
bbox_polygon <- st_as_sfc(st_bbox(c(xmin = lon_min, xmax = lon_max, 
                          ymin = lat_min, ymax = lat_max), crs = 4326))

# get polygons/vector data for states
sf::sf_use_s2(FALSE)

# crop states (CONUS)
states_crop <- ne_states(country = "United States of America", returnclass = "sf") %>% 
  filter(!name %in% c("Alaska", "Hawaii", "Puerto Rico")) %>% st_crop(bbox_polygon)


# reproject to Albers Equal Area Conic (EPSG:5070), the standard CONUS
# projection, fixes map compression: coord_sf() and enforces an accurate aspect ratio 
# compared to lon/lat (EPSG:4326) where degrees of longitude and latitude don't
# represent equal real-world distances
crs_albers <- 5070
states_proj <- st_transform(states_crop, crs_albers)
sites_proj  <- st_transform(sites_sf, crs_albers)


######################
# SHARED STYLING 
######################
# shared colors and theme
shared_colors <- scale_fill_manual(
  values = c("#D55E00", "#F0E442", "#009E73", "#0072B2", "#CC79A7"), 
  breaks = levels_5, name = "Group")

clean_theme <- theme_minimal() +
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank(),
        legend.position = "bottom",
        plot.margin = margin(0, 0, 0, 0),
        plot.title = element_text(size = 12, hjust = 0.5))

coord_limits <- coord_sf(expand = FALSE)

#########################################
# ESA PARTITION MAP (partitions 1 - 10)
#########################################
cluster_colors <- c("#000000", "#999999", "#F0E442", "#56B4E9", "#009E73", 
                    "#CC79A7", "#0072B2", "#7F7F7F", "#D55E00", "#E69F00" )

cluster_shapes <- c(21, 22, 23, 24, 25, 21, 22, 23, 24, 25)

# location map with clusters (747 sites)
partition_map <- ggplot() +
  geom_sf(data = states_proj, fill = NA, color = "black", linewidth = 0.4) +
  geom_sf(data = sites_proj, aes(fill = Partition, shape = Partition),
          color = "black", size = 2, stroke = 0.4, alpha = 0.9) +
  scale_fill_manual(values = cluster_colors, name = "Partition") +
  scale_shape_manual(values = cluster_shapes, name = "Partition") +
  scale_x_continuous(breaks = seq(-100, -70, 10)) +
  scale_y_continuous(breaks = seq(25, 50, 5)) +
  coord_sf(datum = st_crs(4326), expand = FALSE) +
  theme_minimal() +
  theme(axis.title = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  labs(title = "ESA Partitions")

ggsave(here("Figures/ESA_partition_map.png"), plot = partition_map, dpi = 320,
       width = 8, height = 6, bg = "white")

######################################
# ESA RSTD AND TAILS MAPS
######################################
distribution_map <- ggplot() +
  geom_sf(data = states_proj, fill = NA, color = "black", linewidth = 0.4) +
  geom_sf(data = sites_proj, aes(fill = ESA), shape = 21, color = "black",
          size = 2.25, stroke = 0.4, alpha = 0.9) +
  labs(title = expression(ESA[rstd])) +
  shared_colors + clean_theme + coord_limits

tail_map <- ggplot() +
  geom_sf(data = states_proj, fill = NA, color = "black", linewidth = 0.4) +
  geom_sf(data = sites_proj, aes(fill = dESA), shape = 21, color = "black",
          size = 2.25, stroke = 0.4, alpha = 0.9) +
  labs(title = expression(ESA[tails])) +
  shared_colors + clean_theme + coord_limits

final_map <- distribution_map + tail_map +
  plot_annotation(tag_levels = "a", tag_prefix = "(", tag_suffix = ")") +
  plot_layout(ncol = 2, guides = "collect") &
  theme(legend.position = "bottom",
        plot.margin = margin(0, 0, 0, 0),
        panel.spacing = unit(0, "pt"),
        plot.tag.position = c(0, 1),
        plot.tag = element_text(size = 12, hjust = 0, vjust = 0))

ggsave(here("Figures/ESA_groups_map.png"), plot = final_map, dpi = 320,
       width = 10, height = 8, bg = "white")