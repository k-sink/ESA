# ESA manuscript graphics 
# Feb 2026
# Katharine Sink

library(tidyverse)
library(sf)
library(terra)
library(rnaturalearth)
library(rnaturalearthdata)
library(geodata)
library(ggnewscale)
library(patchwork)
library(maps)
library(data.table)

#################################
# LOCATION MAP
#################################
# get group numbers for each ESA clustering
ESA_groups <-  read.table("F:/ESA/ESApargrps.txt", sep = "")
ESA_groups <-  ESA_groups %>% rename(SITENO = V1, Partitions = V2, Slope = V3, ESA = V4, dESA = V5)
ESA_groups$SITENO <- str_pad(as.character(ESA_groups$SITENO), width = 8, side = "left", pad = "0")

# get site info, latitude and longitude
site_info <-  read_csv("F:/ESA/site_info.csv") %>% 
  dplyr::select(SITENO, dec_lat_va, dec_long_va)

ESA_groups <- ESA_groups %>% left_join(site_info, by = "SITENO")

# sites as sf 
sites_sf <- st_as_sf(ESA_groups, coords = c("dec_long_va", "dec_lat_va"), crs = 4326)

# define extents 
lon_min <- -105; lon_max <- -50
lat_min <- 20; lat_max <- 60
bbox_polygon <- st_as_sfc(st_bbox(c(xmin = lon_min, xmax = lon_max, 
                          ymin = lat_min, ymax = lat_max), crs = 4326))

# get polygons/vector data for states and rivers, crop 
sf::sf_use_s2(FALSE)

states_crop <- ne_states(country = "United States of America", returnclass = "sf") %>% 
  filter(!name %in% c("Alaska", "Hawaii", "Puerto Rico")) %>% st_crop(bbox_polygon)

# reproject to remove space
crs_albers <- 5070
states_proj <- st_transform(states_crop, crs_albers)
sites_proj  <- st_transform(sites_sf, crs_albers)

coord_limits <- coord_sf(expand = FALSE)
location_map <- ggplot() +
  geom_sf(data = states_proj, fill = NA, color = "black", linewidth = 0.4) +
  geom_sf(data = sites_proj, shape = 16, color = "black",
          size = 1.5, stroke = 0.4, alpha = 0.9) +
    labs(title = "Watershed Locations") +
  coord_limits +
    theme(axis.title = element_blank(),
    panel.grid.major = element_line(color = "blue", linewidth = 0.3),
    panel.grid.minor = element_line(color = "blue", linewidth = 0.2)) +
  theme_bw()

###########################################################
# LOCATION and ELEVATION MAP #
###########################################################
# get polygons/vector data for states and rivers, crop 
sf::sf_use_s2(FALSE)

states_crop <- ne_states(country = "United States of America", returnclass = "sf") %>% 
  filter(!name %in% c("Alaska", "Hawaii", "Puerto Rico")) %>% st_crop(bbox_polygon)
states_union <- st_union(states_crop)

land <- ne_download(scale = "medium", type = "land", category = "physical", 
                    returnclass = "sf") %>% st_crop(bbox_polygon)
land_union <- st_union(land)

rivers_union <- ne_download(scale = "medium", type = "rivers_lake_centerlines", 
                category = "physical", returnclass = "sf") %>% 
                st_crop(bbox_polygon) %>% st_intersection(states_crop)

lakes <- ne_download(scale = "medium", type = "lakes", category = "physical", 
                     returnclass = "sf") %>% st_crop(bbox_polygon)

# get raster dem data 
dem_us <- geodata::elevation_30s(country = "USA", path = tempdir())
dem_can <- geodata::elevation_30s(country = "CAN",  path = tempdir())
dem_mex <- geodata::elevation_30s(country = "MEX",  path = tempdir())

dem_all <- terra::merge(dem_us, dem_can, dem_mex)

dem_crop <- crop(dem_all, ext(lon_min, lon_max, lat_min, lat_max))
dem_crop <- aggregate(dem_crop, fact = 2)

# mask dem to land 
dem_masked <- mask(dem_crop, vect(states_union))

# terrain from masked DEM
#slope <- terrain(dem_masked, "slope", unit = "radians")
#aspect <- terrain(dem_masked, "aspect", unit = "radians")
#hillshade <- shade(slope, aspect, angle = 40, direction = 315)

slope <- terrain(dem_crop, "slope", unit = "radians")
aspect <- terrain(dem_crop, "aspect", unit = "radians")
hillshade <- shade(slope, aspect, angle = 40, direction = 315)

# convert dem to dataframe for plotting
#dem_df <- as.data.frame(dem_masked, xy = TRUE, na.rm = TRUE)
dem_df <- as.data.frame(dem_crop, xy = TRUE, na.rm = TRUE)
colnames(dem_df) <- c("x","y","elev")

hs_df <- as.data.frame(hillshade, xy = TRUE, na.rm = TRUE)
colnames(hs_df) <- c("x","y","shade")

# reproject to albers equal area 
# aea_crs <- "+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +datum=WGS84 +units=m +no_defs"

location_map <- ggplot() + 
  geom_raster(data = hs_df, aes(x = x, y = y, fill = shade), alpha = 0.35) +
  scale_fill_gradient(low = "black", high = "white", guide = "none") +
  ggnewscale::new_scale_fill() +
  geom_raster(data = dem_df, aes(x = x, y = y, fill = elev), alpha = 0.8) +
  scale_fill_gradientn(name = "Elevation (m)", 
  # colors = c("#2c5f2d","#97bc62","#c2b280","#8b5a2b"),
    colors = c("darkgreen","palegreen","tan","saddlebrown"),
    limits = c(5, 2050), oob = scales::squish, na.value = NA) + 
  geom_sf(data = states_crop, fill = NA, color = "black", linewidth = 0.4) +
  geom_sf(data = rivers_union, color = "blue", linewidth = 0.25) +
  geom_sf(data = sites_sf, shape = 21, fill = "yellow", color = "black", size = 2, stroke = 0.5) +
  geom_sf(data = lakes, color = "blue", linewidth = 0.25) + 
  scale_x_continuous(breaks = seq(-100,-70, 10)) +
  scale_y_continuous(breaks = seq(25,50,5)) +
  coord_sf(xlim = c(lon_min, lon_max),
           ylim = c(lat_min, lat_max),
           expand = FALSE, datum = st_crs(4326)) +
  theme_minimal() +
  theme(axis.title = element_blank(),
   panel.grid.major = element_blank(), 
   panel.grid.minor = element_blank()) +
   # panel.grid.major = element_line(color = "gray75", linewidth = 0.3),
   # panel.grid.minor = element_line(color = "gray85", linewidth = 0.2)) +
  labs(title = "Watershed Locations")

ggsave("F:/ESA/location_map.png", plot = location_map, dpi = 320, width = 6, height = 4, units = "in")

###########################################################
# ESA partition maps #
###########################################################
# sites as sf 
sites_sf$Partitions <- factor(sites_sf$Partitions, levels = 1:10)
cluster_colors <- c("#000000", "#999999", "#F0E442", "#56B4E9", "#009E73", 
                    "#CC79A7", "#0072B2", "#7F7F7F", "#D55E00", "#E69F00" )

cluster_shapes <- c(21, 22, 23, 24, 25, 21, 22, 23, 24, 25)

# location map with clusters (747 sites)
partition <- ggplot() + geom_sf(data = states_crop, fill = NA, color = "black", linewidth = 0.4) +
  geom_sf(data = sites_sf, aes(fill = factor(Partitions), shape = factor(Partitions)), 
          color = "black", size = 2, stroke = 0.4, alpha = 0.9) +
      scale_fill_manual(values = cluster_colors, name = "Partition") +
      scale_shape_manual(values = cluster_shapes, name = "Partition") +
  scale_x_continuous(breaks = seq(-100,-70, 10)) +
  scale_y_continuous(breaks = seq(25,50,5)) +
#  coord_sf(xlim = c(lon_min, lon_max),
#           ylim = c(lat_min, lat_max),
#           expand = FALSE, datum = st_crs(4326)) +
  theme_minimal() +
  theme(axis.title = element_blank(),  
   panel.grid.major = element_blank(), 
   panel.grid.minor = element_blank()) + 
   labs(title = "ESA Partitions")

###############################################
# groupings for distribution (RMSE), Group E (MAD/0.6745)
levels_5 <- c("1", "2", "3", "4", "5")

sites_sf$ESA <- factor(sites_sf$ESA, levels = levels_5) 

shared_colors <- scale_fill_manual(
  values = c("#D55E00", "#F0E442", "#009E73", "#0072B2", "#CC79A7"), 
  breaks = levels_5, name = "Group")

distribution <- ggplot() + geom_sf(data = states_crop, fill = NA, color = "black", linewidth = 0.4) +
  geom_sf(data = sites_sf, aes(fill = ESA), 
          shape = 21, color = "black", size = 1.5, stroke = 0.4, alpha = 0.9) +
     shared_colors +
    coord_sf(xlim = c(lon_min, lon_max),
           ylim = c(lat_min, lat_max),
           expand = FALSE, datum = st_crs(4326)) +
  labs(title = expression(ESA[rstd])) + 
  theme_minimal() + theme(axis.title = element_blank(),  
                          legend.title = element_blank(), 
                          legend.position = "none", 
                          axis.text = element_blank(), 
                          axis.ticks = element_blank(), 
                          panel.grid = element_blank(), 
                          plot.title = element_text(size = 12))

###############################################
# groupings for outliers dESA/dt, Group T
sites_sf$dESA <- factor(sites_sf$dESA, levels = levels_5)  

tail <- ggplot() + geom_sf(data = states_crop, fill = NA, color = "black", linewidth = 0.4) +
  geom_sf(data = sites_sf, aes(fill = dESA), 
          shape = 21, color = "black", size = 1.5, stroke = 0.4, alpha = 0.9) +
      shared_colors + 
  coord_sf(xlim = c(lon_min, lon_max),
           ylim = c(lat_min, lat_max), expand = FALSE, datum = st_crs(4326)) +
   labs(title = expression(ESA[tails])) + 
  theme_minimal() + theme(axis.title = element_blank(),  
                           legend.title = element_blank(), 
                          legend.position = "none", 
                          axis.text = element_blank(), 
                          axis.ticks = element_blank(), 
                          panel.grid.major = element_blank(), 
                          panel.grid.minor = element_blank(),  
                          plot.title = element_text(size = 12))

###############################################
# paneled figure for groupings
# reproject to remove space
crs_albers <- 5070
states_proj <- st_transform(states_crop, crs_albers)
sites_proj  <- st_transform(sites_sf, crs_albers)

# theme
clean_theme <- theme_minimal() +
  theme(axis.title = element_blank(), 
        axis.text = element_blank(), 
        axis.ticks = element_blank(), 
        panel.grid = element_blank(),
        legend.position = "bottom", 
        plot.margin = margin(0,0,0,0), 
        plot.title = element_text(size = 12, hjust = 0.5))

coord_limits <- coord_sf(expand = FALSE)
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
    plot.tag.position = c(0,1), 
    plot.tag = element_text(size = 12, hjust = 0, vjust = 0))  

ggsave("F:/ESA/ESA_groups_2.png", plot = final_map, dpi = 320, width = 10, height = 8, bg = "white")

#########################################################
### ESA TIME SERIES ###
#########################################################
# get time series output file from ESA calculation 
# 88 values for each basin (747 basins)
# modified slightly in excel (qep_output_march.csv)
esa_data <- read_csv("F:/ESA_R/esa_manuscript.csv") 
esa_data$DATE <- as.Date(esa_data$DATE, format = "%m/%d/%Y")        

# rename columns
esa_data <-  esa_data %>% rename(P = PRCP, Ep = PET, E = AET) 
  
# limit to 1 (overestimation of E, greater than P)
esa_data$ESA <- ifelse(esa_data$ESA > 1.0, 1.0, esa_data$ESA)

# set at datatable type 
esa_data <- setDT(esa_data)

######################################################
# paneled time series figure
# get unique basins per partition
n_df <- esa_data %>% 
  group_by(Partitions) %>% 
  summarise(n = n_distinct(SITENO)) %>% 
  arrange(Partitions)

# get partition levels
partitions <- n_df$Partitions

# create labels (a), (b), etc
facet_labels <- setNames(paste0("(", letters[seq_along(partitions)], ") Partition ", partitions, 
                                " (n = ", n_df$n, ")"), partitions)
# time series plot 
p <- ggplot() + 
  geom_line(data = esa_data, aes(x = DATE, y = ESA, group = SITENO), color = "darkgray", alpha = 0.25, linewidth = 0.3) +
#  geom_line(data = median_ESA, aes(x = DATE, y = median), color = "blue", linewidth = 0.8) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.8, color = "red") + 
  scale_y_continuous(limits = c(-1, 1), expand = c(0,0)) + 
  scale_x_date(breaks = seq(as.Date("1980-01-01"), as.Date("2025-01-01"), by = "5 years"),
    labels = scales::date_format("%Y"),
    limits = as.Date(c("1980-01-01", "2025-01-01"))) + 
  labs(y = "ESA", x = NULL) + 
  facet_wrap(~Partitions, nrow = 5, labeller = labeller(Partitions = facet_labels)) +
  theme_bw()+ 
  theme(legend.position = "none", 
        axis.title = element_text(size = 10), 
       # plot.title = element_text(size = 12, hjust = 0.5), 
        axis.text = element_text(size = 10), 
       strip.background = element_blank(), 
       strip.placement = "outside", 
       strip.text = element_text(face = "bold", hjust = 0))

###################################################################
# individual partition plots, filter by partition number
sites_sf_1 <- sites_sf %>% filter(Partitions == 1)
sites_sf_2 <- sites_sf %>% filter(Partitions == 2)
sites_sf_3 <- sites_sf %>% filter(Partitions == 3)
sites_sf_4 <- sites_sf %>% filter(Partitions == 4)
sites_sf_5 <- sites_sf %>% filter(Partitions == 5)
sites_sf_6 <- sites_sf %>% filter(Partitions == 6)
sites_sf_7 <- sites_sf %>% filter(Partitions == 7)
sites_sf_8 <- sites_sf %>% filter(Partitions == 8)
sites_sf_9 <- sites_sf %>% filter(Partitions == 9)
sites_sf_10 <- sites_sf %>% filter(Partitions == 10)

# change based on grouped partitions
sites_sf_2_3 <- sites_sf %>% filter(Partitions %in% c(2,3))

# location plot for each partition to inset into time series plot
# update relative to partition (title, cluster, etc)
map_partition_plot_1 <- ggplot() + geom_sf(data = states_crop, fill = NA, color = "black", linewidth = 0.4) +
  geom_sf(data = sites_sf_1, aes(fill = Partitions, shape = Partitions), 
          color = "black", size = 1.5, stroke = 0.4, alpha = 0.9) +
          scale_fill_manual(values = setNames(cluster_colors, 1:10)) +
          scale_shape_manual(values = setNames(cluster_shapes, 1:10)) +
          theme_minimal() +
          theme(axis.title = element_blank(),  
          axis.text = element_blank(), axis.ticks = element_blank(), 
          legend.position = "none", 
          panel.grid.major = element_blank(), panel.grid.minor = element_blank())

# similar groups (RMSE)
partition_1 <- esa_data %>% filter(Partitions == 1)
partition_2 <- esa_data %>% filter(Partitions == 2)
partition_3 <- esa_data %>% filter(Partitions == 3)
partition_4 <- esa_data %>% filter(Partitions == 4)
partition_5 <- esa_data %>% filter(Partitions == 5)
partition_6 <- esa_data %>% filter(Partitions == 6)
partition_7 <- esa_data %>% filter(Partitions == 7)
partition_8 <- esa_data %>% filter(Partitions == 8)
partition_9 <- esa_data %>% filter(Partitions == 9)
partition_10 <- esa_data %>% filter(Partitions == 10)

# median ESA values for plot by partition
median_ESA <- esa_data %>% group_by(Partitions, DATE) %>% summarise(median = median(ESA))
median_1 <- median_ESA %>% filter(Partitions == 1)
median_2 <- median_ESA %>% filter(Partitions == 2)
median_3 <- median_ESA %>% filter(Partitions == 3)
median_4 <- median_ESA %>% filter(Partitions == 4)
median_5 <- median_ESA %>% filter(Partitions == 5)
median_6 <- median_ESA %>% filter(Partitions == 6)
median_7 <- median_ESA %>% filter(Partitions == 7)
median_8 <- median_ESA %>% filter(Partitions == 8)
median_9 <- median_ESA %>% filter(Partitions == 9)
median_10 <- median_ESA %>% filter(Partitions == 10)

# time series plot with all time series (gray) and median (blue) by partition
# update with partition number, sites number, etc
plot_9 <-  ggplot() +
  geom_line(data = partition_9, aes(x = DATE, y = ESA, group = SITENO), 
            color = "darkgray", alpha = 0.25, linewidth = 0.3) +
  geom_line(data = median_9, aes(x = DATE, y = median), color = "blue", linewidth = 1) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.8, color = "red") +
  scale_y_continuous(limits = c(-1, 1), expand = c(0,0)) + 
  scale_x_date(breaks = seq(as.Date("1980-01-01"), as.Date("2025-01-01"), by = "5 years"),
    labels = scales::date_format("%Y"),
    limits = as.Date(c("1980-01-01", "2025-01-01"))) + 
  labs(title = expression(ESA[partition] ~ "9" ~ "(n = 106)"), y = "ESA", x = NULL) + theme_bw() + 
  theme(legend.position = "none", axis.title = element_text(size = 10), 
        plot.title = element_text(size = 10, hjust = 0.5), axis.text = element_text(size = 10)) 

# ggplot(data = partition_1, aes(x = ESA)) + geom_histogram(aes(y = after_stat(density)), bins = 44)
# histogram plots for distribution 
pdf_10 = ggplot(data = partition_10, aes(x = ESA)) + geom_histogram(aes(y = after_stat(count/sum(count))), bins = 40, 
       # binwidth = \(x) 2 * IQR(x) / length (x)^(1/3), 
         fill = "steelblue", color = "black") + 
  labs(x = "ESA", y = "PDF", 
       title = expression(atop("Partition 10",
                               ESA[rstd]~group~4))) +
         coord_cartesian(xlim = c(-1, 1)) + scale_y_continuous(expand = c(0, 0)) + theme_bw() +
  theme(plot.title = element_text(size = 10, hjust = 0.5), axis.title = element_text(size = 10), 
        axis.text = element_text(size = 10))

# aspect ratios 
map_partition_plot_1 <- map_partition_plot_1 + theme(aspect.ratio = 0.7)
plot_1 <- plot_1 + theme(aspect.ratio = 0.4)
pdf_1 <- pdf_1 + theme(aspect.ratio = 0.5)

# standardize plot margins
base_theme <- theme(plot.margin = margin(5, 5, 5, 5))

plot_1 <- plot_1 + base_theme
pdf_1 <- pdf_1 + base_theme
map_partition_plot_1 <- map_partition_plot_1 + base_theme

# one time series, one pdf
left_column <- plot_1 / plot_spacer() + plot_layout(heights = c(2,2))

right_column <- pdf_1 / plot_spacer() / map_partition_plot_1 + plot_layout(heights = c(2, 2, 3))

# two time series, two pdf
left_column <- plot_8 / plot_9 + plot_layout(heights = c(2,2))

right_column <- pdf_8 / pdf_9 / map_partition_plot_1 + plot_layout(heights = c(2, 2, 3))

# overall layout of 2 columns and panel labeling 
combined <- wrap_plots(left_column, right_column, ncol = 2, widths = c(2, 1)) +
 plot_annotation(tag_levels = 'a', tag_prefix = "", tag_suffix = ") ")
 
row_heights <- c(2,2)

# one partition
left_column <- plot_1 / pdf_1 + plot_layout(heights = row_heights)

map_col <- map_partition_plot_1

combined <- wrap_plots(left_column, map_col, ncol = 2, widths = c(2, 1)) + 
  plot_annotation(tag_levels = "a", tag_suffix = ")")


plot <- plot_1 + plot_2 + plot_3 + plot_4 + plot_5 + plot_6 + plot_7 + plot_8 + plot_9 + plot_10 +
  plot_layout(ncol = 2, axes = "collect") +  plot_annotation(tag_levels = "a", tag_prefix = "(", tag_suffix = ") ") &
  theme(plot.tag = element_text(size = 10))

ggsave(filename = "F:/ESA/time_series.png", 
       plot = plot, device = "png", width = 10, height = 8, dpi = 300)

pdf_combined <- pdf_2 + pdf_3 + pdf_4 + pdf_5 + pdf_6 + pdf_7 + pdf_8 + pdf_9 + pdf_10 + pdf_1 +
  plot_layout(ncol = 4) +  plot_annotation(tag_levels = "a", tag_prefix = "(", tag_suffix = ") ") &
  theme(plot.tag = element_text(size = 10))

ggsave(filename = "F:/ESA/pdf_plots.png", 
       plot = pdf_combined, device = "png", width = 10, height = 8, dpi = 300)

######################################################
### DERIVATIVE PLOTS ###
######################################################
# create 10-panel plot dESAda
esa_deriv_data_pdf <- left_join(esa_deriv_data, table_stats, by = "Partitions")

dESA_pdf <- function(partition_num, data, fence_table) {
  df <- data %>% filter(Partitions == partition_num)
  fences <- fence_table %>% filter(Partitions == partition_num)
  group_label <- df$dESA_group
  
  ggplot(data = df, aes(x = dESA)) +
    geom_histogram(
      aes(y = after_stat(count / sum(count))),
      bins  = 40,
      fill  = "steelblue",
      color = "black") +
    geom_vline(
      xintercept = fences$fenlo,
      color      = "firebrick",
      linetype   = "dashed",
      linewidth  = 0.6) +
    geom_vline(
      xintercept = fences$fenhi,
      color      = "firebrick",
      linetype   = "dashed",
      linewidth  = 0.6) +
    labs(x = "dESA/dt", y = "PDF",
      title = bquote(atop(
        "Partition" ~ .(partition_num),
        ESA[tails] ~ group ~ .(group_label)))) +
    coord_cartesian(xlim = c(-3, 3)) +
    scale_y_continuous(expand = c(0, 0)) +
    theme_bw() +
    theme(
      plot.title = element_text(size = 10, hjust = 0.5),
      axis.title = element_text(size = 10),
      axis.text  = element_text(size = 10))
}

# build a fence lookup table with one row per partition
# using your fndout_R function applied per partition
fence_table <- esa_deriv_data_pdf %>%
  group_by(Partitions) %>%
  group_modify(~{
    out <- fndout_R(.x$dESA)
    tibble(fenlo = out$fenlo, fenhi = out$fenhi)
  }) %>%
  ungroup()

# build all 10 panels in one pass
dESA_pdf_list <- map(1:10, ~ dESA_pdf(.x, esa_deriv_data_pdf, fence_table))

# combine into a 4-column panel, ordered by ESA_tails_group like your ESA PDF figure
panel_order <- c(1, 10, 9, 3, 4, 8, 2, 6, 7, 5)

dESA_pdf_combined <- wrap_plots(dESA_pdf_list[panel_order], ncol = 4) +
  plot_annotation(
    tag_levels = "a", tag_prefix = "(", tag_suffix = ") "
  ) &
  theme(plot.tag = element_text(size = 10))

ggsave(filename = "F:/ESA/pdf_dESA_plots.png", 
       plot = dESA_pdf_combined, device = "png", width = 10, height = 8, dpi = 300)

######################################################
### HEATMAP FOR CLUSTERING ###
######################################################

group_table <- tibble(
  Partition = c(4, 3, 5, 2, 7, 6, 8, 9, 10, 1), 
  ESArstd_grp = c(1, 1, 1, 1, 2, 2, 3, 3, 4, 5), 
  ESAtails_grp = c(3, 3, 5, 4, 4, 4, 3, 2, 1, 1), 
  Median_ESA = c(0.24, 0.32, 0.16, 0.35, -0.20, -0.14, -0.25, -0.38, -0.64, 0.90)) %>% 
  arrange(Partition) %>% 
  mutate(Partition = factor(Partition, levels = 1:10))

group_long <- group_table %>% 
  pivot_longer(
    cols = c(ESArstd_grp, ESAtails_grp), 
    names_to = "Classification", 
    values_to = "Group") %>% 
  mutate(
    Classification = recode(Classification, 
                            "ESArstd_grp" = "ESA[rstd]~Group", 
                            "ESAtails_grp" = "ESA[tails]~Group"), 
    Group = factor(Group, levels = 1:5))

 group_colors <- c(
   "1" = "#D55E00", 
   "2" = "#F0E442", 
   "3" = "#009E73", 
   "4" = "#0072B2", 
   "5" = "#CC79A7") 

grps_plot <- ggplot(group_long,
       aes(x = Classification, y = Partition, fill = Group)) +
  geom_tile(color = "white", linewidth = 1.2) +
  geom_text(aes(label = as.character(Group)),
            color = "white", fontface = "bold", size = 5) +
  scale_fill_manual(values = group_colors,
                    name   = "Group") +
  scale_x_discrete(
    labels = c(
      "ESA[rstd]~Group"  = expression(ESA[rstd]),
      "ESA[tails]~Group" = expression(ESA[tails])
    )
  ) +
  scale_y_discrete(limits = rev) +
  labs(
 #   title = "Second-Stage Classification Summary",
    x     = NULL,
    y     = "ESA Partition\n(1 = most arid, 10 = most humid)"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid    = element_blank(),
    axis.text.x   = element_text(size = 11, face = "bold"),
    axis.text.y   = element_text(size = 10),
    plot.title    = element_text(size = 12, face = "bold",
                                 hjust = 0.5),
    legend.position = "right"
  ) 
 

ggsave(filename = "F:/ESA/grps_plot.png", 
       plot = grps_plot, device = "png", width = 10, height = 8, dpi = 300)
######################################################
### INDICES PLOTS ###
######################################################

# grouped bar chart with r values by partition
cor_long <- cor_results_final %>%
  select(Partitions, r_SPI_12, r_SPEI_12, r_SPEI_24) %>%
  pivot_longer(
    cols      = starts_with("r_"),
    names_to  = "Index",
    values_to = "r") %>%
  mutate(
    Index = recode(Index,
      "r_SPI_12"  = "SPI-12",
      "r_SPEI_12" = "SPEI-12",
      "r_SPEI_24" = "SPEI-24"),
    Partitions = factor(Partitions)
  )

indices_cor_plot <- ggplot(cor_long,
       aes(x = Partitions, y = r, fill = Index)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6) +
  geom_hline(yintercept = 0, color = "black", linewidth = 0.5) +
  scale_fill_manual(
    values = c("SPI-12"  = "#2166ac",
               "SPEI-12" = "#d6604d",
               "SPEI-24" = "#92c5de")) +
  scale_y_continuous(
    limits = c(-0.8, 0),
    breaks = seq(-0.8, 0, 0.2),
    labels = seq(-0.8, 0, 0.2)) +
  labs(
  #  title = "Pearson Correlation between ESA and Drought Indices by Partition",
    x     = "ESA Partition",
    y     = "Pearson r",
    fill  = "Index"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    legend.position  = "bottom",
    panel.grid.minor = element_blank(),
    plot.title       = element_text(size = 10, face = "bold")
  )

ggsave(filename = "F:/ESA/indices_cor_plot.png", 
       plot = indices_cor_plot, device = "png", width = 10, height = 8, dpi = 300)

####
# compute partition median ESA and mean SPEI-12 at each time point
# for representative partitions
rep_partitions <- c(1, 7, 10)

ts_plot_data <- esa_spi12_spei12 %>%
  filter(Partitions %in% rep_partitions) %>%
  group_by(Partitions, DATE) %>%
  summarise(
    median_ESA    = median(ESA,     na.rm = TRUE),
    mean_SPEI_12  = mean(SPEI_12,   na.rm = TRUE),
    mean_SPI_12   = mean(SPI_12,    na.rm = TRUE),
    .groups       = "drop") %>%
  mutate(Partition_label = paste("Partition", Partitions))

# scale SPEI-12 to ESA range for dual axis
# ESA range approximately -1 to 1
# SPEI range approximately -3 to 3
# scale factor: 1/3
scale_factor <- 1/3

# function to build one time series panel for a given partition
make_ts_panel <- function(partition_num, data, scale_factor,
                          show_legend = TRUE) {

  df <- data %>% filter(Partitions == partition_num)

  p <- ggplot(df, aes(x = DATE)) +
    geom_line(aes(y = median_ESA, color = "ESA"),
              linewidth = 0.6) +
    geom_line(aes(y = mean_SPEI_12 * scale_factor, color = "SPEI-12"),
              linewidth = 0.5) +
    geom_hline(yintercept = 0,
               linetype   = "dotted",
               color      = "red",
               linewidth  = 0.8) +
    annotate("rect",
             xmin = as.Date("1988-01-01"),
             xmax = as.Date("1990-12-31"),
             ymin = -Inf, ymax = Inf,
             alpha = 0.08, fill = "orange") +
    annotate("rect",
             xmin = as.Date("2012-01-01"),
             xmax = as.Date("2012-12-31"),
             ymin = -Inf, ymax = Inf,
             alpha = 0.08, fill = "orange") +
    scale_y_continuous(
      name     = "Median ESA",
      limits   = c(-1.2, 1.2),
      breaks   = seq(-1, 1, 0.5),
      sec.axis = sec_axis(
        transform = ~ . / scale_factor,
        name      = "Mean SPEI-12",
        breaks    = seq(-3, 3, 1)
      )
    ) +
    scale_color_manual(
      values = c("ESA"     = "#593196",
                 "SPEI-12" = "#d6604d")
    ) +
    scale_x_date(
      date_breaks = "5 years",
      date_labels = "%Y"
    ) +
    labs(
      title = paste("Partition", partition_num),
      x     = NULL,
      color = NULL
    ) +
    theme_minimal() +
    theme(
      legend.position  = if (show_legend) "bottom" else "none",
      axis.title       = element_text(size = 10),
      axis.text        = element_text(size = 10),
      panel.grid.minor = element_blank(),
      plot.title       = element_text(size = 10, hjust = 0.5)
    )

  p
}

# build all 3 panels
ts_panel_list <- map(rep_partitions,
                      ~ make_ts_panel(.x, ts_plot_data, scale_factor))

# combine with patchwork, shared legend, tagged (a), (b), (c)
ts_combined <- wrap_plots(ts_panel_list, ncol = 1) +
  plot_layout(guides = "collect") +
  plot_annotation(
    tag_levels = "a", tag_prefix = "(", tag_suffix = ") "
  ) &
  theme(
    legend.position = "bottom",
    plot.tag        = element_text(size = 10)
  )

ggsave(filename = "F:/ESA/ts_combined.png", 
       plot = ts_combined, device = "png", width = 10, height = 8, dpi = 300)

