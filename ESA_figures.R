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
# Hodograph
#################################
basin <- read_csv("D:/Journal submissions/Scientific Data/MACH_dataset/MACH_ts/ALLVAR/basin_07340300_MACH.csv")
indices <- basin %>% select(DATE, AET, PET, PRCP)
indices <- indices %>% mutate(MNTH = month(DATE), YR = year(DATE))
semi_annual <- indices %>% mutate(period = ifelse(MNTH <=6, 1, 2)) %>% 
  group_by(YR, period) %>% 
  summarise(PET = sum(PET, na.rm = TRUE), 
            AET = sum(AET, na.rm = TRUE), 
            PRCP = sum(PRCP, na.rm = TRUE))
semi_annual <- semi_annual %>% mutate(EpP = PET/PRCP, EP = AET/PRCP)

annual <- indices %>% group_by(YR) %>% 
  summarise(PRCP = sum(PRCP, na.rm = TRUE), 
            AET = sum(AET, na.rm = TRUE), 
            PET = sum(PET, na.rm = TRUE), 
            EpP = PET/PRCP, EP = AET/PRCP, 
            .groups = "drop") %>% 
  arrange(YR) %>% mutate(decade = paste0(floor(YR / 10) * 10, "s"))

ggplot(data = annual, aes(x = EpP, y = EP)) + 
  geom_line(aes(color = decade), linewidth = 1) +
  geom_point(aes(color = decade), size = 2) + theme_bw()


ggplot(annual, aes(x = EpP, y = EP)) + 
  geom_path(aes(group = 1), linewidth = 1.2, color = "grey70") +
  geom_point(aes(color = decade), size = 3) + theme_bw()

ggplot(data = annual, aes(x = EpP, y = EP)) + 
  geom_path(aes(color = YR, group = 1), linewidth = 1.2) +
  geom_point(aes(color = YR), size = 2) +
  scale_color_viridis_c(option = "mako") + 
  coord_cartesian(xlim = c(0, 2), ylim = c(0, 1)) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.2)) +
  geom_hline(yintercept = 0) + geom_vline(xintercept = 0)

budyko_function = function(PET_P) {
  ifelse(PET_P == 0,0, 
         (PET_P * tanh(1/PET_P) * (1 - exp(-PET_P)))^0.5)
}

ggplot(data = annual, aes(x = EpP, y = EP)) + 
  stat_function(fun = budyko_function, color = "red", linewidth = 1, linetype = 1) +
  geom_segment(aes(x = 0, y = 0, xend = 1, yend = 1), color = "red", linetype = 2, linewidth = 1) +
  geom_segment(aes(x = 1, y = 1, xend = 2, yend = 1), color = "red", linetype = 2, linewidth = 1) +
  coord_cartesian(xlim = c(0,2)) + 
  geom_vline(xintercept = 0, linetype = 2, linewidth = 1) + 
  geom_hline(yintercept = 0, linetype = 5, linewidth = 1) + 
  geom_segment(aes(x = 1, y = 1, xend = 1.4142, yend = 0), color = "blue", linetype = 5, linewidth = 0.8) +
  geom_path(aes(color = YR, group = 1), linewidth = 1.2) +
  geom_point(aes(color = YR), size = 2) +
  scale_color_viridis_c(option = "mako") + 
  coord_fixed(ratio = 1, xlim = c(0, 2), ylim = c(0, 1)) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.2)) +
  theme_classic() +
  labs(x = "Ep/P", y = "E/P", color = "Year", title = "07340300 (34.38, -94.24)") +
  theme(plot.title = element_text(hjust = 0.5))



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

# site_locations <- ESA_groups %>% dplyr::select(SITENO, dec_lat_va, dec_long_va)

# sites as sf 
sites_sf <- st_as_sf(ESA_groups, coords = c("dec_long_va", "dec_lat_va"), crs = 4326)

# define extents 
# lon_min <- -105; lon_max <- -67
# lat_min <- 25; lat_max <- 49
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
# ESA clusters #
###########################################################
# sites as sf 
sites_sf$Partitions <- factor(sites_sf$Partitions, levels = 1:10)
cluster_colors <- c("#000000", "#999999", "#F0E442", "#56B4E9", "#009E73", 
                    "#CC79A7", "#0072B2", "#7F7F7F", "#D55E00", "#E69F00" )

cluster_shapes <- c(21, 22, 23, 24, 25, 21, 22, 23, 24, 25)


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
# groupings for distribution, Group E (MAD/0.6745)
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
# groupings for slope, Group S
sites_sf$Slope <- factor(sites_sf$Slope, levels = levels_5) 

slope <- ggplot() + geom_sf(data = states_crop, fill = NA, color = "black", linewidth = 0.4) +
  geom_sf(data = sites_sf, aes(fill = Slope), 
          shape = 21, color = "black", size = 1.5, stroke = 0.4, alpha = 0.9) +
  shared_colors +    
  # scale_fill_manual(values = cluster_colors, name = "Slope") +
  coord_sf(xlim = c(lon_min, lon_max),
           ylim = c(lat_min, lat_max),
           expand = FALSE, datum = st_crs(4326)) +
    labs(title = "Slope") + 
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
          size = 1.5, stroke = 0.4, alpha = 0.9) +
    labs(title = expression(ESA[rstd])) +
  shared_colors + clean_theme + coord_limits

tail_map <- ggplot() +
  geom_sf(data = states_proj, fill = NA, color = "black", linewidth = 0.4) +
  geom_sf(data = sites_proj, aes(fill = dESA), shape = 21, color = "black",
          size = 1.5, stroke = 0.4, alpha = 0.9) +
     labs(title = expression(ESA[tails])) +
  shared_colors + clean_theme + coord_limits

slope_map <- ggplot() +
  geom_sf(data = states_proj, fill = NA, color = "black", linewidth = 0.4) +
  geom_sf(data = sites_proj, aes(fill = Slope), shape = 21, color = "black",
          size = 1.5, stroke = 0.4, alpha = 0.9) +
   labs(title = "Slope (S)") +
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
###############################################
# plot with dem background

ggplot() +
  geom_raster(data = dem_df, aes(x = x, y = y, fill = elev), alpha = 0.8) +
  scale_fill_gradientn(
    name = "Elevation (m)",
    colors = c("darkgreen","palegreen","tan","saddlebrown"),
    limits = c(5, 2050),
    oob = scales::squish,
    na.value = NA) +
  # reset fill scale so clusters can use a discrete one
  ggnewscale::new_scale_fill() +
  geom_sf(data = states_crop, fill = NA, color = "black", linewidth = 0.4) +
  geom_sf(data = rivers_union, color = "blue", linewidth = 0.25) +
  geom_sf(data = sites_sf,
    aes(fill = Partitions, shape = Partitions),
    color = "black", size = 2, stroke = 0.4, alpha = 0.9) +
  scale_fill_manual(values = cluster_colors, name = "Partition") +
  scale_shape_manual(values = cluster_shapes, name = "Partition") +
  scale_x_continuous(breaks = seq(-100,-70, 10)) +
  scale_y_continuous(breaks = seq(25,50,5)) +
  coord_sf(
    xlim = c(lon_min, lon_max),
    ylim = c(lat_min, lat_max),
    expand = FALSE,
    datum = st_crs(4326)) +
  theme_minimal() +
  theme(axis.title = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()) #+ 
 # labs(title = "ESA Partitions on DEM")


#########################################################
partition_summary <- ESA_groups %>% 
   group_by(Partitions) %>% 
   summarise(ESA = mean(ESA), dESA = mean(dESA))

partition_long <- partition_summary %>% 
  pivot_longer(cols = c("ESA", "dESA"), 
               names_to = "Metric", values_to = "Group") %>% 
  mutate(Metric = factor(Metric, levels = c("ESA", "dESA"))) # remove alpha ordering

ggplot(partition_long, aes(x = Partitions, y = Metric, fill = factor(Group))) +
  geom_tile(color = NA) +   # removes white borders
  scale_fill_manual(values = c("#D55E00", "#F0E442", "#009E73", "#0072B2", "#CC79A7"),name = "Group") +
# scale_x_reverse(breaks = 1:10, expand = c(0,0), name = "Partition") +
  scale_x_continuous(breaks = 1:10, expand = c(0,0), name = expression(ESA[partition])) +
  scale_y_discrete(labels = c(expression(ESA[rstd]), expression(ESA[tails]))) +
  coord_fixed() +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(size = 10), 
    axis.title.y = element_blank(), 
    panel.grid = element_blank(),
    legend.position = "right")
  

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
  
 # mutate(DATE = as.Date("1980-01-01") + (row_number() - 1) * 183, 
  #       YR = year(DATE), 
        # Period = floor((YR - 1980) / 5) + 1), 
   #      Time_Index = row_number()) %>%  ungroup()

# limit to 1 (overestimation of E, greater than P)
esa_data$ESA <- ifelse(esa_data$ESA > 1.0, 1.0, esa_data$ESA)

# set at datatable type 
esa_data <- setDT(esa_data)

all_plots = list() # remove if saving individual cluster plots

# loop through clusters to create plots 
for (cl in 1:10) {
  # subset data for cluster
  cluster_data = esa_data[Partitions == cl]
  cluster_data[, Date := as.Date("1980-01-01") + (Time_Index - 1) * 183]

  # get median ESA value for each time period
  cluster_mean = cluster_data[, .(Med_ESA = median(ESA, na.rm = TRUE)), by = Time_Index]
  cluster_mean[, Date := as.Date("1980-01-01") + (Time_Index - 1) * 183]
  
  # Plot A: All basins' ESA time series
plot_a = ggplot(cluster_data, aes(x = Date, y = ESA, group = SITENO)) +
  geom_line(alpha = 0.3, linewidth = 0.3) +
  ylim(-1, 1.25) + 
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  scale_x_date(breaks = seq(as.Date("1980-01-01"), as.Date("2025-01-01"), by = "5 years"),
    labels = scales::date_format("%Y"),
    limits = as.Date(c("1980-01-01", "2025-01-01"))) +
  labs(title = paste("Partition", cl, ": ESA Time Series"), y = "ESA", x = NULL) + theme_minimal() + 
  theme(legend.position = "none", title = element_text(size = 8), axis.title = element_text(size = 6), 
        axis.text = element_text(size = 6), axis.text.x = element_text(angle = 45))
  
  # Plot B: Median ESA time series
  plot_b = ggplot(cluster_mean, aes(x = Date, y = Med_ESA)) +
    geom_line(color = "blue", linewidth = 0.75) + ylim(-1, 1.25) + 
    geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
      scale_x_date(breaks = seq(as.Date("1980-01-01"), as.Date("2025-01-01"), by = "5 years"),
    labels = scales::date_format("%Y"),
    limits = as.Date(c("1980-01-01", "2025-01-01"))) + theme_minimal() + 
   # labs(title = paste("Cluster", cl, ": Mean ESA Time Series"), y = "Mean ESA", x = NULL) +
    labs(title = NULL, y = "ESA", x = NULL) + 
  theme(legend.position = "none", axis.title = element_text(size = 6), 
        axis.text = element_text(size = 6), axis.text.x = element_text(angle = 45))
  
  # Combine plots
  # combined_plot = plot_a / plot_b  # for individual plots
  combined_plot = plot_a / plot_b + plot_layout(heights = c(2, 1))  # for all plots combined
  
  all_plots[[cl]] = wrap_elements(combined_plot)  # for all plots combined only
  
  # Save as PNG, individual plots
 # ggsave(filename = sprintf("cluster_plots_west/cluster_%02d_plots.png", cl),
    #     plot = combined_plot, width = 8, height = 8, dpi = 300)
}

# all plots together as one (12 plots)
final_clusters = wrap_plots(all_plots, ncol = 3, nrow = 4) +
  plot_annotation(tag_levels = 'a', tag_prefix = "", tag_suffix = ") ")

ggsave(filename = "F:/Maps/final_clusters.png", 
       plot = final_clusters, device = "png", width = 10, height = 16, dpi = 300)

######################################################
# all partitions in grid plot

# get unique basins per partition
n_df <- esa_data %>% 
  group_by(Partitions) %>% 
  summarise(n = n_distinct(SITENO)) %>% 
  arrange(Partitions)

# get partition levels
# partitions <- sort(unique(esa_data$Partitions))
partitions <- n_df$Partitions

# create labels (a), (b), etc
facet_labels <- setNames(paste0("(", letters[seq_along(partitions)], ") Partition ", partitions, 
                                " (n = ", n_df$n, ")"), partitions)

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

# tag_facet <- function(p, open = "(", close = ")", tag_fun = function(i) letters[i]) {
#  gb <- ggplot_build(p)
#  lay <- gb$layout$layout
  
#  tags <- paste0(open, tag_fun(seq_len(nrow(lay))), close)
  
#  p + geom_text(data = lay, aes(x = -Inf, y = Inf, label = tags), 
#                hjust = -0.3, vjust = 1.3, size = 3.5, inherit.aes = FALSE)
#}

# tag_facet(p)


###################################################################
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


sites_sf_2_3 <- sites_sf %>% filter(Partitions %in% c(2,3))

# location plot for each partition to inset into time series plot
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

# similar groups 
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

# two partitions
layout <- "
ABM
CDM
"

combined <- wrap_plots(
  A = plot_2,  B = plot_3,
  C = pdf_2,   D = pdf_3,
  M = map_partition_plot_2_3,
  design = layout,
  widths = c(2, 2, 1),
  heights = row_heights
) +
  plot_annotation(tag_levels = "a", tag_suffix = ") ")



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
# create 10-panel plot with Mean_ESA and dESAda
plot <- ggplot(deriv_switch_data, aes(x = Date)) +
  geom_line(aes(y = Med_ESA, color = "Median ESA"), linetype = "dashed") +
  geom_line(aes(y = dESAda, color = "dESAda")) +  # Scale dESAda for visibility
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  scale_x_date(breaks = seq(as.Date("1980-01-01"), as.Date("2023-07-01"), by = "5 years"),
               labels = scales::date_format("%Y"),
               limits = as.Date(c("1980-01-01", "2023-07-01"))) +
  scale_y_continuous(name = "ESA") +
  scale_color_manual(values = c("Median ESA" = "blue", "dESAda" = "black")) +
  facet_wrap(~ Cluster, ncol = 3, nrow = 4) +
  labs(title = "Median ESA and Derivative ESA Time Series by Partition (1980–2023)",
       x = "Year", color = "Series") +
  theme_minimal() +
  theme(legend.position = "bottom",
        axis.text = element_text(size = 6),
        axis.title = element_text(size = 6),
        strip.text = element_text(size = 6),
        axis.text.x = element_text(angle = 45, hjust = 1))


######################################################
### ESA plots ###
######################################################
# overall median value among each partition 
median_esa <- esa_data %>% group_by(Partitions, DATE) %>% 
  summarise(median = median(ESA))

median_partition <- esa_data %>% group_by(Partitions) %>% 
  summarise(median = median(ESA))


# determine median absolute deviation of ESA, robust standard deviation  
rstd_partition <- esa_data %>% group_by(Partitions) %>% 
  summarise(rstd = mad(ESA, center = median(ESA), constant = 1)/0.6745)


# derivative 
# use 88 median ESA values 
deriv13 <- function(f, dx) {
  n <- length(f)
  dfdx <- numeric(n)
  
  # Forward difference at start
  dfdx[1] <- (-3 * f[1] + 4 * f[2] - f[3])
  
  # Centered difference for interior points
  for (i in 2:(n-1)) {
    dfdx[i] <- (-f[i-1] + f[i+1])
  }
  
  # Backward difference at end
  dfdx[n] <- (f[n-2] - 4 * f[n-1] + 3 * f[n])
  
  # Scale by 2*dx
  dfdx <- dfdx / (2 * dx)
  
  return(dfdx)
}


# calculate derivative 
esa_desa <- median_esa %>% group_by(Partitions) %>% 
  mutate(dESA = deriv13(median, dx = 0.5))

dESA_partitions <- esa_desa %>% group_by(Partitions) %>% 
   summarise(desa_rstd = mad(dESA, center = median(dESA), constant = 1)/0.6745)

pct_tails_dESA <- esa_desa %>% group_by(Partitions) %>% 
  summarise(tails = mean(dESA < quantile(dESA, 0.25) - 1.5*IQR(dESA) |
                         dESA > quantile(dESA, 0.75) + 1.5*IQR(dESA)) *100)