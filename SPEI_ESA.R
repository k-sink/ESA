#############################################
# SPEI INDEX
#############################################
# Katharine Sink

# calculate SPEI for 12 month period, sample at same points as 88 ESA time series
# compare ESA time series to SPEI index

library(tidyverse)
library(SPEI)
library(stats)

#######################
# MONTHLY TOTALS
#######################
# daily data files 
data_dir <- "F:/ESA/ESA_2026"
files <- list.files(data_dir, pattern = "^basin_.*_MACH\\.csv$", full.names = TRUE)

# function to process one file, calculate monthly totals for prcp and pet
process_file <- function(file) {
  
  # extract SITENO from filename
  siteno <- str_extract(basename(file), "\\d{8}")

  # monthly values for climate variables   
  read_csv(file, show_col_types = FALSE) %>%
    mutate(DATE = as.Date(DATE),
           YR = year(DATE),
           MNTH = month(DATE)) %>% 
    group_by(YR, MNTH) %>%
    summarise(PRCP = sum(PRCP, na.rm = TRUE), PET  = sum(PET,  na.rm = TRUE),  .groups = "drop") %>% 
    mutate(SITENO = siteno) %>% select(SITENO, YR, MNTH, PRCP, PET)
}

df_monthly <- map_dfr(files, process_file)

# monthly balance
df_monthly <- df_monthly %>% mutate(BAL = PRCP - PET, DATE = as.Date(paste(YR, MNTH, "01", sep = "-")))

# get record check for completeness (528 months expected)
record_check <- df_monthly %>%
  group_by(SITENO) %>%
  summarise(
    n_months = n(),
    start_date = min(DATE),
    end_date = max(DATE),
    .groups = "drop"
  )

print(summary(record_check$n_months))

######################
# SPEI CALCULATION
######################
# calculate SPEI for each site 
# timescales selected to match ESA equivalent width of ~ 1 year
# SPEI-12: 12 month water balance ~annual water-energy deficit
# SPEI-24: 24 month accumulation ~multi-year drought signal 
# SPEI-3 and SPEI-6 timescales are not compatible with seminannual ESA 
# (FIR filter EW ~ 1 yr suppresses variability at these scales) 

calc_spei <- function(df_site) {

  start_year <- min(df_site$YR)
  start_month <- min(df_site$MNTH[df_site$YR == start_year])
  n_months <- nrow(df_site)

  # create time series objects
  prcp_ts <- ts(df_site$PRCP, start = c(start_year, start_month), frequency = 12)
  balance_ts <- ts(df_site$BAL, start = c(start_year, start_month), frequency = 12)
  spei_12 <- tryCatch(
    as.numeric(spei(balance_ts, scale = 12, na.rm = TRUE)$fitted),
    error = function(e) rep(NA_real_, n_months)
  )
  tibble(
    DATE = df_site$DATE,
    YR = df_site$YR,
    MNTH = df_site$MNTH,
    PRCP = df_site$PRCP,
    PET = df_site$PET,
    BALANCE = df_site$BAL,
    SPEI_12 = spei_12
  )
}

indices_monthly <- df_monthly %>%
  group_by(SITENO) %>%
  group_modify(~ calc_spei(.x)) %>%
  ungroup()

######################################################
# read in ESA data 
esa_data <- read_csv("esa_data.csv") %>%
  mutate(
    DATE = as.Date(DATE, format = "%m/%d/%Y"),
    MNTH = month(DATE),
    # cap ESA at 1.0 (overestimation where E > P)
    ESA = pmin(ESA, 1.0)
  )

print(esa_data %>%
        group_by(SITENO) %>%
        summarise(n = n(), .groups = "drop") %>%
        pull(n) %>%
        summary())

# join time series with partition assignment for each watershed
esa_data <- esa_data %>% left_join(kmeans_r_results %>% dplyr::select(SITENO, Partition), by = "SITENO")

# join SPEI to ESA at semiannual time points
# sample indices only at ESA decimated time points (approx Jan and July)
# ensures both signals represent same temporal window and avoid mixing timescales
esa_indices <- esa_data %>%
  select(SITENO, DATE, YR, MNTH, ESA, Partition) %>%
  left_join(
    indices_monthly %>%
      select(SITENO, YR, MNTH, SPEI_12),
    by = c("SITENO", "YR", "MNTH")
  )
# check join 
join_check <- esa_indices %>%
  summarise(
    n_total = n(),
    n_spei_na = sum(is.na(SPEI_12))
  )

# remove initial NA values 
# SPEI 12 undefined for first 11 months
esa_spei12 <- esa_indices %>%
  select(SITENO, DATE, YR, MNTH, ESA, Partition, SPEI_12) %>% drop_na()

# write csv file 
write.csv(esa_spei12, here("esa_spei12.csv"), row.names = FALSE)

#######################################
# CORRELATION BETWEEN ESA AND SPEI
#######################################
# partition level correlations
# pearson R between ESA and index across all basins and semiannual time points 
# within each partition 

cor_12 <- esa_spei12 %>%
  group_by(Partition) %>%
  summarise(
    n = n(),
    n_watersheds = n_distinct(SITENO), 
    r_SPEI_12 = cor(ESA, SPEI_12, use = "complete.obs"),
    .groups = "drop")

write.csv(cor_12, here("spei_esa_partion_cor.csv"), row.names = FALSE)

#########################################
# watershed level correlation 
# one correlation calculated separately for each watershed, 
# then summarized within each partition, watersheds are independent units and 
# summarize variability among sites

site_cor_12 <- esa_spei12 %>% 
  group_by(Partition, SITENO) %>% 
  summarise(
    r_SPEI_12 = cor(ESA, SPEI_12, use = "complete.obs"), 
    .groups = "drop"
  )

write.csv(site_cor_12, here("spei_esa_watershed_cor.csv"), row.names = FALSE)

# summarize watershed-level correlations by partition 
cor_watershed <- site_cor_12 %>% 
  group_by(Partition) %>% 
  summarise(
    
    # number of independent watersheds
    n_watersheds = n(), 
    
    # descriptive stats for watershed correlations
    mean_r_SPEI_12 = mean(r_SPEI_12, na.rm = TRUE), 
    median_r_SPEI_12 = median(r_SPEI_12, na.rm = TRUE), 
    sd_r_SPEI_12 = sd(r_SPEI_12, na.rm = TRUE), 
    
    # Fisher z transformation because correlation coeff bounded between +/- 1
    mean_z = mean(
      atanh(r_SPEI_12), 
      na.rm = TRUE), 
    
    sd_z = sd(
      atanh(r_SPEI_12), 
      na.rm = TRUE), 
    
    # backtransformed mean Fisher z
    mean_r_from_z = tanh(mean_z), 
    
    # one sample t-test of Fisher z values against zero
    t_stat = mean_z / (sd_z / sqrt(n_watersheds)), 
    df = n_watersheds - 1, 
    p_value = 2 * pt(-abs(t_stat), df = df), 
    
    # 95% confidence interval for mean Fisher z
    z_lower = mean_z - qt(0.975, df = df) * sd_z / sqrt(n_watersheds), 
    z_upper = mean_z + qt(0.975, df = df) * sd_z / sqrt(n_watersheds),
    .groups = "drop") %>% 
  
  # transform confidence limits to correlation scale 
  mutate(
    r_lower = tanh(z_lower), 
    r_upper = tanh(z_upper), 
    
    significance = case_when(
      p_value < 0.001 ~ "***", 
      p_value < 0.01 ~ "**", 
      p_value < 0.05 ~ "*", 
      TRUE  ~ "")) %>% 
 
  arrange(Partition)

# summary table with partition and watershed level correlations between ESA and SPEI
supplement_table <- cor_12 %>%
  left_join(
    cor_watershed %>%
      select(Partition, 
             mean_r_SPEI_12, median_r_SPEI_12, sd_r_SPEI_12,
             r_lower, r_upper, p_value, significance), by = "Partition") %>% 
  arrange(Partition)

write.csv(supplement_table, here("spei_cor_summary.csv"), row.names = FALSE)

#########################################
# LAG COMPARISON BETWEEN SPEI AND ESA
#########################################
# confirm direct comparison, not that SPEI lags ESA or vice versa
# positive lag: shift SPEI forward (SPEI leads ESA)
# negative lag: shift SPEI backward (ESA leads SPEI)

lag_months <- c(-6, -3, 0, 3, 6)

lag_cor_results <- map_dfr(lag_months, function(lag) {

  lagged_df <- esa_spei12 %>%
    group_by(SITENO) %>%
    arrange(DATE) %>%
    mutate(
      SPEI_12_lagged = lag(SPEI_12, n = abs(lag), default = NA) *
                       ifelse(lag >= 0, 1, 1),
      SPEI_12_lagged = if (lag < 0) {
        lead(SPEI_12, n = abs(lag), default = NA)
      } else {
        lag(SPEI_12,  n = abs(lag), default = NA)
      }
    ) %>%
    ungroup()

  lagged_df %>%
    group_by(Partition) %>%
    summarise(
      lag_months = lag,
      r_SPEI_12  = cor(ESA, SPEI_12_lagged, use = "complete.obs"),
      .groups = "drop"
    )
})

# find lag of maximum absolute correlation per partition
best_lag <- lag_cor_results %>%
  group_by(Partition) %>%
  slice_max(abs(r_SPEI_12), n = 1) %>%
  ungroup()

######################################################
### INDICES PLOTS ###
######################################################

# grouped bar chart with r values by partition
cor_long <- site_cor_12 %>%
  select(Partition, r_SPEI_12) %>%
  pivot_longer(
    cols = starts_with("r_"),
    names_to = "Index",
    values_to = "r") %>%
  mutate(
    Index = recode(Index,
      "r_SPEI_12" = "SPEI-12"),
    Partition = factor(Partition)
  )

indices_cor_plot <- ggplot(cor_long,
       aes(x = Partition, y = r, fill = Index)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6) +
  geom_hline(yintercept = 0, color = "black", linewidth = 0.5) +
  scale_fill_manual(
    values = c("SPEI-12"  = "steelblue")) +
  scale_y_continuous(
    limits = c(-0.8, 0),
    breaks = seq(-0.8, 0, 0.2),
    labels = seq(-0.8, 0, 0.2)) +
  labs(
  #  title = "Pearson Correlation between ESA and SPEI by Partition",
    x = "ESA Partition",
    y = "Pearson r",
    fill  = "Index") +
  theme_minimal(base_size = 11) +
  theme(
    legend.position  = "bottom",
    panel.grid.minor = element_blank(),
    plot.title = element_text(size = 10, face = "bold")
  )

ggsave(filename = "F:/ESA/indices_cor_plot.png", 
       plot = indices_cor_plot, device = "png", width = 10, height = 8, dpi = 300)

###############################################################
# compute partition median ESA and mean SPEI-12 at each time point
# for representative partitions
rep_partitions <- c(1, 7, 10)

ts_plot_data <- esa_spi12_spei12 %>%
  filter(Partitions %in% rep_partitions) %>%
  group_by(Partitions, DATE) %>%
  summarise(
    median_ESA = median(ESA, na.rm = TRUE),
    mean_SPEI_12 = mean(SPEI_12, na.rm = TRUE),
    mean_SPI_12 = mean(SPI_12, na.rm = TRUE),
    .groups = "drop") %>%
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
               color = "red",
               linewidth = 0.8) +
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
      name = "Median ESA",
      limits = c(-1.2, 1.2),
      breaks = seq(-1, 1, 0.5),
      sec.axis = sec_axis(
        transform = ~ . / scale_factor,
        name = "Mean SPEI-12",
        breaks = seq(-3, 3, 1))) +
    scale_color_manual(
      values = c("ESA"     = "#593196",
                 "SPEI-12" = "#d6604d")) +
    scale_x_date(
      date_breaks = "5 years",
      date_labels = "%Y") +
    labs(
      title = paste("Partition", partition_num),
      x = NULL,
      color = NULL) +
    theme_minimal() +
    theme(
      legend.position  = if (show_legend) "bottom" else "none",
      axis.title = element_text(size = 10),
      axis.text = element_text(size = 10),
      panel.grid.minor = element_blank(),
      plot.title = element_text(size = 10, hjust = 0.5))

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
