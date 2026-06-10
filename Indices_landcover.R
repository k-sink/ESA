library(tidyverse)
library(broom)
library(patchwork)
library(SPEI)

# get group numbers for each ESA clustering
ESA_groups <-  read.table("F:/ESA/ESApargrps.txt", sep = "")
ESA_groups <-  ESA_groups %>% rename(SITENO = V1, Partitions = V2, Slope = V3, ESA = V4, dESA = V5)
ESA_groups$SITENO <- str_pad(as.character(ESA_groups$SITENO), width = 8, side = "left", pad = "0")

# get annual indices data from MACH
indices <- read.csv("D:/Journal submissions/Scientific Data/MACH_dataset/annual_climate.csv")
indices$SITENO <- str_pad(as.character(indices$SITENO), width = 8, side = "left", pad = "0")

# get the site numbers 
sites <- unique(ESA_groups$SITENO)

# filter indices for 747 sites
indices_filtered <- indices %>% filter(SITENO %in% sites)

# add grouping column 
indices_filtered <- indices_filtered %>% left_join(ESA_groups, by = "SITENO")

indices_yr <- indices_filtered %>% group_by(Partitions, YR) %>% 
  summarise(CDD = mean(CDD, na.rm = TRUE), CSDI = mean(CSDI, na.rm = TRUE), CWD = mean(CWD), 
            FD = mean(FD), GSL = mean(GSL, na.rm = TRUE), ID = mean(ID), R00mm = mean(R00mm), 
            R01mm = mean(R01mm), R10mm = mean(R10mm), R20mm = mean(R20mm), 
            R95pTOT = mean(R95pTOT), R99pTOT = mean(R99pTOT), 
            Rx1day = mean(Rx1day), Rx5day = mean(Rx5day), SDII = mean(SDII, na.rm = TRUE), 
            SPEI = mean(SPEI, na.rm = TRUE), SPI = mean(SPI, na.rm = TRUE), SU = mean(SU), TR = mean(TR), 
            WSDI = mean(WSDI))

# saved dataframe as csv, read file
indices_yr <- read_csv("F:/ESA_indices_climate.csv")

indices_long <- indices_yr %>%
  pivot_longer(cols = -c(Partitions, YR), names_to = "variable", values_to = "value")

trends <- indices_long %>%
  filter(!is.na(value)) %>%
  group_by(Partitions, variable) %>%
  summarise(
    slope = coef(lm(value ~ YR))[2],
    p_value = summary(lm(value ~ YR))$coefficients[2,4],
    .groups = "drop")

indices_plot <- indices_long %>% 
  left_join(trends, by = c("Partitions", "variable"))

partition_colors <- c(
  "1" = "#D55E00", "2" = "#E69F00", "3" = "#CC79A7", "4" = "#F0E442", "5" = "#009E73", 
  "6" =  "#56B4E9", "7" = "#0072B2", "8" = "#7F7F7F", "9" = "#999999", "10" =  "#000000")


plot_variable <- function(var_name) {
  data_var <- indices_plot %>% filter(variable == var_name)
  trends_var <- trends %>% filter(variable == var_name)

  ggplot(data_var, aes(x = YR, y = value, color = factor(Partitions))) +
    
    # time series
    geom_line(alpha = 0.6) +
    
    # trend lines
    geom_smooth(method = "lm", se = FALSE, linewidth = 1) +
    
    # label significant slopes
    geom_text(data = trends_var %>% filter(p_value < 0.05),
      aes(x = max(indices_plot$YR), y = Inf, label = paste0("m=", round(slope, 3))), 
      hjust = 1.1, vjust = 1.5,
      size = 3, inherit.aes = FALSE) + 
    
    scale_color_manual(values = partition_colors) +
     labs(title = paste("Trend for", var_name),
      color = "Partition",  y = var_name, x = "Year") + 
    theme_minimal()
}


plots <- lapply(unique(indices_plot$variable), plot_variable)

plots[[1]]

##
plot_variable <- function(var_name) {
  data_var <- indices_plot %>% filter(variable == var_name)
  trends_var <- trends %>% filter(variable == var_name)
  
  # Main time series plot
  p_main <- ggplot(data_var, aes(x = YR, y = value, color = factor(Partitions))) +
    geom_line(alpha = 0.6) +
    geom_smooth(method = "lm", se = FALSE, linewidth = 1) +
    scale_color_manual(values = partition_colors) +
    labs(title = paste("Trend for", var_name),
         color = "Partition", y = var_name, x = "Year") +
    theme_minimal()
  
  # Slope panel below
  p_slopes <- ggplot(trends_var, aes(x = factor(Partitions), y = slope, label = round(slope, 3), fill = factor(Partitions))) +
    geom_col(alpha = 0.6, show.legend = FALSE) +
    geom_text(vjust = -0.5, size = 3) +
    scale_fill_manual(values = partition_colors) +
    labs(x = "Partition", y = "Slope") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  # Combine vertically
  p_main / p_slopes + plot_layout(heights = c(3, 1))
}

# Generate plots for all variables
plots <- lapply(unique(indices_plot$variable), plot_variable)

# Show first plot
plots[[1]]

# get esa values
esa_values <- read.csv("F:/ESA_R/esa_manuscript.csv")
# get median annual value for each partition
median_esa <- esa_values %>% group_by(Partitions, YR) %>% summarise(ESA = median(ESA))

indices_partition <- left_join(indices_yr, median_esa, by = c("Partitions", "YR"))
slopes <- indices_partition %>% group_by(Partitions) %>% 
  summarise(model = list(lm(CDD ~ YR)), 
            slope = coef(model[[1]])[2],
            p_value = summary(model[[1]])$coefficients[2,4])
slopes <- slopes %>% mutate(label = paste0("Slope = ", round(slope, 3), 
                                           "\np = ", signif(p_value, 2)))

ggplot(data = indices_partition, aes(x = YR, y = CDD)) + geom_line() + 
  geom_smooth(method = "lm", se = FALSE, linewidth = 1) + 
  geom_text(data = slopes, aes(x = Inf, y = Inf, label = label), 
            hjust = 1.1, vjust = 1.5, size = 3, inherit.aes = FALSE) +
   theme_bw() + labs(title = "Consecutive dry days") +
  theme(axis.title.x = element_blank()) + facet_wrap(~Partitions, nrow = 5) 


###################################################################
## SPI ##
# standardized precipitation index
# drought index that captures how observed precipitation deviates from climatological average over given time period
# drought indices >0 wetter than normal, <0 drier than normal 
# <-2 extreme drought, <-1.5 to -2 severe drought
month = MACH_data %>% group_by(SITENO, YR, MNTH) %>% reframe(PRCP = sum(PRCP))
#month = month %>% mutate(date = as.Date(paste(YR, MNTH, "1", sep = "-")))
#month = month %>% dplyr::select(-c(YR, MNTH))

# spi and spei functions are identical except for distribution function and data, 
# time ordered values of precipitation (SPI) and gamma or climatic balance
# precipitation minus potential evapotranspiration (SPEI) and log-logistic
calc_spi = function(precip) {
  spi_result = spi(precip, scale = 1)
  return(spi_result$fitted)
}

spi_1 = month %>% group_by(SITENO) %>% mutate(SPI = calc_spi(PRCP)) %>% ungroup()
# remove -inf values which occur when month is zero total precipitation 
spi = spi_1 %>% filter(!is.infinite(SPI)) %>% group_by(SITENO, YR) %>% reframe(SPI = mean(SPI))
spi = pivot_wider(spi, names_from = YR, values_from = SPI)
write.csv(spi, "D:/University of Texas at Dallas/CombinedDataset/AttributesMACH/SPI.csv")

## SPEI (standardized precipitation evapotranspiration index) ##

# use water balance, wb = precip - pet for calculating index instead of precipitation 
month2 = MACH_data %>% group_by(SITENO, YR, MNTH) %>% reframe(PRCP = sum(PRCP), PET = sum(PET))
month2 = month2 %>% mutate(WB = PRCP - PET)

# remove site because no PET data available 
month2 = month2 %>% filter(SITENO != "10336740")

calc_spei = function(wb) {
  spei_result = spei(wb, scale = 1)
  return(spei_result$fitted)
}

spei_1 = month2 %>% group_by(SITENO) %>% mutate(SPEI = calc_spei(WB)) %>% ungroup()
# remove -inf values which occur when month is zero total precipitation 
spei = spei_1 %>% filter(!is.infinite(SPEI)) %>% group_by(SITENO, YR) %>% reframe(SPEI = mean(SPEI))
spei = pivot_wider(spei, names_from = YR, values_from = SPEI)
write.csv(spei, "D:/University of Texas at Dallas/CombinedDataset/AttributesMACH/SPEI.csv")

##############################
# Comparison with ESA

data_dir <- "D:/Journal submissions/Scientific Data/MACH_dataset/MACH_ts/ALLVAR/"
files <- list.files(data_dir, pattern = "^basin_.*_MACH\\.csv$", full.names = TRUE)

# function to process one file
process_file <- function(file) {
  
  # extract SITENO from filename
  siteno <- str_extract(basename(file), "\\d{8}")

  # monthly values for climate variables   
  read_csv(file, show_col_types = FALSE) %>%
    mutate(DATE = as.Date(DATE),
           YEAR = year(DATE),
           MONTH = month(DATE)) %>% 
    group_by(YEAR, MONTH) %>%
    summarise(PRCP = sum(PRCP, na.rm = TRUE), PET  = sum(PET,  na.rm = TRUE), 
              AET = sum(AET, na.rm = TRUE), TAIR = mean(TAIR, na.rm = TRUE), 
              VP = mean(VP, na.rm = TRUE), SRAD = mean(SRAD, na.rm = TRUE),  .groups = "drop") %>% 
    mutate(SITENO = siteno) %>% select(SITENO, YEAR, MONTH, PRCP, PET, AET, TAIR, SRAD, VP)
}

df_monthly <- map_dfr(files, process_file)

###################################
# calculate SPEI/SPI

calc_indices <- function(df_site) {
  start_year <- min(df_site$YEAR)

  prcp_ts <- ts(df_site$PRCP, start = c(start_year, 1), frequency = 12)
  bal_ts  <- ts(df_site$PRCP - df_site$PET, start = c(start_year, 1), frequency = 12)

  tibble(
    date = seq.Date(from = as.Date(paste0(start_year, "-01-01")),
                    by = "month",
                    length.out = length(prcp_ts)),
    SPI_12  = as.numeric(spi(prcp_ts, 12)$fitted),
    SPEI_3  = as.numeric(spei(bal_ts, 3)$fitted),
    SPEI_6  = as.numeric(spei(bal_ts, 6)$fitted),
    SPEI_12 = as.numeric(spei(bal_ts, 12)$fitted),
    SPEI_24 = as.numeric(spei(bal_ts, 24)$fitted)
  )
}

indices_df <- df_monthly %>%
  group_by(SITENO) %>% group_modify(~ calc_indices(.x)) %>% ungroup()

# create lagged spi/spei variables
indices_lagged <- indices_df %>%
  group_by(SITENO) %>% arrange(date) %>%
  mutate(
    # 3 month lags
    SPI12_lag3   = lag(SPI_12, 3),
    SPEI3_lag3   = lag(SPEI_3, 3),
    SPEI6_lag3   = lag(SPEI_6, 3),
    SPEI12_lag3  = lag(SPEI_12, 3),
    SPEI24_lag3  = lag(SPEI_24, 3),
    # 6 month lags
    SPI12_lag6   = lag(SPI_12, 6),
    SPEI3_lag6   = lag(SPEI_3, 6),
    SPEI6_lag6   = lag(SPEI_6, 6),
    SPEI12_lag6  = lag(SPEI_12, 6),
    SPEI24_lag6  = lag(SPEI_24, 6)) %>% ungroup()




esa_data <- read_csv("F:/ESA_R/esa_manuscript.csv") 
esa_data$DATE <- as.Date(esa_data$DATE, format = "%m/%d/%Y")        
# limit to 1 (overestimation of E, greater than P)
esa_data$ESA <- ifelse(esa_data$ESA > 1.0, 1.0, esa_data$ESA)

# rename columns
esa_data <-  esa_data %>% rename(P = PRCP, Ep = PET, E = AET) 
# get month and year components to match with SPEI, SPI indices 
esa_df <- esa_data %>% mutate(MNTH = month(DATE, label = FALSE))  
indices_df <- indices_df %>% mutate(MNTH = month(date, label = FALSE), YR = year(date))

esa_indices <- esa_df %>% left_join(indices_df, by = c("SITENO", "MNTH", "YR"))
esa_indices <- esa_indices %>% select(-c(P, Ep, E, EP, EpP, Time_Index, Slope, ESA_group, dESA))

# multi-month accumulation windows of SPI, SPEI result in undefined index values for initial portion
# of record (remove periods with NA, spei_3 and spei_6 will not have Jan 1980, spi_12 and spei_12 will not have
# any of 1980, spei_24 will not have 1980 or 1981)
spei_3_6 <- esa_indices %>% select(SITENO, ESA, DATE, Partitions, SPEI_3, SPEI_6) %>% drop_na()
spei_spi_12 <- esa_indices %>% select(SITENO, ESA, DATE, Partitions, SPI_12, SPEI_12) %>% drop_na()
spei_24 <- esa_indices %>% select(SITENO, ESA, DATE, Partitions, SPEI_24) %>% drop_na()

cor_spei_3 <- spei_3_6 %>% group_by(Partitions) %>% summarise(cor_spei3 = cor(ESA, SPEI_3, use = "complete.obs"))
cor_spei_6 <- spei_3_6 %>% group_by(Partitions) %>% summarise(cor_spei6 = cor(ESA, SPEI_6, use = "complete.obs"))
cor_spei_12 <- spei_spi_12 %>% group_by(Partitions) %>% summarise(cor_spei12 = cor(ESA, SPEI_12, use = "complete.obs"))
cor_spi_12 <- spei_spi_12 %>% group_by(Partitions) %>% summarise(cor_spi12 = cor(ESA, SPI_12, use = "complete.obs"))
cor_spei_24 <- spei_24 %>% group_by(Partitions) %>% summarise(cor_spei24 = cor(ESA, SPEI_24, use = "complete.obs"))

#######
