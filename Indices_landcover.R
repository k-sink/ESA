library(tidyverse)

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
  summarise(CDD = mean(CDD), CSDI = mean(CSDI), CWD = mean(CWD), 
            FD = mean(FD), GSL = mean(GSL, na.rm = TRUE), ID = mean(ID), R00mm = mean(R00mm), 
            R01mm = mean(R01mm), R10mm = mean(R10mm), R20mm = mean(R20mm), 
            R95pTOT = mean(R95pTOT), R99pTOT = mean(R99pTOT), 
            Rx1day = mean(Rx1day), Rx5day = mean(Rx5day), SDII = mean(SDII), 
            SPEI = mean(SPEI), SPI = mean(SPI), SU = mean(SU), TR = mean(TR), 
            WSDI = mean(WSDI))

# get esa values
esa_values <- read.csv("F:/ESA_R/esa_manuscript.csv")
# get median annual value for each partition
median_esa <- esa_values %>% group_by(Partitions, YR) %>% summarise(ESA = median(ESA))

indices_partition <- left_join(indices_yr, median_esa, by = c("Partitions", "YR"))

ggplot(data = indices_partition) + geom_line(aes(x = YR, y = CDD)) + facet_wrap(~Partitions)
