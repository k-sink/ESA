#################################################################
### PARTITION ASSIGNMENT AND STATISTICS FOR ESA TIME SERIES ###
#################################################################
# Katharine Sink

# each watershed (n = 747, 88 observations) is treated as a single vector in 88 dimensional space
# k-means groups watersheds whose entire ESA time series are similar, using squared Euclidean distance
# as measure of similarity, each partition's centroid is the element-wise mean ESA trajectory of 
# all watersheds assigned to it (average value at each of the 88 time steps across cluster members)
# similarity of pattern grouping, not location (spatial coherence is expected but no geographic information included)

library(tidyverse)
library(ClusterR) # kmeans
library(spdep) # spatial coherence
library(geosphere)
library(here)
library(patchwork)

######################################################
### K-MEANS CLUSTERING ###
######################################################  
# Matlab code = clustanalkm.m 
# caveat: results will not be identical due to inherent implementation differences between R and matlab
# MATLAB (Mathworks) vs R (ClusterR)
# https://www.mathworks.com/help/stats/kmeans.html
# https://www.rdocumentation.org/packages/ClusterR/topics/KMeans_rcpp

# Feature             MATLAB                           R (ClusterR::KMeans_rcpp)
# Default distance    Squared Euclidean                Squared Euclidean
# Initialization      k-means++ (default)              k-means++
# Replicates          1 (default)                      1 (default)
# Iterations          100 (default)                    10 (default)
# Implementation      Compiled                         Compiled (RcppArmadillo)  
# Output centers      Yes (C)                          Yes ($centroids)
# Distances (D)       n x k, squared euclidean         No direct equivalent
# Empty clusters      Handles                          Handles differently

# Distinguish nomenclature:
# 1. ITERATION (max_iters): happens inside a single kmeans attempt, starting from one set of centroids
#   the algorithm (a) assigns every watershed to its nearest current centroid and (b) recomputes each centroid
#   as the mean of its newly assigned watersheds, until assignments stop changing (convergence) or max_iters is 
#   reached, happens the same way every single time KMeans_rcpp() is called
# 2. INITIALIZATION/REPLICATION (num_init = replicates in Matlab): because WHERE you start determines which of several 
#   possible stopping points (local optima) you land on, KMeans tries num_init different starting points, lets each 
#   one run through the iteration independently, to its own convergence, and automatically keeps the attempt that 
#   ended with the lowest total WCSS, happens internally and the final best result is yielded 
# 3. SEED: fixes which specific random numbers get generated, so num_init starting points are reproducible, changing 
#   seed does not run more replicates, it just changes which starting points get tried 

#################################
# FUNCTIONS  
#################################
# helper function to build ESA matrix 
# converts long format data frame into wide matrix 
# rows = watersheds (SITENO), columns = time steps (DATE)

build_ESA_matrix <- function(ESA_df) {
  
  ESA_df <- ESA_df %>% dplyr::arrange(SITENO, DATE)
  
  # check for duplicate siteno/date combination
  dupes <- ESA_df %>% dplyr::count(SITENO, DATE) %>% dplyr::filter(n > 1)
  if (nrow(dupes) > 0) {
    print(head(dupes, 20))
    stop(nrow(dupes), "duplicate SITENO/DATE found")
  }
  
  # check for same number of observations for each site
  site_counts <- ESA_df %>% dplyr::count(SITENO)
  if (length(unique(site_counts$n)) > 1) {
    print(table(site_counts$n))
    stop("Differing number of observations")
  }
  
  S <- ESA_df %>%
    dplyr::select(SITENO, DATE, ESA) %>%
    tidyr::pivot_wider(names_from = DATE, values_from = ESA) %>%
    tibble::column_to_rownames("SITENO") %>%
    as.matrix()
  
  # check for NA
  if (any(is.na(S))) {
    na_rows <- rownames(S)[rowSums(is.na(S)) > 0]
    stop("NA values found")
  }
  S
}

########################################################
# function to calculate centroids, squared distances, RMSE
# centroid: mean ESA trajectory of whichever watersheds are actually assigned to each partition
# centroid is the mean of its members

compute_partition_geometry <- function(S, par, npar) {

  mnpar <- t(sapply(1:npar, function(i) colMeans(S[par == i, , drop = FALSE])))

  # squared Euclidean distance from every stite to every centroid
  ss_S <- rowSums(S^2)
  ss_C <- rowSums(mnpar^2)
  D <- outer(ss_S, ss_C, "+") - 2 * S %*% t(mnpar)
  D[D < 0] <- 0

  rmse <- sapply(1:npar, function(i) {
    idx <- which(par == i)
    sqrt(mean(D[idx, i]^2))
  })

  list(mnpar = mnpar, D = D, rmse = rmse)
}

#################################
# KMEANS PARTITIONING FUNCTION
#################################
# runs k-means++ seeded clustering (ClusterR::KMeans_rcpp) on the ESA matrix, 
# reorders partitions by decreasing median ESA (median of every ESA value across all watersheds and time steps), 
# and computes per-partition RMSE 

clustanalkm <- function(ESA_df, npar, num_init = 1000, seed = 1) {
 
  S <- build_ESA_matrix(ESA_df)
  dates <- sort(unique(ESA_df$DATE))

  km <- KMeans_rcpp(
    data = S,
    clusters = npar,
    num_init = num_init,   
    max_iters = 100,       
    initializer = "kmeans++",
    seed = seed)
 
  par <- km$clusters # initial unordered partition assignments 

  # reorder by decreasing pooled median ESA
  raw_partition_lookup <- data.frame(SITENO = rownames(S), raw_par = par)
  pooled_median_cluster <- ESA_df %>%
    dplyr::left_join(raw_partition_lookup, by = "SITENO") %>%
    dplyr::group_by(raw_par) %>%
    dplyr::summarise(pooled_ESA_median = median(ESA, na.rm = TRUE), .groups = "drop")
  
  mdesa <- pooled_median_cluster$pooled_ESA_median[match(1:npar, pooled_median_cluster$raw_par)]
  ord <- order(mdesa, decreasing = TRUE)
 
  oldpar <- par
  for (i in seq_along(ord)) {
    par[oldpar == ord[i]] <- i
  }
 
  geom <- compute_partition_geometry(S, par, npar)
 
  list(
    SITENO = rownames(S),
    dates = dates,
    par = par,
    mnpar = geom$mnpar,
    D = geom$D,
    rmse = geom$rmse,
    ESA_matrix = S,
    WCSS = sum(km$WCSS_per_cluster)
  )
}

##########################################
# load ESA data (747 watersheds, 88 ESA observations per watershed)
# output from ESA_calculation.R 
# ESA is theoretically bounded but can slightly exceed 1.0 due to overestimated AET
# limit ESA to range of -1 to +1

esa_data <- read_csv(here("results", "esa_data.csv"))
esa_data$ESA <-  ifelse(esa_data$ESA > 1.0, 1.0, esa_data$ESA)

esa_df <- esa_data %>% select(SITENO, DATE, ESA) 

##########################################
# JUSTIFY NUMBER OF PARTITIONS (k)
##########################################
# evaluate partition values (k = 2, 4, ... 20) using three metrics
# WCSS (total_SSE): within cluster sum of squares, decreases as k increases, the elbow
#   indicates the point past which additional partitions yield little improvement
# Silhouette width: measures how well separated and internally cohesive clusters are, ranges
#   from -1 to 1, higher is better, can decrease or increase with k so helps to identify peak/plateau 
# Marginal gain in cohesion: percent reduction in WCSS from previous k, reflects gain per two 
#   additional partitions (using even number of k), large early drops followed by small, flattening drops
#   indicate little improvement

# create matrix for diagnostics, clustanalkm() builds copy 
S <- build_ESA_matrix(esa_df)

# partition values tested
k_vals <- seq(2, 20, by = 2)

# distance object for silhouette scoring, computed once and reused across all k values tested
# uses standard (non-squared) Euclidean distance, consistent with typical silhouette-score practice.
dist_S <- dist(S)
 
# use num_init = 100 for selecting k, not quality stability 
k_selection <- lapply(k_vals, function(k) {
  km_k <- KMeans_rcpp(S, clusters = k, num_init = 100, max_iters = 100,
                       initializer = "kmeans++", seed = 1)
  sil <- cluster::silhouette(km_k$clusters, dist_S)
  data.frame(
    k = k,
    WCSS = sum(km_k$WCSS_per_cluster),
    avg_silhouette = mean(sil[, "sil_width"])
  )
})

k_selection <- do.call(rbind, k_selection)
 
# marginal gain: percent reduction in WCSS relative to the previous k tested
k_selection$marginal_gain_pct <- c(NA, -diff(k_selection$WCSS) / head(k_selection$WCSS, -1) * 100)

# SUPPLEMENTARY TABLE (S1)  
print(k_selection)

# SUPPLEMENTARY PLOT (FIGURE S1)
# elbow plot: WCSS vs k
wcss_plot <- ggplot(k_selection, aes(x = k, y = WCSS)) +
  geom_line() + geom_point() +
  geom_vline(xintercept = 10, linetype = "dashed", color = "red") +
  labs(title = "Within-Cluster Sum of Squares (WCSS)",
       x = "Number of Partitions (k)", y = "Total WCSS (squared error)") +
  theme_classic()
 
# silhouette plot: average silhouette width vs k
silo_plot <- ggplot(k_selection, aes(x = k, y = avg_silhouette)) +
  geom_line() + geom_point() +
  geom_vline(xintercept = 10, linetype = "dashed", color = "red") +
  labs(title = "Average Silhouette Width by Number of Partitions",
       x = "Number of Partitions (k)", y = "Average Silhouette Width") +
  theme_classic()

# combined plots for panel figure
ESA_cluster_validation <- wcss_plot + silo_plot + 
  plot_layout(ncol = 2) + plot_annotation(tag_levels = "a", tag_prefix = "(", tag_suffix = ") ") &
  theme(plot.tag = element_text(size = 10))

ggsave(filename = here("Figures/ESA_cluster_validation.png"), 
       plot = ESA_cluster_validation, device = "png", width = 10, height = 6, dpi = 300)

####################################################
# JUSTIFY num_init and CONFIRM CLUSTERING STABILITY
####################################################
# sensitivity of clustering to number of initializations (num_init) at k = 10
# for each candidate num_init (replicates), five independent clustering runs (random seed) were compared
# sd_WCSS near 0 means every seed's best attempt finds same objective value
# mean_ARI / min_ARI (Adjusted Rand Index) near 1 means every seed finds same actual cluster assignments

# 10 partitions (k) selected 
npar <- 10

n_reps <- 5  
num_init_vals <- c(10, 25, 50, 100, 250, 500, 1000)
 
num_init_diagnostic <- lapply(num_init_vals, function(n) {
  runs <- lapply(1:n_reps, function(r) {
    KMeans_rcpp(S, clusters = npar, num_init = n, max_iters = 100,
                initializer = "kmeans++", seed = r)
  })
 
  wcss <- sapply(runs, function(km) sum(km$WCSS_per_cluster))
  clusters_list <- lapply(runs, function(km) km$clusters)
 
  # mean and minimum pairwise ARI across all n_reps runs at this num_init
  pairs <- combn(1:n_reps, 2)
  ari_vals <- apply(pairs, 2, function(ij) {
    mclust::adjustedRandIndex(clusters_list[[ij[1]]], clusters_list[[ij[2]]])
  })
 
  data.frame(
    num_init = n,
    best_WCSS = min(wcss),
    sd_WCSS = sd(wcss),
    mean_ARI = mean(ari_vals),
    min_ARI = min(ari_vals)
  )
})
num_init_diagnostic <- do.call(rbind, num_init_diagnostic)
print(num_init_diagnostic)

# sd_WCSS typically flattens before min_ARI does, gaps reflect genuine near ties between watersheds 
# and two competing centroids, not insufficient num_init

chosen_num_init <- 1000

##############################################
# FINAL PARTITION ASSIGNMENT - MAJORITY VOTE
#############################################
# final partition assignment for each watershed using most frequent assignment across independent runs
# each site receives confidence score with 1 for unanimous, lower for tied sites
# does not depend on results from a single seed (starting point)
# consistency check (majority vote across many seeds)

# independent full clustering runs, directly identifies any contested watersheds (located near-equidistant
# between two centroids), far more reliably than comparing two runs against each other
# higher than n_reps because vote_confidence near 0.5 needs more samples to be trustworthy 
n_votes <- 200 

# run clustering 200 times, using seeds 1 through 200
vote_runs <- lapply(1:n_votes, function(r) {
  clustanalkm(esa_df, npar = 10, num_init = 1000, seed = r)$par
})
 
# stack all 200 runs into one table (one row per watershed, one column per run)
vote_matrix <- do.call(cbind, vote_runs)
rownames(vote_matrix) <- rownames(S)
 
# mode for each watershed (partition assignment)
get_mode <- function(x) as.integer(names(sort(table(x), decreasing = TRUE))[1])
 
# determine majority assignment and confidence in assignment
final_partition <- apply(vote_matrix, 1, get_mode)
top_votes <- apply(vote_matrix, 1, function(x) max(table(x)))
vote_confidence <- top_votes / n_votes

# 95% wilson score confidence interval on each site
wilson_ci <- t(sapply(top_votes, function(x) {
  ci <- prop.test(x, n_votes, correct = FALSE)$conf.int
  c(lower = ci[1], upper = ci[2])
}))

# centroids, distance, RMSE recomputed from FINAL partition assignment (majority vote)
final_geom <- compute_partition_geometry(S, final_partition, npar)

# summary data frame for each watershed and partition assignment
kmeans_results <- data.frame(
  SITENO = rownames(S),
  Partition = final_partition,
  Vote_confidence = vote_confidence,
  CI_lower = wilson_ci[, "lower"],
  CI_upper = wilson_ci[, "upper"], 
  Distance = final_geom$D[cbind(seq_len(nrow(final_geom$D)), final_partition)], 
  RMSE = final_geom$rmse[final_partition]
)
 
# for each contested site, find its two most-voted partitions across the 200 runs
# anything without unanimous agreement across all n_votes
contested_pairs <- kmeans_results %>%
  filter(Vote_confidence < 1) %>%
  rowwise() %>%
  mutate(
    vote_counts = list(sort(table(vote_matrix[SITENO, ]), decreasing = TRUE)),
    top_partition = as.integer(names(vote_counts)[1]),
    second_partition = as.integer(names(vote_counts)[2])
  ) %>%
  ungroup() %>%
  select(SITENO, Partition, Vote_confidence, top_partition, second_partition)

write.csv(contested_pairs, here("contested_pairs.csv"), row.names = FALSE)

###########################
# PARTITION LEVEL RMSE 
###########################
# RMSE within a partition
# how far each watershed within a partition (88 ESA time series) sits
# from partition centroid (mean time series of all watersheds within the partition)
# kmeans assigns every watershed to centroid that squared distance is closest to and iterates to reduce
# total sum of squared distances (WCSS) across every watershed and every partition simultaneously

partition_rmse_summary <- data.frame(
  Partition = 1:npar,
  RMSE = final_geom$rmse,
  n_watersheds = as.integer(table(factor(final_partition, levels = 1:npar)))
)
 
print(partition_rmse_summary)

##############################
# WATERSHED LEVEL TABLE
##############################
# table with 747 SITENO, partition (1 - 10), coordinates

site_info <- read_csv("F:/ESA/site_info.csv")

# table with attributes for each watershed and assigned partition
partition_attributes <- kmeans_results %>% 
  left_join(site_info %>% dplyr::select(SITENO, huc_cd, station_name, state, 
                                        dec_lat_va, dec_long_va, area_sqkm, elev_mean))

write.csv(partition_attributes, here("partition_attributes.csv"), row.names = FALSE)

# kmeans results with coordinates and partitions
kmeans_r_results <- kmeans_results %>% 
  left_join(site_info %>% dplyr::select(SITENO, dec_lat_va, dec_long_va), by = "SITENO")

# format siteno column as character with leading zero if needed
kmeans_r_results <- kmeans_r_results %>% mutate(SITENO = sprintf("%08d", SITENO))

write.csv(kmeans_r_results, here("kmeans_r_results.csv"), row.names = FALSE)

###############################################
### PARTITION LEVEL ESA STATISTICS ###
###############################################
# calculate statistics by partition using 88 ESA time series (pooled)
# includes ESA, dESA (rate of change) 
# Matlab code = clustchar.m 

# data frame with SITENO, DATE, ESA and dESA values, join with partition assignment 
partition_data <- esa_data %>% dplyr::select(SITENO, DATE, ESA, dESA) %>% 
  dplyr::left_join(kmeans_r_results %>% dplyr::select(SITENO, Partition), by = "SITENO") %>% 
  arrange(SITENO, DATE)

#######################################################################
# function (Matlab code =  fndout.m) to find outliers in the distribution of random variable in array y
# using Tukey (1977), outliers are improbable values (5% tails of Gaussian distribution) \
# definition based on integer depth into the sorted data, not interpolated percentiles
fndout_R <- function(y){
 
  # remove missing values so there is no interference with the sorting or hinge calculation but
  # remember where they are so final output can be same length as original input y 
  pnotnan <- !is.na(y)  # TRUE/FALSE flag 
  yy <- y[pnotnan]  # yy = y with all NA removed (shorter vector if applicable)
  n <- length(yy) # length of valid (non-NA) values
  
  # get rank order of data, not sorted values
  # vector of indices (p[1] indicates which position in yy holds smallest value, then p[2], etc)
  p <- order(yy) 
  
  # compute depth positions (how many steps in from either of the sorted data to look for median and each hinge)
  # Tukey's integer based depth definition, not interpolated percentile like quantile() would use
  k <- floor((n + 1) / 2) # median depth: position of the middle value (from either end)
  j <- floor(k / 2) + 1 # hinge depth: position of the lower/upper quartile like hinge
  
  # find the two hinges (Tukey's version of Q1 and Q3)
  # each hinge is the average of two specific order statistics, the value at depth j counted in from the low end, 
  # and the value at the mirror image depth counted in from the high end
  # averaging two values makes this well-defined, even when the depth j does not land exactly on single data point
  hlo <- (yy[p[j]] + yy[p[k + 1 - j]]) / 2    # lower hinge (~ Q1)
  hhi <- (yy[p[n - k + j]] + yy[p[n + 1 - j]]) / 2  # upper hinge (~Q3)
  
  # spread between the hinges is Tukey's version of IQR, data values delimiting inner 50% of values
  hsprd <- hhi - hlo  
  
  # step is how far beyond each hinge you have to go before a value counts as an outlier
  # 1.5x the hinge spread is Tukey's standard outlier rule (same 1.5 multiplier as familiar boxplot/IQR rule, 
  # just applied to hinges instead of quantiles)
  step <- 1.5 * hsprd
  
  # fences are actual outlier cutoff boundaries, anything beyond is flagged as an outlier
  fenlo <- hlo - step  # low fence 
  fenhi <- hhi + step  # high fence 
  
  # build TRUE/FALSE flag vectors, one entry per original value in y (including NA), all starting as FALSE by default
  poutlo <- rep(FALSE, length(y))  # array of logical pointers to low outliers
  pouthi <- rep(FALSE, length(y))  # array of logical pointers to high outliers 
  
  # fill in real answers, but only at non-NA positions (pnotnan), this is what lets poutlo/pouthi come back the 
  # same length as original y, with NA positions safely left as FALSE, while actual 
  # comparison (yy < fenlo) uses NA-free data
  poutlo[pnotnan] <- yy < fenlo # TRUE wherever a value falls below the low fence
  pouthi[pnotnan] <- yy > fenhi # TRUE wherever a value falls above the high fence

  # return list with flags and numeric fence values, quantities, etc
  list(
    poutlo = poutlo,
    pouthi = pouthi,
    fenlo = fenlo,
    fenhi = fenhi,
    step = step,
    hsprd = hsprd
  )
}

###############################################################
# summary statistics by partition 
# used to create TABLE 1

partition_stats <- partition_data %>%
  group_by(Partition) %>%
  group_modify(~{

    # ESA
    esa <- .x$ESA
    out_esa <- fndout_R(esa)
    n_esa <- sum(!is.na(esa))

    ESA_low_fence  <- out_esa$fenlo
    ESA_high_fence <- out_esa$fenhi
    ESA_low_out_pct  <- 100 * sum(out_esa$poutlo, na.rm = TRUE) / n_esa
    ESA_high_out_pct <- 100 * sum(out_esa$pouthi, na.rm = TRUE) / n_esa

    esa_tail_flag <- out_esa$poutlo | out_esa$pouthi
    esa_tails <- esa[esa_tail_flag]
    ESA_oma <- if (length(esa_tails) > 0) mean(abs(esa_tails), na.rm = TRUE) else max(abs(esa), na.rm = TRUE)

    ESA_min <- min(esa, na.rm = TRUE)
    ESA_max <- max(esa, na.rm = TRUE)

    # dESA 
    desa <- .x$dESA
    out_desa <- fndout_R(desa)
    n_desa <- sum(!is.na(desa))

    dESA_low_fence  <- out_desa$fenlo
    dESA_high_fence <- out_desa$fenhi
    dESA_low_out_pct  <- 100 * sum(out_desa$poutlo, na.rm = TRUE) / n_desa
    dESA_high_out_pct <- 100 * sum(out_desa$pouthi, na.rm = TRUE) / n_desa

    desa_tail_flag <- out_desa$poutlo | out_desa$pouthi
    desa_tails <- desa[desa_tail_flag]
    dESA_oma <- if (length(desa_tails) > 0) mean(abs(desa_tails), na.rm = TRUE) else max(abs(desa), na.rm = TRUE)

    dESA_min <- min(desa, na.rm = TRUE)
    dESA_max <- max(desa, na.rm = TRUE)

    tibble(
      n_watersheds = n_distinct(.x$SITENO),
      n_ESA  = n_esa,
      n_dESA = n_desa,

      ESA_median = median(esa, na.rm = TRUE),
      ESA_rstd = mad(esa, na.rm = TRUE),
      ESA_min = ESA_min,
      ESA_max = ESA_max,
      ESA_range = ESA_max - ESA_min,
      ESA_midrange = (ESA_min + ESA_max) / 2,
      ESA_low_fence = ESA_low_fence,
      ESA_high_fence = ESA_high_fence,
      ESA_low_out_pct = ESA_low_out_pct,
      ESA_high_out_pct = ESA_high_out_pct,
      ESA_oma = ESA_oma,     # mean Abs ESA outlier

      dESA_median = 10 * median(desa, na.rm = TRUE),
      dESA_rstd = mad(desa, na.rm = TRUE),
      dESA_min = dESA_min,
      dESA_max = dESA_max,
      dESA_range = dESA_max - dESA_min,
      dESA_low_fence = dESA_low_fence,
      dESA_high_fence = dESA_high_fence,
      dESA_low_out_pct = dESA_low_out_pct,
      dESA_high_out_pct = dESA_high_out_pct,
      dESA_oma  = dESA_oma     # mean Abs dESA/dt outlier
    )
  })

write.csv(partition_stats, here("partition_stats.csv"), row.names = FALSE)

##########################################################
# attribute statistics within each partition, summarise among all watersheds
# min, max, mean area and elevation 
summary_attribute <- partition_attributes %>% 
  group_by(Partition) %>% 
  summarise(min_area = min(area_sqkm), 
            max_area = max(area_sqkm), 
            mean_area = mean(area_sqkm), 
            min_elev = min(elev_mean), 
            max_elev = max(elev_mean), 
            mean_elev = mean(elev_mean))

##############################################################
# IQR of ESA values within each partition across all basins and time points
within_partition_iqr <- esa_data %>%
  group_by(Partition, DATE) %>%
  summarise(iqr_esa = IQR(ESA, na.rm = TRUE), .groups = "drop") %>%
  group_by(Partition) %>%
  summarise(median_iqr = median(iqr_esa, na.rm = TRUE), .groups = "drop")

# median across 10 partitions
overall_median_iqr <- median(within_partition_iqr$median_iqr)


###################################################
# WATER LIMITED VS ENERGY LIMITED TRANSITIONS
##################################################
# proportion of individual observations that are water limited (ESA < 0) vs energy limited (ESA > 0)
transitions <- partition_data %>% group_by(Partition) %>% 
  summarise(pct_above_0 = mean(ESA > 0) * 100, 
            pct_below_0 = mean(ESA < 0) * 100)

# median ESA of each partition using all member watersheds, 88 time steps
median_esa <- partition_data %>%
  group_by(Partition, DATE) %>%
  summarise(median_esa = median(ESA, na.rm = TRUE), .groups = "drop")

# zero crossings 
# partition's median ESA trajectory crosses transition (zero), indicative of regime switching frequency 
count_crossings <- function(x) {sum(diff(sign(x)) != 0, na.rm = TRUE)}

crossings <- median_esa %>%
  group_by(Partition) %>%
  summarise(n_crossings = count_crossings(median_esa), .groups = "drop") %>%
  arrange(Partition)

###########################################
### ESA TIME SERIES PARTITION PLOTS ###
###########################################
# paneled time series figure 

# ensure date format 
partition_data$DATE <- as.Date(partition_data$DATE, format = "%m/%d/%Y")

# get unique watersheds per partition for plot title
n_df <- partition_data %>% 
  group_by(Partition) %>% 
  summarise(n = n_distinct(SITENO)) %>% 
  arrange(Partition)

# creates individual plots per partition
make_partition_plot <- function(p) {
 
  partition_p <- partition_data %>% filter(Partition == p)
  median_p <- median_esa %>% filter(Partition == p)
  n_p <- n_df$n[n_df$Partition == p]
 
  ggplot() +
    geom_line(data = partition_p, aes(x = DATE, y = ESA, group = SITENO),
              color = "darkgray", alpha = 0.25, linewidth = 0.3) +
    geom_line(data = median_p, aes(x = DATE, y = median_esa),
              color = "blue", linewidth = 1) +
    geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.8, color = "red") +
    scale_y_continuous(limits = c(-1, 1), expand = c(0, 0)) +
    scale_x_date(breaks = seq(as.Date("1980-01-01"), as.Date("2025-01-01"), by = "5 years"),
                 labels = scales::date_format("%Y"),
                 limits = as.Date(c("1980-01-01", "2025-01-01"))) +
    # bquote lets p and n_p be inserted dynamically into the plotmath
    # expression -- .( ) substitutes the current variable's value in
    labs(title = bquote(Partition ~ .(p) ~ .(paste0("(n = ", n_p, ")"))),
         y = "ESA", x = NULL) +
    theme_bw() +
    theme(legend.position = "none",
          axis.title = element_text(size = 10),
          plot.title = element_text(size = 10, hjust = 0.5),
          axis.text = element_text(size = 10))
}
 
# create 10 individual plots, in partition order, into one list
partition_plots <- lapply(1:10, make_partition_plot)

# combine into one 2-column figure with shared axes tags (a)
combined_plot <- wrap_plots(partition_plots, ncol = 2) +
  plot_layout(axes = "collect") +
  plot_annotation(tag_levels = "a", tag_prefix = "(", tag_suffix = ") ") &
  theme(plot.tag = element_text(size = 10))
 
ggsave(filename = here("Figures/ESA_partition_timeseries.png"),
       plot = combined_plot, device = "png", width = 10, height = 8, dpi = 300)

###########################################
### dESA TIME SERIES PARTITION PLOTS ###
###########################################
# paneled time series figure for derivatives 

# median dESA of each partition using all member watersheds, 88 time steps
median_desa <- partition_data %>%
  group_by(Partition, DATE) %>%
  summarise(median_desa = median(dESA, na.rm = TRUE), .groups = "drop")

# for shared axes, use same code from ESA time series plot above (make_partition_plot), modify names

# creates individual plots per partition, allowing y-axes values to depend on partition limits of dESA
make_varied_plot <- function(p) {
    partition_p <- partition_data %>% filter(Partition == p)
    median_p <- median_desa %>% filter(Partition == p)
    n_p <- n_df$n[n_df$Partition == p]
  
  # determine y-axis range for this partition
 # y_range <- range(partition_p$dESA, na.rm = TRUE)
  
  # round outward to the nearest integer
 # y_min <- floor(y_range[1])
 # y_max <- ceiling(y_range[2])
  
  # whole number breaks
 # y_breaks <- seq(y_min, y_max, by = 1)
  
  ggplot() +
    geom_line(data = partition_p, aes(x = DATE, y = dESA, group = SITENO),
              color = "darkgray", alpha = 0.25, linewidth = 0.3) +
    geom_line(data = median_p, aes(x = DATE, y = median_desa),
              color = "blue", linewidth = 1) +
    geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.8, color = "red") +
  #  scale_y_continuous(
   #   limits = c(y_min, y_max), breaks = y_breaks, expand = c(0, 0)) + 
    scale_y_continuous(limits = c(-3, 3), expand = c(0, 0)) +
    scale_x_date(
      breaks = seq(as.Date("1980-01-01"), as.Date("2025-01-01"), by = "5 years"),
      labels = scales::date_format("%Y"),
      limits = as.Date(c("1980-01-01", "2025-01-01"))) +
    labs(
      title = bquote(
        Partition ~ .(p) ~ .(paste0("(n = ", n_p, ")"))
      ),
      y = bquote("dESA/dt "(year^-1)), x = NULL) +
    theme_bw() +
    theme(
      legend.position = "none",
      axis.title = element_text(size = 10),
      plot.title = element_text(size = 10, hjust = 0.5),
      axis.text = element_text(size = 10)
    )
}

# create 10 individual plots, in partition order, into one list
partition_plots <- lapply(1:10, make_varied_plot)

# combine into one 2-column figure with shared axes tags (a)
combined_plot <- wrap_plots(partition_plots, ncol = 2) +
 # plot_layout(axes = "keep", axis_titles = "keep") +
  plot_layout(axes = "collect") +
  plot_annotation(tag_levels = "a", tag_prefix = "(", tag_suffix = ") ") &
  theme(plot.tag = element_text(size = 10), 
        axis.text.y = element_text(size = 10), 
        axis.ticks.y = element_line(), 
        axis.line.y = element_line())
 
ggsave(filename = here("Figures/dESA_partition_timeseries.png"),
       plot = combined_plot, device = "png", width = 10, height = 8, dpi = 300)

######################################################
### SPATIAL COHERENCE WITHIN PARTITIONS ###
######################################################
# test whether partition assignment shows significant spatial autocorrelation
# whether basins in the same partition tend to be geographically close, despite
# no geographic information used in clustering process

# site level table with one row per basin, partitions, coordinates
# partition number is spatial variable

# build spatial weights matrix using k-nearest neighbors
coords <- as.matrix(kmeans_r_results[, c("dec_long_va", "dec_lat_va")])
knn <- spdep::knearneigh(coords, k = 8)
nb  <- spdep::knn2nb(knn)
lw  <- spdep::nb2listw(nb, style = "W")

# Moran's I statistical measure of spatial autocorrelation, range -1 to +1 
# values near +1 indicate strong positive spatial autocorrelation (clustering of similar values)
# values near -1 indicate strong negative autocorrelation (dispersion), values near 0 mean random pattern
# Moran's I requires numeric variable, partition number is ordered by aridity
# spatial autocorrelation in partition number is meaningful
moran_result <- spdep::moran.test(kmeans_r_results$Partition, lw, alternative = "greater")

print(moran_result)

###########################################################
### WITHIN PARTITION GEOGRAPHIC CONTIGUITY ###
###########################################################
# quantifies how geographically spread out basins area
# use mean pairwise distance (lat/lon)

# compute mean pairwise great-circle distance (km) within each partition
# smaller mean distance = more geographically contiguous
contiguity_by_partition <- kmeans_r_results %>%
  group_by(Partition) %>%
  group_modify(~{
    if (nrow(.x) < 2) {
      return(tibble(mean_pairwise_km = NA, n_basins = nrow(.x)))
    }
    coords_p <- as.matrix(.x[, c("dec_long_va", "dec_lat_va")])
    dist_mat <- geosphere::distm(coords_p, fun = distHaversine) / 1000  # km
    # exclude diagonal (self-distance = 0)
    dist_vals <- dist_mat[upper.tri(dist_mat)]
    tibble(
      mean_pairwise_km = mean(dist_vals),
      median_pairwise_km = median(dist_vals),
      n_basins = nrow(.x))
  }) %>%
  ungroup() %>%
  arrange(Partition)

# SUPPLEMENTARY TABLE (S2)
print(contiguity_by_partition)

# identify max and min distances between basins in partition (i.e. most contiguous and most dispersed partitions )
min(contiguity_by_partition$mean_pairwise_km)
max(contiguity_by_partition$mean_pairwise_km)

###########################################################
### WITHIN VS BETWEEN PARTITION TIME SERIES SIMILARITY ###
###########################################################
# evaluate the pairwise Pearson correlation between ESA time series within the same ESA partition
# and between randomly sampled basin pairs from different partitions

# build full basin and time series matrix
esa_wide <- esa_df %>% left_join(kmeans_r_results %>% 
  dplyr::select(SITENO, Partition), by = "SITENO") %>% 
  arrange(SITENO, DATE) %>% 
  group_by(SITENO) %>% 
  mutate(t = row_number()) %>% ungroup()

wide_mat <- esa_wide %>%
  select(SITENO, t, ESA) %>%
  pivot_wider(names_from = SITENO, values_from = ESA) %>%
  select(-t) %>%
  as.matrix()

# get partitions for reference 
partition_lookup <- esa_wide %>% distinct(SITENO, Partition)

# correlation within partitions 
full_cor_mat <- cor(wide_mat, use = "pairwise.complete.obs")
basin_ids <- colnames(full_cor_mat)

# extract all within-partition pairwise correlations
within_cors <- numeric(0)
within_partition_id <- character(0)

for (p in unique(partition_lookup$Partition)) {
  basins_p <- partition_lookup$SITENO[partition_lookup$Partition == p]
  if (length(basins_p) < 2) next
  sub_mat <- full_cor_mat[basins_p, basins_p]
  pairwise <- sub_mat[upper.tri(sub_mat)]
  within_cors <- c(within_cors, pairwise)
  within_partition_id <- c(within_partition_id,
                           rep(p, length(pairwise)))
}

# total within partition pairs sampled
length(within_cors)

# per partition breakdown
within_summary_by_partition <- tibble(
  Partition = within_partition_id,
  cor = within_cors) %>%
  mutate(Partition = as.numeric(Partition)) %>%
  group_by(Partition) %>%
  summarise(
    mean_within_cor   = mean(cor, na.rm = TRUE),
    median_within_cor = median(cor, na.rm = TRUE),
    sd_within_cor      = sd(cor, na.rm = TRUE),
    n_pairs            = n(),
    .groups = "drop") %>%
  arrange(Partition)

# SUPPLEMENTARY TABLE (S3)
print(within_summary_by_partition)

write.csv(within_summary_by_partition,
          here("within_partition_similarity.csv"), row.names = FALSE)

####################################################
# extract matched-size random sample of BETWEEN-partition pairwise correlations
# baseline similarity across different partitions

set.seed(123)
n_between_sample <- length(within_cors)  # match sample size

between_cors <- numeric(n_between_sample)

i <- 1
while (i <= n_between_sample) {
  pair <- sample(basin_ids, 2)
  p1 <- partition_lookup$Partition[partition_lookup$SITENO == pair[1]]
  p2 <- partition_lookup$Partition[partition_lookup$SITENO == pair[2]]
  if (p1 != p2) {
    between_cors[i] <- full_cor_mat[pair[1], pair[2]]
    i <- i + 1
  }
}

# total between partition pairs sampled
length(between_cors)

#####################################################
# descriptive summary of WITHIN and BETWEEN partition correlation
cat(sprintf("Within-partition:  mean = %.3f, median = %.3f, SD = %.3f\n",
            mean(within_cors), median(within_cors), sd(within_cors)))
cat(sprintf("Between-partition: mean = %.3f, median = %.3f, SD = %.3f\n",
            mean(between_cors), median(between_cors), sd(between_cors)))
cat(sprintf("Difference in means: %.3f\n",
            mean(within_cors) - mean(between_cors)))

########################################################
# formal significance test
# Wilcoxon rank-sum test (Mann-Whitney U), does not assume normality, 
# appropriate for correlation coefficients which are bounded and often non-normally distributed
# within partition ESA more similar than between 

wilcox_result <- wilcox.test(
  within_cors, between_cors,
  alternative = "greater")  # testing within > between

print(wilcox_result)

# effect size (rank-biserial correlation, a common companion to
# Wilcoxon test for reporting effect magnitude)
# 0.1 small, 0.3 moderate, 0.5+ strong separation
n1 <- length(within_cors)
n2 <- length(between_cors)
r_effect <- 1 - (2 * wilcox_result$statistic) / (n1 * n2)
cat(sprintf("Effect size (rank-biserial r) = %.3f\n", r_effect))

#######################################
# visualization for supplementary figure
# density distribution for pairwise Pearson correlation coefficients between and within partitions

comparison_df <- bind_rows(
  tibble(correlation = within_cors, group = "Within-Partition"),
  tibble(correlation = between_cors, group = "Between-Partition")
)

within_between <- ggplot(comparison_df, aes(x = correlation, fill = group)) +
  geom_density(alpha = 0.6) +
  geom_vline(xintercept = mean(within_cors),
             color = "#1b9e77", linetype = "dashed", linewidth = 0.8) +
  geom_vline(xintercept = mean(between_cors),
             color = "#d95f02", linetype = "dashed", linewidth = 0.8) +
  scale_fill_manual(values = c("Within-Partition" = "#1b9e77",
                               "Between-Partition" = "#d95f02")) +
  labs(
    title = "Within- vs. Between-Partition ESA Time Series Similarity",
    x = "Pairwise Pearson Correlation", y = "Density", fill = NULL) +
  theme_minimal(base_size = 11) +
  theme(legend.position = "bottom")

ggsave(here("Figures/within_between.png"), within_between, dpi = 360, width = 10, height = 8)

###################################################################################
# evaluate relationship between geographic contiguity and within partition temporal correlation
# use mean pairwise distance and within partition correlation 

geo_cor_check <- contiguity_by_partition %>%
  select(Partition, mean_pairwise_km) %>%
  left_join(within_summary_by_partition %>% select(Partition, mean_within_cor), by = "Partition")

cor.test(geo_cor_check$mean_pairwise_km, geo_cor_check$mean_within_cor, method = "pearson")
# checked with pearson and spearmen
# pearson -0.853, rho -0.830