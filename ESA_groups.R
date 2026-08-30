###########################################
### GROUPING PARTITIONS ###
###########################################
# Katharine Sink

# group ESA partitions (n = 10) based on summary statistics for ESA and dESA
# second-stage k-means classification applied at the PARTITION level (not the watershed level)
# each of the 10 ESA partitions is treated as a single observation, described by a feature vector of
# pooled distributional statistics computed across all watersheds and all 88 time steps within that partition
# (i.e. the columns already produced) in `partition_stats` by the earlier "PARTITION LEVEL ESA STATISTICS"
# section of the ESA_partitions.R
# two independent 5-group classifications are produced:
#   ESArstd groups:  based on the central distribution of ESA values
#   ESAtails groups: based on the extreme-value (dESA/dt) structure

# requires "partition_stats" to already exist (from ESA_partitions.R) with one row per ESA partition (n = 10) 
# and the relevant columns

library(tidyverse)
library(ClusterR)

##################################
## KMEANS ON PARTITION STATS ###
##################################
partition_level_kmeans <- function(feature_df, order_cols, k = 5, num_init = 50, seed = 1) {
 
  X <- as.matrix(feature_df)   # raw units, no scaling (standardization or z score)
 
  km <- KMeans_rcpp(
    data = X,
    clusters = k,
    num_init = num_init,
    max_iters = 100,
    initializer = "kmeans++",
    seed = seed
  )
 
  grp <- km$clusters
 
  # order groups by decreasing mean of (possibly summed) raw ordering variable(s)
  raw_order_val <- rowSums(feature_df[, order_cols, drop = FALSE])
  group_means <- tapply(raw_order_val, grp, mean)
  ord <- order(group_means, decreasing = TRUE)
 
  old_grp <- grp
  for (i in seq_along(ord)) {
    grp[old_grp == ord[i]] <- i
  }
 
  setNames(grp, rownames(feature_df))
}

############################################################
# classification of partitions based on central variability 
# ESA groups (ESArstd)
# Matlab code = ESAgroup.m 

# select statistics for central variability 
ESArstd_features <- partition_stats %>%
  dplyr::select(Partition, ESA_median, ESA_rstd, ESA_midrange, ESA_range,
                ESA_min, ESA_max, ESA_low_fence, ESA_high_fence,
                ESA_low_out_pct, ESA_high_out_pct) %>%
  tibble::column_to_rownames("Partition")
 
# kmeans central variability groups
ESArstd_group <- partition_level_kmeans(
  ESArstd_features, order_cols = "ESA_rstd", k = 5, num_init = 50, seed = 1)

###############################################################
# classification of partitions based on derivative tail behavior
# Tails groups (ESAtails)
# Matlab code = dESAdtgroupt.m

# select statistics for volatility 
ESAtails_features <- partition_stats %>%
  dplyr::select(Partition, dESA_min, dESA_max, dESA_low_fence, dESA_high_fence,
                dESA_low_out_pct, dESA_high_out_pct) %>%
  tibble::column_to_rownames("Partition")
 
# kmeans volatility groups 
ESAtails_group <- partition_level_kmeans(
  ESAtails_features, order_cols = c("dESA_low_out_pct", "dESA_high_out_pct"), 
  k = 5, num_init = 50, seed = 1)

#######################################
# summarize groups

# add group numbers to partition stats table, as character
partition_groups <- partition_stats %>% 
  dplyr::mutate(
    dESA_tail = dESA_low_out_pct + dESA_high_out_pct, 
    ESArstd_group = ESArstd_group[as.character(Partition)], 
    ESAtails_group = ESAtails_group[as.character(Partition)])

write.csv(partition_groups, here("partition_groups.csv"), row.names = FALSE)

# TABLE 1 (MANUSCRIPT)
# order rows by descending rstd value 
table1 <- partition_groups %>% 
  dplyr::transmute(
    Label = paste0(Partition, "/", ESArstd_group, "/", ESAtails_group), 
    `ESA rstd` = round(ESA_rstd, 2), 
    `ESA median` = round(ESA_median, 2),
    `dESA/dt rstd` = round(dESA_rstd, 2),
    `dESA/dt tail %` = round(dESA_tail, 2),
    `Mean |dESA/dt| outlier` = round(dESA_oma, 2)
  ) %>%
  dplyr::arrange(dplyr::desc(`ESA rstd`))  


print(table1)

#################################
## HEATMAP FOR GROUPS ##
#################################
# visual display matrix for partitions and group assignments 
# using designated group colors

group_table <- partition_groups %>%
  dplyr::select(Partition, ESA_median, ESArstd_group, ESAtails_group) %>%
  arrange(Partition) %>%
  mutate(Partition = factor(Partition, levels = 1:10))

group_long <- group_table %>%
  pivot_longer(
    cols = c(ESArstd_group, ESAtails_group),
    names_to = "Classification",
    values_to = "Group") %>%
  mutate(Group = factor(Group, levels = 1:5))

group_colors <- c(
  "1" = "#D55E00",
  "2" = "#F0E442",
  "3" = "#009E73",
  "4" = "#0072B2",
  "5" = "#CC79A7"
)

grps_plot <- ggplot(group_long,
       aes(x = Classification, y = Partition, fill = Group)) +
  geom_tile(color = "white", linewidth = 1.2) +
  geom_text(aes(label = as.character(Group)),
            color = "white", fontface = "bold", size = 5) +
  scale_fill_manual(values = group_colors, name = "Group") +
  scale_x_discrete(
    labels = c(
      "ESArstd_group" = expression(ESA[rstd]),
      "ESAtails_group"= expression(ESA[tails])
    )
  ) +
  scale_y_discrete(limits = rev) +
  labs(
    x = NULL,
    y = "ESA Partition\n(1 = most arid, 10 = most humid)"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(size = 11, face = "bold"),
    axis.text.y = element_text(size = 10),
    plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
    legend.position = "right"
  )
 

ggsave(filename = here("Figures/grps_plot.png"), 
       plot = grps_plot, device = "png", width = 10, height = 8, dpi = 300)

#################################
## CORRELATION BETWEEN GROUPS ##
#################################
# relationship between variability (ESArstd) and volatility (dESAtail %)
# quantify correlation between variability and volatility
cor_test <- cor.test(
  partition_groups$ESA_rstd, 
  partition_groups$dESA_tail, 
  method = "spearman")

cat("Spearman r:", cor_test$estimate)
cat("p-value:", cor_test$p.value)

############################################
# spearman bootstrap confidence interval
# bootstrap resampling of 10 partitions to estimate stability of ESArstd vs tail
# frequency trade off

set.seed(123)
n_boot <- 10000

# one row per partition
boot_data <- partition_groups %>%
  select(Partition, ESA_rstd, dESA_tail) %>%
  drop_na()

n_partitions <- nrow(boot_data)

boot_rhos <- numeric(n_boot)
for (i in 1:n_boot) {
  boot_idx <- sample(1:n_partitions, n_partitions, replace = TRUE)
  boot_sample <- boot_data[boot_idx, ]

  # skip degenerate resamples with no variance (rare but possible
  # with small n and replacement sampling)
  if (length(unique(boot_sample$ESA_rstd)) < 3 ||
      length(unique(boot_sample$dESA_tail)) < 3) {
    boot_rhos[i] <- NA
    next
  }

  boot_rhos[i] <- cor(
    boot_sample$ESA_rstd,
    boot_sample$dESA_tail,
    method = "spearman"
  )
}

boot_rhos <- boot_rhos[!is.na(boot_rhos)]

ci_lower <- quantile(boot_rhos, 0.025)
ci_upper <- quantile(boot_rhos, 0.975)
boot_mean <- mean(boot_rhos)
boot_se   <- sd(boot_rhos)

# SUMMARY
cat(sprintf("Original Spearman rho = %.3f\n", cor(boot_data$ESA_rstd, boot_data$dESA_tail, method = "spearman")))
cat(sprintf("Bootstrap mean rho = %.3f (SE = %.3f)\n", boot_mean, boot_se))
cat(sprintf("95%% Bootstrap CI: [%.3f, %.3f]\n", ci_lower, ci_upper))
cat(sprintf("Valid bootstrap replicates: %d / %d\n", length(boot_rhos), n_boot))

# SUPPLEMENTARY FIGURE (S3)
# histogram of bootstrap distribution 
rstd_tails_cor <- ggplot(data.frame(rho = boot_rhos), aes(x = rho)) +
  geom_histogram(bins = 50, fill = "steelblue", color = "white") +
  geom_vline(xintercept = ci_lower, linetype = "dashed", color = "red") +
  geom_vline(xintercept = ci_upper, linetype = "dashed", color = "red") +
  geom_vline(xintercept = cor(boot_data$ESA_rstd, boot_data$dESA_tail,
                              method = "spearman"), color = "black", linewidth = 0.8) +
  labs(title = "Bootstrap Distribution of ESArstd–ESAtails Spearman Correlation",
    x = "Spearman \u03C1 (resampled)", y = "Count") +
  theme_minimal(base_size = 11)

ggsave(here("Figures/rstd_tails_cor.png"), rstd_tails_cor, dpi = 360, width = 10, height = 8)


##################################
## ESA DISTRIBUTION PLOTS ###
##################################
# one probability density histogram per partition, 
# title with partition number and ESArstd group assignment 
# creates individual plots per partition
# requires the partition_data from ESA_partitions.R

make_esapdf_plot <- function(p) {
 
  partition_p <- partition_data %>% filter(Partition == p)

  # partition ESArstd group number, from partition_groups
  group_p <- partition_groups$ESArstd_group[partition_groups$Partition == p]
  
  ggplot(data = partition_p, aes(x = ESA)) +
    geom_histogram(aes(y = after_stat(count/sum(count))), bins = 40, 
              fill = "steelblue", color = "black") + 
    labs(x = "ESA", y = "PDF", 
      title = bquote(atop(.(paste0("Partition ", p)), 
                     ESA[rstd] ~ group ~ .(group_p)))) + 
    coord_cartesian(xlim = c(-1, 1)) + 
    scale_y_continuous(expand = c(0, 0)) + 
    theme_bw() +
    theme(plot.title = element_text(size = 10, hjust = 0.5),
          axis.title = element_text(size = 10),
          axis.text = element_text(size = 10))
}
 
# sort partitions by ESArstd group (ascending), then by partition number within each group
plot_order <- partition_groups %>% 
  arrange(ESArstd_group, Partition) %>% 
  pull(Partition)

# create 10 individual plots, in partition order, into one list
esa_pdf_plots <- lapply(plot_order, make_esapdf_plot)

# combine into one 2-column figure with shared axes tags (a)
combined_esapdf_plot <- wrap_plots(esa_pdf_plots, ncol = 4) +
  plot_annotation(tag_levels = "a", tag_prefix = "(", tag_suffix = ") ") &
  theme(plot.tag = element_text(size = 10))
 
ggsave(filename = here("Figures/pdf_esa_plot.png"),
       plot = combined_esapdf_plot, device = "png", width = 10, height = 8, dpi = 300)

#####################################
## DERIVATIVE DISTRIBUTION PLOTS ###
#####################################
# one probability density histogram per partition for the dESA/dt derivative values
# with low/high Tukey-hinge fence lines, groups based on ESAtails 


make_dESApdf_plot <- function(p) {
 
  partition_p <- partition_data %>% filter(Partition == p)
 
  # ESAtails group number, looked up from partition_groups
  group_p <- partition_groups$ESAtails_group[partition_groups$Partition == p]
 
  # low/high dESA fences 
  fence_lo <- partition_groups$dESA_low_fence[partition_groups$Partition == p]
  fence_hi <- partition_groups$dESA_high_fence[partition_groups$Partition == p]
 
  ggplot(data = partition_p, aes(x = dESA)) +
    geom_histogram(
      aes(y = after_stat(count / sum(count))),
      bins = 40,
      fill = "steelblue",
      color = "black"
    ) +
    geom_vline(xintercept = fence_lo, color = "firebrick", linetype = "dashed", linewidth = 0.6) +
    geom_vline(xintercept = fence_hi, color = "firebrick", linetype = "dashed", linewidth = 0.6) +
    labs(
      x = "dESA/dt",
      y = "PDF",
      title = bquote(atop(.(paste0("Partition ", p)),
                           ESA[tails] ~ group ~ .(group_p)))
    ) +
    coord_cartesian(xlim = c(-3, 3)) +
    scale_y_continuous(expand = c(0, 0)) +
    theme_bw() +
    theme(
      plot.title = element_text(size = 10, hjust = 0.5),
      axis.title = element_text(size = 10),
      axis.text  = element_text(size = 10)
    )
}
 
# sort partitions by ESAtails_group (ascending), then by partition number
plot_order <- partition_groups %>%
  arrange(ESAtails_group, Partition) %>%
  pull(Partition)
 
dESA_pdf_plots <- lapply(plot_order, make_dESApdf_plot)
 
dESA_pdf_combined <- wrap_plots(dESA_pdf_plots, ncol = 4) +
  plot_annotation(tag_levels = "a", tag_prefix = "(", tag_suffix = ") ") &
  theme(plot.tag = element_text(size = 10))
 
ggsave(filename = here("Figures/pdf_dESA_plots.png"),
       plot = dESA_pdf_combined, device = "png", width = 10, height = 8, dpi = 300)