# 11_high_prevalence_LOH_regions.R
# MSc Cancer Dissertation reconstruction
# Purpose: identify chromosomal regions with high prevalence of LOH categories

# Load packages
library(dplyr)
library(readr)

# Load LOH segments with category information
loh_segments <- readRDS("data/processed/loh_segments_categories.rds")

# Count LOH events by chromosome, start position, end position, and category
loh_region_counts <- loh_segments %>%
  count(chr, startpos, endpos, LOH_category, name = "LOH_count")

# Calculate thresholds for high-prevalence regions
median_threshold <- median(loh_region_counts$LOH_count, na.rm = TRUE)
upper_quartile_threshold <- quantile(loh_region_counts$LOH_count, 0.75, na.rm = TRUE)

# Identify regions above median threshold
regions_above_median <- loh_region_counts %>%
  filter(LOH_count > median_threshold)

# Identify regions above upper-quartile threshold
regions_above_upper_quartile <- loh_region_counts %>%
  filter(LOH_count > upper_quartile_threshold)

# Save outputs
write_csv(
  loh_region_counts,
  "results/loh_region_counts.csv"
)

write_csv(
  regions_above_median,
  "results/loh_regions_above_median.csv"
)

write_csv(
  regions_above_upper_quartile,
  "results/loh_regions_above_upper_quartile.csv"
)

saveRDS(
  regions_above_upper_quartile,
  "data/processed/high_prevalence_LOH_regions.rds"
)
