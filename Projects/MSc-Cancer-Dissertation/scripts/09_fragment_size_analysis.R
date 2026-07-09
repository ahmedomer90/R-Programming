# 09_fragment_size_analysis.R
# MSc Cancer Dissertation reconstruction
# Purpose: calculate and summarise LOH fragment sizes

# Load packages
library(dplyr)
library(readr)

# Load LOH segments
loh_segments <- readRDS("data/processed/loh_segments_categories.rds")

# Calculate fragment size (base pairs)
loh_segments <- loh_segments %>%
  mutate(
    fragment_size_bp = endpos - startpos,
    fragment_size_mb = fragment_size_bp / 1e6
  )

# Summary statistics
fragment_size_summary <- loh_segments %>%
  summarise(
    Number_of_fragments = n(),
    Mean_size_mb = mean(fragment_size_mb, na.rm = TRUE),
    Median_size_mb = median(fragment_size_mb, na.rm = TRUE),
    Minimum_size_mb = min(fragment_size_mb, na.rm = TRUE),
    Maximum_size_mb = max(fragment_size_mb, na.rm = TRUE)
  )

# Display summary
print(fragment_size_summary)

# Save outputs
saveRDS(
  loh_segments,
  "data/processed/loh_segments_fragment_sizes.rds"
)

write_csv(
  fragment_size_summary,
  "results/fragment_size_summary.csv"
)
