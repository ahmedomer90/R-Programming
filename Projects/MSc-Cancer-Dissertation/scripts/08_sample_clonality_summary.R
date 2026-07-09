# 08_sample_clonality_summary.R
# MSc Cancer Dissertation reconstruction
# Purpose: summarise clonality patterns by sample

# Load packages
library(dplyr)
library(readr)

# Load LOH segments with category and clonality information
loh_segments <- readRDS("data/processed/loh_segments_categories.rds")

# Count LOH fragments by sample and clonality
sample_clonality_counts <- loh_segments %>%
  count(sample, clonality, name = "count")

# Calculate proportions within each sample
sample_clonality_summary <- sample_clonality_counts %>%
  group_by(sample) %>%
  mutate(
    proportion = count / sum(count)
  ) %>%
  ungroup()

# View summary
print(sample_clonality_summary)

# Save outputs
saveRDS(
  sample_clonality_summary,
  "data/processed/sample_clonality_summary.rds"
)

write_csv(
  sample_clonality_summary,
  "results/sample_clonality_summary.csv"
)
