# 10_sample_category_summary.R
# MSc Cancer Dissertation reconstruction
# Purpose: summarise LOH categories by sample

# Load packages
library(dplyr)
library(readr)

# Load LOH segments with category information
loh_segments <- readRDS("data/processed/loh_segments_categories.rds")

# Count LOH categories in each sample
sample_category_summary <- loh_segments %>%
  count(sample, LOH_category, name = "count") %>%
  group_by(sample) %>%
  mutate(
    proportion = count / sum(count)
  ) %>%
  ungroup()

# View summary
print(sample_category_summary)

# Save outputs
saveRDS(
  sample_category_summary,
  "data/processed/sample_category_summary.rds"
)

write_csv(
  sample_category_summary,
  "results/sample_category_summary.csv"
)
