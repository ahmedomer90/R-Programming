# 05_category_clonality_summary.R
# MSc Cancer Dissertation reconstruction
# Purpose: summarise LOH categories by clonality

# Load packages
library(dplyr)
library(readr)

# Load LOH segments with category and clonality information
loh_segments <- readRDS("data/processed/loh_segments_categories.rds")

# Count LOH fragments by category and clonality
category_clonality_counts <- loh_segments %>%
  count(LOH_category, clonality, name = "count")

# Calculate proportions within each LOH category
category_clonality_summary <- category_clonality_counts %>%
  group_by(LOH_category) %>%
  mutate(
    proportion = count / sum(count)
  ) %>%
  ungroup()

# View summary
print(category_clonality_summary)

# Save outputs
saveRDS(
  category_clonality_summary,
  "data/processed/category_clonality_summary.rds"
)

write_csv(
  category_clonality_summary,
  "results/category_clonality_summary.csv"
)
