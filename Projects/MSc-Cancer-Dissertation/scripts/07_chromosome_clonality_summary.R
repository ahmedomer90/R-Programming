# 07_chromosome_clonality_summary.R
# MSc Cancer Dissertation reconstruction
# Purpose: summarise clonality patterns by chromosome

# Load packages
library(dplyr)
library(readr)

# Load LOH segments with category and clonality information
loh_segments <- readRDS("data/processed/loh_segments_categories.rds")

# Count LOH fragments by chromosome and clonality
chromosome_clonality_counts <- loh_segments %>%
  count(chr, clonality, name = "count")

# Calculate proportions within each chromosome
chromosome_clonality_summary <- chromosome_clonality_counts %>%
  group_by(chr) %>%
  mutate(
    proportion = count / sum(count)
  ) %>%
  ungroup()

# View summary
print(chromosome_clonality_summary)

# Save outputs
saveRDS(
  chromosome_clonality_summary,
  "data/processed/chromosome_clonality_summary.rds"
)

write_csv(
  chromosome_clonality_summary,
  "results/chromosome_clonality_summary.csv"
)
