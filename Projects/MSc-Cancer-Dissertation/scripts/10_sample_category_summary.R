# -----------------------------------------------------------------------------
# Original MSc Cancer Dissertation analysis (2021)
# Refactored for reproducibility and GitHub portfolio (2026)
# -----------------------------------------------------------------------------

# 10_sample_category_summary.R
# Purpose: summarise the number of LOH fragments in each category
# for every tumour-region sample

library(dplyr)
library(tidyr)
library(readr)

# Load validated LOH dataset
loh_segments <- readRDS(
  "data/processed/loh_segments_fragment_sizes.rds"
)

# Confirm required columns exist
required_columns <- c("sample", "category")

if (!all(required_columns %in% names(loh_segments))) {
  stop(
    "Required columns 'sample' and/or 'category' were not found."
  )
}

# Count fragments by sample and category
sample_category_summary <- loh_segments %>%
  count(sample, category, name = "count")

# Create all sample/category combinations
all_combinations <- expand.grid(
  sample = unique(loh_segments$sample),
  category = sort(unique(loh_segments$category)),
  stringsAsFactors = FALSE
)

# Include combinations with zero fragments
sample_category_complete <- all_combinations %>%
  left_join(
    sample_category_summary,
    by = c("sample", "category")
  ) %>%
  mutate(
    count = replace_na(count, 0)
  ) %>%
  arrange(sample, as.integer(category)) %>%
  tibble::as_tibble()

# Validation
print(sample_category_complete, n = 20)

category_totals <- sample_category_complete %>%
  group_by(category) %>%
  summarise(
    total_fragments = sum(count),
    .groups = "drop"
  )

print(category_totals)

# Create output directories
dir.create("data/processed", recursive = TRUE, showWarnings = FALSE)
dir.create("results", recursive = TRUE, showWarnings = FALSE)

# Save outputs
saveRDS(
  sample_category_complete,
  "data/processed/sample_category_summary.rds"
)

write_csv(
  sample_category_complete,
  "results/sample_category_summary.csv"
)