# -----------------------------------------------------------------------------
# Original MSc Cancer Dissertation analysis (2021)
# Refactored for reproducibility and GitHub portfolio (2026)
# -----------------------------------------------------------------------------

# 11_sample_category_fragment_size_summary.R
# Purpose: summarise total LOH fragment size for every
# sample-category combination

library(dplyr)
library(tidyr)
library(readr)

# Load validated LOH dataset containing fragment sizes
loh_segments <- readRDS(
  "data/processed/loh_segments_fragment_sizes.rds"
)

# Confirm required columns exist
required_columns <- c(
  "sample",
  "category",
  "fragment_size_bp"
)

if (!all(required_columns %in% names(loh_segments))) {
  stop(
    "One or more required columns were not found: ",
    paste(required_columns, collapse = ", ")
  )
}

# Calculate total fragment size for observed sample-category combinations
sample_category_fragment_size_summary <- loh_segments %>%
  group_by(sample, category) %>%
  summarise(
    fragment_count = n(),
    total_fragment_size_bp = sum(
      fragment_size_bp,
      na.rm = TRUE
    ),
    .groups = "drop"
  )

# Add missing sample-category combinations and assign zero values
sample_category_fragment_size_complete <-
  sample_category_fragment_size_summary %>%
  complete(
    sample = unique(loh_segments$sample),
    category = sort(unique(loh_segments$category)),
    fill = list(
      fragment_count = 0,
      total_fragment_size_bp = 0
    )
  ) %>%
  arrange(sample, as.integer(category))

# Calculate fragment size in megabases
sample_category_fragment_size_complete <-
  sample_category_fragment_size_complete %>%
  mutate(
    total_fragment_size_mb =
      total_fragment_size_bp / 1e6
  )

# Summarise total fragment size by category
category_total_fragment_size <-
  sample_category_fragment_size_complete %>%
  group_by(category) %>%
  summarise(
    fragment_count = sum(fragment_count),
    total_fragment_size_bp = sum(total_fragment_size_bp),
    total_fragment_size_mb = sum(total_fragment_size_mb),
    .groups = "drop"
  ) %>%
  arrange(as.integer(category))

# Validation output
print(sample_category_fragment_size_complete)
print(category_total_fragment_size, n = Inf)

cat(
  "Number of sample-category combinations:",
  nrow(sample_category_fragment_size_complete),
  "\n"
)

cat(
  "Total fragment count:",
  sum(sample_category_fragment_size_complete$fragment_count),
  "\n"
)

cat(
  "Total fragment size in base pairs:",
  sum(sample_category_fragment_size_complete$total_fragment_size_bp),
  "\n"
)

# Create output directories if needed
dir.create(
  "data/processed",
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  "results",
  recursive = TRUE,
  showWarnings = FALSE
)

# Save processed summary
saveRDS(
  sample_category_fragment_size_complete,
  "data/processed/sample_category_fragment_size_summary.rds"
)

# Export readable result tables
write_csv(
  sample_category_fragment_size_complete,
  "results/sample_category_fragment_size_summary.csv"
)

write_csv(
  category_total_fragment_size,
  "results/category_total_fragment_size.csv"
)

# -----------------------------------------------------------------------------
# Validation note:
# Missing sample-category combinations are retained with fragment counts and
# total fragment sizes set to zero, following the intention of the original
# MSc workflow.
# -----------------------------------------------------------------------------
