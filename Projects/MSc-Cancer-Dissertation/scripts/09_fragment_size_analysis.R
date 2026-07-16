# -----------------------------------------------------------------------------
# Original MSc Cancer Dissertation analysis (2021)
# Refactored for reproducibility and GitHub portfolio (2026)
# -----------------------------------------------------------------------------

# 09_fragment_size_analysis.R
# Purpose: calculate and summarise LOH fragment sizes

library(dplyr)
library(readr)

# Load validated LOH dataset from Script 4
loh_segments <- readRDS(
  "data/processed/loh_segments_categories.rds"
)

# Confirm required columns exist
required_columns <- c(
  "startpos",
  "endpos",
  "category",
  "clonality"
)

if (!all(required_columns %in% names(loh_segments))) {
  stop(
    "One or more required columns were not found: ",
    paste(required_columns, collapse = ", ")
  )
}

# Calculate fragment size using the original MSc method
loh_segments <- loh_segments %>%
  mutate(
    fragment_size_bp = endpos - startpos,
    fragment_size_mb = fragment_size_bp / 1e6
  )

# Check for invalid or unusual fragment sizes
negative_fragment_count <- sum(
  loh_segments$fragment_size_bp < 0,
  na.rm = TRUE
)

zero_length_fragment_count <- sum(
  loh_segments$fragment_size_bp == 0,
  na.rm = TRUE
)

if (negative_fragment_count > 0) {
  warning(
    negative_fragment_count,
    " fragment(s) have an end position smaller than the start position."
  )
}

# Overall fragment-size summary
fragment_size_summary <- loh_segments %>%
  summarise(
    number_of_fragments = n(),
    minimum_size_bp = min(fragment_size_bp, na.rm = TRUE),
    first_quartile_size_bp = quantile(
      fragment_size_bp,
      0.25,
      na.rm = TRUE
    ),
    median_size_bp = median(fragment_size_bp, na.rm = TRUE),
    mean_size_bp = mean(fragment_size_bp, na.rm = TRUE),
    third_quartile_size_bp = quantile(
      fragment_size_bp,
      0.75,
      na.rm = TRUE
    ),
    maximum_size_bp = max(fragment_size_bp, na.rm = TRUE),
    zero_length_fragments = zero_length_fragment_count
  )

# Summarise fragment sizes by LOH category
category_fragment_size_summary <- loh_segments %>%
  group_by(category) %>%
  summarise(
    number_of_fragments = n(),
    total_size_bp = sum(fragment_size_bp, na.rm = TRUE),
    mean_size_bp = mean(fragment_size_bp, na.rm = TRUE),
    median_size_bp = median(fragment_size_bp, na.rm = TRUE),
    minimum_size_bp = min(fragment_size_bp, na.rm = TRUE),
    maximum_size_bp = max(fragment_size_bp, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(as.integer(category))

# Summarise fragment sizes by category and clonality
category_clonality_fragment_size_summary <- loh_segments %>%
  group_by(category, clonality) %>%
  summarise(
    number_of_fragments = n(),
    total_size_bp = sum(fragment_size_bp, na.rm = TRUE),
    mean_size_bp = mean(fragment_size_bp, na.rm = TRUE),
    median_size_bp = median(fragment_size_bp, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(as.integer(category), clonality)

# Display validation summaries
print(fragment_size_summary)
print(category_fragment_size_summary, n = Inf)
print(category_clonality_fragment_size_summary, n = Inf)

cat(
  "Negative fragment sizes:",
  negative_fragment_count,
  "\n"
)

cat(
  "Zero-length fragments:",
  zero_length_fragment_count,
  "\n"
)

# Create output directories if needed
dir.create("data/processed", recursive = TRUE, showWarnings = FALSE)
dir.create("results", recursive = TRUE, showWarnings = FALSE)

# Save dataset containing calculated fragment sizes
saveRDS(
  loh_segments,
  "data/processed/loh_segments_fragment_sizes.rds"
)

# Save summary tables
saveRDS(
  fragment_size_summary,
  "data/processed/fragment_size_summary.rds"
)

write_csv(
  fragment_size_summary,
  "results/fragment_size_summary.csv"
)

write_csv(
  category_fragment_size_summary,
  "results/category_fragment_size_summary.csv"
)

write_csv(
  category_clonality_fragment_size_summary,
  "results/category_clonality_fragment_size_summary.csv"
)

# -----------------------------------------------------------------------------
# Validation note:
# Fragment size is calculated as endpos - startpos to preserve the original
# MSc method. Consequently, single-position intervals where startpos equals
# endpos have a calculated fragment size of zero.
# -----------------------------------------------------------------------------