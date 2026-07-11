# -----------------------------------------------------------------------------
# Original MSc Cancer Dissertation analysis (2021)
# Refactored for reproducibility and GitHub portfolio (2026)
# -----------------------------------------------------------------------------

# 05_category_clonality_summary.R
# Purpose: summarise the count and proportion of clonality states
# within each LOH category

library(dplyr)
library(readr)

# Load validated output from Script 4
loh_segments <- readRDS(
  "data/processed/loh_segments_categories.rds"
)

# Confirm required columns exist
required_columns <- c("category", "clonality")

if (!all(required_columns %in% names(loh_segments))) {
  stop("Required columns 'category' and/or 'clonality' were not found.")
}

# Count fragments by LOH category and clonality
category_clonality_summary <- loh_segments %>%
  count(category, clonality, name = "count") %>%
  group_by(category) %>%
  mutate(
    proportion = count / sum(count)
  ) %>%
  ungroup() %>%
  arrange(as.numeric(category), clonality)

# Display results
print(category_clonality_summary)

# Check that proportions sum to 1 within each category
proportion_check <- category_clonality_summary %>%
  group_by(category) %>%
  summarise(
    total_proportion = sum(proportion),
    .groups = "drop"
  )

print(proportion_check)

# Create output directories if needed
dir.create("data/processed", recursive = TRUE, showWarnings = FALSE)
dir.create("results", recursive = TRUE, showWarnings = FALSE)

# Save outputs
saveRDS(
  category_clonality_summary,
  "data/processed/category_clonality_summary.rds"
)

write_csv(
  category_clonality_summary,
  "results/category_clonality_summary.csv"
)
