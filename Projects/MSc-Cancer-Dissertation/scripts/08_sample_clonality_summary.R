# -----------------------------------------------------------------------------
# Original MSc Cancer Dissertation analysis (2021)
# Refactored for reproducibility and GitHub portfolio (2026)
# -----------------------------------------------------------------------------

# 08_sample_clonality_summary.R
# Purpose: summarise LOH clonality counts and proportions by sample

library(dplyr)
library(readr)

# Load validated LOH dataset from Script 4
loh_segments <- readRDS(
  "data/processed/loh_segments_categories.rds"
)

# Confirm required columns exist
required_columns <- c("sample", "clonality")

if (!all(required_columns %in% names(loh_segments))) {
  stop(
    "Required columns 'sample' and/or 'clonality' were not found."
  )
}

# Count LOH fragments by sample and clonality
sample_clonality_summary <- loh_segments %>%
  count(sample, clonality, name = "count") %>%
  group_by(sample) %>%
  mutate(
    proportion = count / sum(count)
  ) %>%
  ungroup() %>%
  arrange(sample, clonality)

# Check that proportions sum to 1 within every sample
sample_proportion_check <- sample_clonality_summary %>%
  group_by(sample) %>%
  summarise(
    total_proportion = sum(proportion),
    .groups = "drop"
  )

# Produce separate clonal and subclonal datasets
# These follow the data-splitting approach described in the original analysis.
clonal_loh_segments <- loh_segments %>%
  filter(clonality == "clonal")

subclonal_loh_segments <- loh_segments %>%
  filter(clonality == "subclonal")

# Retain and display the rare "none" cases transparently
none_loh_segments <- loh_segments %>%
  filter(clonality == "none")

# Display summaries
print(sample_clonality_summary, n = Inf)
print(sample_proportion_check, n = Inf)

table(loh_segments$clonality, useNA = "ifany")

# Create output directories if needed
dir.create("data/processed", recursive = TRUE, showWarnings = FALSE)
dir.create("results", recursive = TRUE, showWarnings = FALSE)

# Save processed objects
saveRDS(
  sample_clonality_summary,
  "data/processed/sample_clonality_summary.rds"
)

saveRDS(
  clonal_loh_segments,
  "data/processed/clonal_loh_segments.rds"
)

saveRDS(
  subclonal_loh_segments,
  "data/processed/subclonal_loh_segments.rds"
)

saveRDS(
  none_loh_segments,
  "data/processed/none_loh_segments.rds"
)

# Export readable summary table
write_csv(
  sample_clonality_summary,
  "results/sample_clonality_summary.csv"
)

# -----------------------------------------------------------------------------
# Validation note:
# The 10 fragments classified as "none" during Script 3 are retained in the
# sample summary to preserve the complete validated LOH dataset. Separate
# clonal and subclonal datasets exclude these cases by definition.
# -----------------------------------------------------------------------------