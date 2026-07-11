# -----------------------------------------------------------------------------
# Original MSc Cancer Dissertation analysis (2021)
# Refactored for reproducibility and GitHub portfolio (2026)
# -----------------------------------------------------------------------------

# 07_chromosome_clonality_summary.R
# Purpose: summarise clonality counts and proportions by chromosome
# for both all ASCAT fragments and LOH fragments

library(dplyr)
library(readr)

# Load validated outputs from Scripts 3 and 4
tracerx_ascat_seg <- readRDS(
  "data/processed/tracerx_ascat_seg_categories.rds"
)

loh_segments <- readRDS(
  "data/processed/loh_segments_categories.rds"
)

# Confirm required columns exist
required_columns <- c("chr", "clonality")

if (!all(required_columns %in% names(tracerx_ascat_seg))) {
  stop(
    "Required columns 'chr' and/or 'clonality' were not found ",
    "in the full ASCAT dataset."
  )
}

if (!all(required_columns %in% names(loh_segments))) {
  stop(
    "Required columns 'chr' and/or 'clonality' were not found ",
    "in the LOH dataset."
  )
}

# -------------------------------------------------------------------------
# 1. All ASCAT fragments
# -------------------------------------------------------------------------

all_chromosome_clonality_summary <- tracerx_ascat_seg %>%
  count(chr, clonality, name = "count") %>%
  group_by(chr) %>%
  mutate(
    proportion = count / sum(count)
  ) %>%
  ungroup() %>%
  arrange(chr, clonality)

# Check proportions for all fragments
all_proportion_check <- all_chromosome_clonality_summary %>%
  group_by(chr) %>%
  summarise(
    total_proportion = sum(proportion),
    .groups = "drop"
  )

# -------------------------------------------------------------------------
# 2. LOH fragments only
# -------------------------------------------------------------------------

loh_chromosome_clonality_summary <- loh_segments %>%
  count(chr, clonality, name = "count") %>%
  group_by(chr) %>%
  mutate(
    proportion = count / sum(count)
  ) %>%
  ungroup() %>%
  arrange(chr, clonality)

# Check proportions for LOH fragments
loh_proportion_check <- loh_chromosome_clonality_summary %>%
  group_by(chr) %>%
  summarise(
    total_proportion = sum(proportion),
    .groups = "drop"
  )

# Display summaries and validation checks
print(all_chromosome_clonality_summary, n = Inf)
print(all_proportion_check, n = Inf)

print(loh_chromosome_clonality_summary, n = Inf)
print(loh_proportion_check, n = Inf)

# Create output directories if needed
dir.create("data/processed", recursive = TRUE, showWarnings = FALSE)
dir.create("results", recursive = TRUE, showWarnings = FALSE)

# Save processed summaries
saveRDS(
  all_chromosome_clonality_summary,
  "data/processed/all_chromosome_clonality_summary.rds"
)

saveRDS(
  loh_chromosome_clonality_summary,
  "data/processed/loh_chromosome_clonality_summary.rds"
)

# Export readable result tables
write_csv(
  all_chromosome_clonality_summary,
  "results/all_chromosome_clonality_summary.csv"
)

write_csv(
  loh_chromosome_clonality_summary,
  "results/loh_chromosome_clonality_summary.csv"
)
