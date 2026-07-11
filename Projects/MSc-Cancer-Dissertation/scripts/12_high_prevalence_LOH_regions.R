# -----------------------------------------------------------------------------
# Original MSc Cancer Dissertation analysis (2021)
# Refactored for reproducibility and GitHub portfolio (2026)
# -----------------------------------------------------------------------------

# 12_high_prevalence_LOH_regions.R
# Purpose: identify recurrent (high-prevalence) LOH regions across samples

library(dplyr)
library(readr)

# Load validated LOH dataset
loh_segments <- readRDS(
  "data/processed/loh_segments_fragment_sizes.rds"
)

# Confirm required columns exist
required_columns <- c(
  "chr",
  "startpos",
  "endpos",
  "category",
  "sample"
)

if (!all(required_columns %in% names(loh_segments))) {
  stop(
    "One or more required columns were not found: ",
    paste(required_columns, collapse = ", ")
  )
}

# Count how many tumour samples contain each identical LOH region
high_prevalence_regions <- loh_segments %>%
  group_by(
    chr,
    category,
    startpos,
    endpos
  ) %>%
  summarise(
    sample_count = n_distinct(sample),
    fragment_count = n(),
    total_fragment_size_bp = first(fragment_size_bp),
    .groups = "drop"
  ) %>%
  arrange(
    desc(sample_count),
    chr,
    startpos
  )

# Validation summaries
cat(
  "Unique genomic regions:",
  nrow(high_prevalence_regions),
  "\n"
)

cat(
  "Maximum sample prevalence:",
  max(high_prevalence_regions$sample_count),
  "\n"
)

cat(
  "Regions present in two or more samples:",
  sum(high_prevalence_regions$sample_count >= 2),
  "\n"
)

print(high_prevalence_regions, n = 20)

# Create output directories
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

# Save processed object
saveRDS(
  high_prevalence_regions,
  "data/processed/high_prevalence_LOH_regions.rds"
)

# Export readable table
write_csv(
  high_prevalence_regions,
  "results/high_prevalence_LOH_regions.csv"
)

# -----------------------------------------------------------------------------
# Validation note:
# LOH regions are grouped by chromosome, category, genomic start position and
# genomic end position. Each region is summarised by the number of distinct
# tumour-region samples containing that genomic interval.
# -----------------------------------------------------------------------------
