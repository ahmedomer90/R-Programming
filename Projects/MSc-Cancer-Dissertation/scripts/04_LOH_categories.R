# -----------------------------------------------------------------------------
# Original MSc Cancer Dissertation analysis (2021)
# Refactored for reproducibility and GitHub portfolio (2026)
# -----------------------------------------------------------------------------

# 04_LOH_categories.R
# Purpose: classify ASCAT fragments into eight LOH categories based on nMajor

library(dplyr)

# Load the full ASCAT dataset with LOH and clonality classifications
tracerx_ascat_seg <- readRDS(
  "data/processed/tracerx_ascat_seg_clonality.rds"
)

# Confirm required columns exist
required_columns <- c("nMajor", "LOH")

if (!all(required_columns %in% names(tracerx_ascat_seg))) {
  stop("Required columns 'nMajor' and/or 'LOH' were not found.")
}

# Create LOH categories using the original MSc logic
tracerx_ascat_seg <- tracerx_ascat_seg %>%
  mutate(
    category = nMajor + 1,
    category = if_else(category >= 8, 8, category),
    category = if_else(LOH, as.character(category), "none")
  )

# Keep LOH fragments only
loh_segments <- tracerx_ascat_seg %>%
  filter(LOH == TRUE)

# Create separate dataframes for each category
loh_category_1 <- loh_segments %>% filter(nMajor == 0)
loh_category_2 <- loh_segments %>% filter(nMajor == 1)
loh_category_3 <- loh_segments %>% filter(nMajor == 2)
loh_category_4 <- loh_segments %>% filter(nMajor == 3)
loh_category_5 <- loh_segments %>% filter(nMajor == 4)
loh_category_6 <- loh_segments %>% filter(nMajor == 5)
loh_category_7 <- loh_segments %>% filter(nMajor == 6)
loh_category_8 <- loh_segments %>% filter(nMajor >= 7)

# Summary checks
table(tracerx_ascat_seg$category, useNA = "ifany")
table(loh_segments$category, useNA = "ifany")

# Create processed data directory if needed
dir.create("data/processed", recursive = TRUE, showWarnings = FALSE)

# Save outputs
saveRDS(
  tracerx_ascat_seg,
  "data/processed/tracerx_ascat_seg_categories.rds"
)

saveRDS(
  loh_segments,
  "data/processed/loh_segments_categories.rds"
)
