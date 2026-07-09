# 04_LOH_categories.R
# MSc Cancer Dissertation reconstruction
# Purpose: classify LOH fragments into copy-number categories based on nMajor

# Load packages
library(dplyr)

# Load LOH segments with clonality classification
loh_segments <- readRDS("data/processed/loh_segments_clonality.rds")

# Classify LOH categories based on nMajor
# Category 1: nMajor = 0
# Category 2: nMajor = 1
# Category 3: nMajor = 2
# Category 4: nMajor = 3
# Category 5: nMajor = 4
# Category 6: nMajor = 5
# Category 7: nMajor = 6
# Category 8: nMajor >= 7

loh_segments <- loh_segments %>%
  mutate(
    LOH_category = case_when(
      nMajor == 0 ~ "Category 1",
      nMajor == 1 ~ "Category 2",
      nMajor == 2 ~ "Category 3",
      nMajor == 3 ~ "Category 4",
      nMajor == 4 ~ "Category 5",
      nMajor == 5 ~ "Category 6",
      nMajor == 6 ~ "Category 7",
      nMajor >= 7 ~ "Category 8",
      TRUE ~ NA_character_
    )
  )

# Check category distribution
table(loh_segments$LOH_category, useNA = "ifany")

# Save output
saveRDS(
  loh_segments,
  "data/processed/loh_segments_categories.rds"
)
