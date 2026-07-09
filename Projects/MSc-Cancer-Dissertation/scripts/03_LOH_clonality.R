# 03_LOH_clonality.R
# MSc Cancer Dissertation reconstruction
# Purpose: classify LOH fragments according to clonality

# Load packages
library(dplyr)

# Load LOH segments
loh_segments <- readRDS("data/processed/loh_segments.rds")

# Classify clonality
# Clonal = all tumour cells
# Subclonal = subset of tumour cells

loh_segments <- loh_segments %>%
  mutate(
    clonality = ifelse(frac1_A == 1, "Clonal", "Subclonal")
  )

# Summary
table(loh_segments$clonality)

# Save output
saveRDS(
  loh_segments,
  "data/processed/loh_segments_clonality.rds"
)
