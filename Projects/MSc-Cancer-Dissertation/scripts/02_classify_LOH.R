# -----------------------------------------------------------------------------
# Original MSc Cancer Dissertation analysis (2021)
# Refactored for reproducibility and GitHub portfolio (2026)
# -----------------------------------------------------------------------------

# 02_classify_LOH.R
# Purpose: classify ASCAT copy-number fragments as LOH or non-LOH

library(dplyr)

# Load imported ASCAT segment data
tracerx_ascat_seg <- readRDS("data/processed/tracerx_ascat_seg_imported.rds")

# Classify fragments as LOH
# LOH is defined as nMinor == 0
tracerx_ascat_seg <- tracerx_ascat_seg %>%
  mutate(
    LOH = nMinor == 0
  )

# Basic check
table(tracerx_ascat_seg$LOH)

# Create separate LOH and non-LOH dataframes
loh_segments <- tracerx_ascat_seg %>%
  filter(LOH == TRUE)

non_loh_segments <- tracerx_ascat_seg %>%
  filter(LOH == FALSE)

# Create processed data directory if it does not exist
dir.create("data/processed", recursive = TRUE, showWarnings = FALSE)

# Save outputs
saveRDS(tracerx_ascat_seg, "data/processed/tracerx_ascat_seg_LOH_classified.rds")
saveRDS(loh_segments, "data/processed/loh_segments.rds")
saveRDS(non_loh_segments, "data/processed/non_loh_segments.rds")
