# 02_classify_LOH.R
# MSc Cancer Dissertation reconstruction
# Purpose: classify ASCAT copy-number fragments as LOH or non-LOH

# Load packages
library(dplyr)
library(readr)

# Import processed ASCAT segment data
tracerx_ascat_seg <- readRDS("data/processed/tracerx_ascat_seg_imported.rds")

# Classify fragments as LOH or non-LOH
# LOH is defined as nMinor == 0
tracerx_ascat_seg <- tracerx_ascat_seg %>%
  mutate(
    LOH = nMinor == 0
  )

# Check LOH classification
table(tracerx_ascat_seg$LOH)

# Create separate dataframes
loh_segments <- tracerx_ascat_seg %>%
  filter(LOH == TRUE)

non_loh_segments <- tracerx_ascat_seg %>%
  filter(LOH == FALSE)

# Basic checks
nrow(tracerx_ascat_seg)
nrow(loh_segments)
nrow(non_loh_segments)

# Save outputs
saveRDS(tracerx_ascat_seg, "data/processed/tracerx_ascat_seg_LOH_classified.rds")
saveRDS(loh_segments, "data/processed/loh_segments.rds")
saveRDS(non_loh_segments, "data/processed/non_loh_segments.rds")
