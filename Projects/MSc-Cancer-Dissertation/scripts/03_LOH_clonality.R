# 03_LOH_clonality.R
# MSc Cancer Dissertation reconstruction
# Purpose: classify LOH fragments according to clonality using overlapping fragments

# Load packages
library(dplyr)

# Load ASCAT segment data with LOH classification
tracerx_ascat_seg <- readRDS("data/processed/tracerx_ascat_seg_LOH_classified.rds")

# Function to find fragments from the same patient and chromosome
# that overlap the current fragment
get_overlapping_fragments <- function(row_number, ascat_data) {
  
  input_chr <- ascat_data$chr[row_number]
  input_start <- ascat_data$startpos[row_number]
  input_end <- ascat_data$endpos[row_number]
  input_patient <- substr(ascat_data$sample[row_number], 3, 8)
  
  overlapping_fragments <- ascat_data %>%
    mutate(patient = substr(sample, 3, 8)) %>%
    filter(chr == input_chr) %>%
    filter(patient == input_patient) %>%
    filter(startpos < input_end) %>%
    filter(endpos > input_start)
  
  return(overlapping_fragments)
}

# Function to classify clonality using the original MSc logic
classify_LOH_clonality <- function(overlapping_fragments) {
  
  minor_allele_values <- overlapping_fragments$nMinor
  LOH_status <- minor_allele_values == 0
  
  clonality <- ifelse(
    all(LOH_status),
    "clonal",
    ifelse(any(LOH_status), "subclonal", "none")
  )
  
  return(clonality)
}

# Apply clonality classification to every fragment
clonality_vector <- character(nrow(tracerx_ascat_seg))

for (i in seq_len(nrow(tracerx_ascat_seg))) {
  clonality_vector[i] <- classify_LOH_clonality(
    get_overlapping_fragments(i, tracerx_ascat_seg)
  )
}

# Add clonality column
tracerx_ascat_seg$clonality <- clonality_vector

# Keep LOH fragments only for downstream LOH analyses
loh_segments <- tracerx_ascat_seg %>%
  filter(LOH == TRUE)

# Summary
table(loh_segments$clonality)

# Create processed data directory if it does not exist
dir.create("data/processed", recursive = TRUE, showWarnings = FALSE)

# Save outputs
saveRDS(tracerx_ascat_seg, "data/processed/tracerx_ascat_seg_clonality.rds")
saveRDS(loh_segments, "data/processed/loh_segments_clonality.rds")
