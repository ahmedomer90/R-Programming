# -----------------------------------------------------------------------------
# Original MSc Cancer Dissertation analysis (2021)
# Refactored for reproducibility and GitHub portfolio (2026)
# -----------------------------------------------------------------------------

# 01_data_import.R
# MSc Cancer Dissertation reconstruction
# Purpose: import TRACERx ASCAT segment data and clinical data

library(tidyverse)
library(readxl)

ascat_file <- "data/raw/tracerx.ascat.seg.129samples.r002.2.20160818.RData"
clinical_file <- "data/raw/trx_clin_stab2.xls"

# Load ASCAT RData file
load(ascat_file)

# Display objects loaded from the .RData file
print(ls())

# Import clinical data
clinical_data <- read_excel(clinical_file)

# Inspect data
glimpse(tracerx.ascat.seg)
glimpse(clinical_data)

# Basic checks
nrow(tracerx.ascat.seg)
ncol(tracerx.ascat.seg)

nrow(clinical_data)
ncol(clinical_data)

names(tracerx.ascat.seg)
names(clinical_data)

# Rename object for cleaner script naming
tracerx_ascat_seg <- tracerx.ascat.seg

# Create processed data directory if it does not exist
dir.create("data/processed", recursive = TRUE, showWarnings = FALSE)

# Save imported objects
saveRDS(tracerx_ascat_seg, "data/processed/tracerx_ascat_seg_imported.rds")
saveRDS(clinical_data, "data/processed/clinical_data_imported.rds")
