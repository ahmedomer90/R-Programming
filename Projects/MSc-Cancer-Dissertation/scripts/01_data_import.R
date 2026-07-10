# 01_data_import.R
# MSc Cancer Dissertation reconstruction
# Purpose: import TRACERx ASCAT segment data and clinical data

library(tidyverse)
library(readxl)

ascat_file <- "data/raw/tracerx.ascat.seg.129samples.r002.2.20160818.RData"
clinical_file <- "data/raw/trx_clin_stab2.xls"

# Load ASCAT RData file
load(ascat_file)

# Check which objects were loaded
ls()

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

dir.create("data/processed", recursive = TRUE, showWarnings = FALSE)

# Save imported objects
saveRDS(tracerx_ascat_seg, "data/processed/tracerx_ascat_seg_imported.rds")
saveRDS(clinical_data, "data/processed/clinical_data_imported.rds")
