# 01_data_import.R
# MSc Cancer Dissertation reconstruction
# Purpose: import TRACERx ASCAT segment data and clinical data

# Load packages
library(tidyverse)
library(readxl)

# Define file paths
ascat_file <- "data/raw/tracerx_ascat_seg.csv"
clinical_file <- "data/raw/tracerx_clinical_data.xlsx"

# Import data
tracerx_ascat_seg <- read_csv(ascat_file)
clinical_data <- read_excel(clinical_file)

# Inspect data
glimpse(tracerx_ascat_seg)
glimpse(clinical_data)

# Basic checks
nrow(tracerx_ascat_seg)
ncol(tracerx_ascat_seg)

nrow(clinical_data)
ncol(clinical_data)

# Check expected columns
names(tracerx_ascat_seg)
names(clinical_data)

# Save imported objects
saveRDS(tracerx_ascat_seg, "data/processed/tracerx_ascat_seg_imported.rds")
saveRDS(clinical_data, "data/processed/clinical_data_imported.rds")
