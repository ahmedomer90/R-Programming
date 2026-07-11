# 12_gene_annotation.R
# MSc Cancer Dissertation reconstruction
# Purpose: annotate high-prevalence LOH regions with gene information

# Load packages
library(dplyr)
library(readr)
library(biomaRt)

# Load high-prevalence LOH regions
high_prevalence_regions <- readRDS(
  "data/processed/high_prevalence_LOH_regions.rds"
)

# Connect to Ensembl
ensembl <- useEnsembl(
  biomart = "genes",
  dataset = "hsapiens_gene_ensembl"
)

# Function to retrieve genes overlapping a genomic region
get_region_genes <- function(chr, startpos, endpos) {
  getBM(
    attributes = c(
      "hgnc_symbol",
      "ensembl_gene_id",
      "chromosome_name",
      "start_position",
      "end_position",
      "gene_biotype"
    ),
    filters = c("chromosome_name", "start", "end"),
    values = list(chr, startpos, endpos),
    mart = ensembl
  )
}

# Annotate each high-prevalence LOH region
gene_annotations <- high_prevalence_regions %>%
  rowwise() %>%
  mutate(
    genes = list(get_region_genes(chr, startpos, endpos))
  ) %>%
  tidyr::unnest(genes) %>%
  ungroup()

# Remove rows without gene symbols
gene_annotations <- gene_annotations %>%
  filter(!is.na(hgnc_symbol), hgnc_symbol != "")

# Save outputs
saveRDS(
  gene_annotations,
  "data/processed/high_prevalence_LOH_gene_annotations.rds"
)

write_csv(
  gene_annotations,
  "results/high_prevalence_LOH_gene_annotations.csv"
)
