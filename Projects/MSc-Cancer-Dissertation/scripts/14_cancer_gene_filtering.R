# 13_cancer_gene_filtering.R
# MSc Cancer Dissertation reconstruction
# Purpose: filter LOH-affected genes using a cancer driver gene list

# Load packages
library(dplyr)
library(readr)

# Load annotated genes from high-prevalence LOH regions
gene_annotations <- readRDS(
  "data/processed/high_prevalence_LOH_gene_annotations.rds"
)

# Load cancer driver gene list
# Expected columns may include: hgnc_symbol, role_in_cancer, tumour_types
cancer_gene_list <- read_csv("data/raw/cancer_driver_gene_list.csv")

# Filter annotated genes to retain cancer driver genes
loh_cancer_genes <- gene_annotations %>%
  inner_join(
    cancer_gene_list,
    by = "hgnc_symbol"
  )

# Remove duplicate gene entries if present
loh_cancer_genes_unique <- loh_cancer_genes %>%
  distinct(hgnc_symbol, .keep_all = TRUE)

# Save outputs
saveRDS(
  loh_cancer_genes,
  "data/processed/loh_cancer_genes.rds"
)

write_csv(
  loh_cancer_genes,
  "results/loh_cancer_genes.csv"
)

write_csv(
  loh_cancer_genes_unique,
  "results/loh_cancer_genes_unique.csv"
)
