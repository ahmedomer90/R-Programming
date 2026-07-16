# -----------------------------------------------------------------------------
# Original MSc Cancer Dissertation analysis (2021)
# Refactored for reproducibility and GitHub portfolio (2026)
# -----------------------------------------------------------------------------

# 13_gene_annotation.R
# Purpose: annotate high-prevalence LOH regions with overlapping genes
# using chromosome-level batched Ensembl queries

library(dplyr)
library(readr)
library(biomaRt)
library(IRanges)

# Load validated high-prevalence LOH regions
high_prevalence_regions <- readRDS(
  "data/processed/high_prevalence_LOH_regions.rds"
)

# Confirm required columns exist
required_columns <- c(
  "chr",
  "startpos",
  "endpos"
)

if (!all(required_columns %in% names(high_prevalence_regions))) {
  stop(
    "One or more required columns were not found: ",
    paste(required_columns, collapse = ", ")
  )
}

# Confirm that genomic coordinates are valid
if (
  any(is.na(high_prevalence_regions$startpos)) ||
  any(is.na(high_prevalence_regions$endpos))
) {
  stop("Missing genomic start or end positions were detected.")
}

if (
  any(
    high_prevalence_regions$endpos <
    high_prevalence_regions$startpos
  )
) {
  stop("One or more regions have endpos smaller than startpos.")
}

# Assign a unique identifier to every LOH region
# Convert ASCAT chromosomes 23 and 24 to Ensembl chromosomes X and Y
high_prevalence_regions <- high_prevalence_regions %>%
  mutate(
    region_id = row_number(),
    ensembl_chr = case_when(
      chr == 23 ~ "X",
      chr == 24 ~ "Y",
      TRUE ~ as.character(chr)
    )
  )

# Connect to the Ensembl human gene dataset
ensembl <- useEnsembl(
  biomart = "genes",
  dataset = "hsapiens_gene_ensembl"
)

# Chromosomes represented in the LOH dataset
chromosomes_to_query <- unique(
  high_prevalence_regions$ensembl_chr
)

# Store chromosome-level annotation results
chromosome_annotation_list <- vector(
  mode = "list",
  length = length(chromosomes_to_query)
)

names(chromosome_annotation_list) <- chromosomes_to_query

# Query Ensembl once per chromosome
for (i in seq_along(chromosomes_to_query)) {
  
  current_chr <- chromosomes_to_query[i]
  
  message(
    "Querying chromosome ",
    current_chr,
    " (",
    i,
    " of ",
    length(chromosomes_to_query),
    ")"
  )
  
  chromosome_genes <- getBM(
    attributes = c(
      "hgnc_symbol",
      "ensembl_gene_id",
      "chromosome_name",
      "start_position",
      "end_position",
      "gene_biotype"
    ),
    filters = "chromosome_name",
    values = current_chr,
    mart = ensembl
  ) %>%
    mutate(
      chromosome_name = as.character(chromosome_name),
      start_position = as.numeric(start_position),
      end_position = as.numeric(end_position)
    ) %>%
    filter(
      !is.na(start_position),
      !is.na(end_position),
      start_position <= end_position
    ) %>%
    distinct()
  
  chromosome_regions <- high_prevalence_regions %>%
    filter(ensembl_chr == current_chr)
  
  # Continue safely if no genes or regions are available
  if (
    nrow(chromosome_genes) == 0 ||
    nrow(chromosome_regions) == 0
  ) {
    chromosome_annotation_list[[i]] <- tibble()
    next
  }
  
  # Construct genomic interval objects
  region_ranges <- IRanges(
    start = chromosome_regions$startpos,
    end = chromosome_regions$endpos
  )
  
  gene_ranges <- IRanges(
    start = chromosome_genes$start_position,
    end = chromosome_genes$end_position
  )
  
  # Identify all region-gene overlaps
  overlap_hits <- findOverlaps(
    query = region_ranges,
    subject = gene_ranges,
    type = "any"
  )
  
  if (length(overlap_hits) == 0) {
    chromosome_annotation_list[[i]] <- tibble()
    next
  }
  
  # Match each LOH region with every overlapping gene
  region_rows <- chromosome_regions[
    queryHits(overlap_hits),
    ,
    drop = FALSE
  ]
  
  gene_rows <- chromosome_genes[
    subjectHits(overlap_hits),
    ,
    drop = FALSE
  ]
  
  chromosome_annotation_list[[i]] <- bind_cols(
    region_rows,
    gene_rows
  ) %>%
    mutate(
      chromosome_name = as.character(chromosome_name)
    )
  
  message(
    "  Found ",
    length(overlap_hits),
    " region-gene overlap records"
  )
}

# Standardise chromosome_name type across every list element
chromosome_annotation_list <- lapply(
  chromosome_annotation_list,
  function(x) {
    
    if (nrow(x) == 0) {
      return(x)
    }
    
    x %>%
      mutate(
        chromosome_name = as.character(chromosome_name)
      )
  }
)

# Combine annotations from all chromosomes
gene_annotations <- bind_rows(
  chromosome_annotation_list
) %>%
  filter(
    !is.na(hgnc_symbol),
    hgnc_symbol != ""
  ) %>%
  distinct() %>%
  arrange(
    chr,
    startpos,
    endpos,
    hgnc_symbol
  )

# Identify regions with no annotated HGNC gene
annotated_region_ids <- unique(
  gene_annotations$region_id
)

unannotated_regions <- high_prevalence_regions %>%
  filter(
    !region_id %in% annotated_region_ids
  )

# Validation summaries
cat(
  "Input LOH regions:",
  nrow(high_prevalence_regions),
  "\n"
)

cat(
  "Regions with at least one HGNC gene:",
  n_distinct(gene_annotations$region_id),
  "\n"
)

cat(
  "Regions without an HGNC gene:",
  nrow(unannotated_regions),
  "\n"
)

cat(
  "Annotated region-gene records:",
  nrow(gene_annotations),
  "\n"
)

cat(
  "Unique HGNC gene symbols:",
  n_distinct(gene_annotations$hgnc_symbol),
  "\n"
)

annotation_by_chromosome <- gene_annotations %>%
  count(
    ensembl_chr,
    name = "annotation_count"
  ) %>%
  mutate(
    chromosome_order = match(
      ensembl_chr,
      c(
        as.character(1:22),
        "X",
        "Y"
      )
    )
  ) %>%
  arrange(chromosome_order) %>%
  dplyr::select(-chromosome_order)

print(
  annotation_by_chromosome,
  n = Inf
)

# Create output directories if necessary
dir.create(
  "data/processed",
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  "results",
  recursive = TRUE,
  showWarnings = FALSE
)

# Save outputs
saveRDS(
  gene_annotations,
  "data/processed/high_prevalence_LOH_gene_annotations.rds"
)

saveRDS(
  unannotated_regions,
  "data/processed/unannotated_high_prevalence_LOH_regions.rds"
)

write_csv(
  gene_annotations,
  "results/high_prevalence_LOH_gene_annotations.csv"
)

write_csv(
  unannotated_regions,
  "results/unannotated_high_prevalence_LOH_regions.csv"
)

# -----------------------------------------------------------------------------
# Validation note:
# Ensembl genes are downloaded once per chromosome rather than once per LOH
# region. Region-gene overlaps are calculated locally using IRanges.
# ASCAT chromosomes 23 and 24 are converted to Ensembl chromosomes X and Y.
# The chromosome_name column is converted to character after every Ensembl
# query so that chromosome-level results can be combined safely.
# Regions without an HGNC gene annotation are retained in separate outputs.
# -----------------------------------------------------------------------------