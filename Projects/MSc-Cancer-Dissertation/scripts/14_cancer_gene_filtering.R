# -----------------------------------------------------------------------------
# Original MSc Cancer Dissertation analysis (2021)
# Refactored for reproducibility and GitHub portfolio (2026)
# -----------------------------------------------------------------------------

# 14_cancer_gene_filtering.R
# Purpose: identify LOH-affected cancer driver genes and add
# Cancer Gene Census annotations

library(dplyr)
library(readr)
library(stringr)

# -----------------------------------------------------------------------------
# Input files
# -----------------------------------------------------------------------------

gene_annotation_file <-
  "data/processed/high_prevalence_LOH_gene_annotations.rds"

driver_gene_file <-
  "data/raw/driverGenes_120516_unique.txt"

cancer_gene_census_file <-
  "data/raw/cancer_gene_census.csv"

# Confirm that all required files exist
required_files <- c(
  gene_annotation_file,
  driver_gene_file,
  cancer_gene_census_file
)

missing_files <- required_files[
  !file.exists(required_files)
]

if (length(missing_files) > 0) {
  stop(
    "The following required file(s) were not found:\n",
    paste(missing_files, collapse = "\n")
  )
}

# -----------------------------------------------------------------------------
# Load gene annotations created by Script 13
# -----------------------------------------------------------------------------

gene_annotations <- readRDS(
  gene_annotation_file
)

if (!"hgnc_symbol" %in% names(gene_annotations)) {
  stop(
    "The Script 13 output does not contain ",
    "a column named 'hgnc_symbol'."
  )
}

gene_annotations <- gene_annotations %>%
  dplyr::mutate(
    hgnc_symbol = stringr::str_trim(
      as.character(hgnc_symbol)
    )
  ) %>%
  dplyr::filter(
    !is.na(hgnc_symbol),
    hgnc_symbol != ""
  )

# -----------------------------------------------------------------------------
# Load the curated cancer driver-gene list
# -----------------------------------------------------------------------------

# The file contains comment lines beginning with "#", followed by
# one gene symbol per line.
driver_gene_list <- readr::read_lines(
  driver_gene_file
) %>%
  stringr::str_trim()

driver_gene_list <- driver_gene_list[
  driver_gene_list != "" &
    !stringr::str_starts(driver_gene_list, "#")
]

driver_gene_reference <- tibble::tibble(
  hgnc_symbol = driver_gene_list
) %>%
  dplyr::mutate(
    hgnc_symbol = stringr::str_trim(hgnc_symbol)
  ) %>%
  dplyr::filter(
    !is.na(hgnc_symbol),
    hgnc_symbol != ""
  ) %>%
  dplyr::distinct(hgnc_symbol)

# -----------------------------------------------------------------------------
# Load and standardise the Cancer Gene Census table
# -----------------------------------------------------------------------------

cancer_gene_census <- readr::read_csv(
  cancer_gene_census_file,
  show_col_types = FALSE
)

required_census_columns <- c(
  "Gene Symbol",
  "Name",
  "Genome Location",
  "Tier",
  "Hallmark",
  "Tumour Types(Somatic)",
  "Tumour Types(Germline)",
  "Cancer Syndrome",
  "Tissue Type",
  "Molecular Genetics",
  "Role in Cancer",
  "Mutation Types"
)

missing_census_columns <- setdiff(
  required_census_columns,
  names(cancer_gene_census)
)

if (length(missing_census_columns) > 0) {
  stop(
    "The Cancer Gene Census file is missing the following column(s):\n",
    paste(missing_census_columns, collapse = "\n")
  )
}

cancer_gene_census_clean <- cancer_gene_census %>%
  dplyr::transmute(
    hgnc_symbol = stringr::str_trim(
      as.character(`Gene Symbol`)
    ),
    gene_name = Name,
    census_genome_location = `Genome Location`,
    census_tier = Tier,
    census_hallmark = Hallmark,
    tumour_types_somatic = `Tumour Types(Somatic)`,
    tumour_types_germline = `Tumour Types(Germline)`,
    cancer_syndrome = `Cancer Syndrome`,
    tissue_type = `Tissue Type`,
    molecular_genetics = `Molecular Genetics`,
    role_in_cancer = `Role in Cancer`,
    mutation_types = `Mutation Types`
  ) %>%
  dplyr::filter(
    !is.na(hgnc_symbol),
    hgnc_symbol != ""
  ) %>%
  dplyr::distinct()

# -----------------------------------------------------------------------------
# Identify LOH-affected genes present in the curated driver-gene list
# -----------------------------------------------------------------------------

loh_driver_genes <- gene_annotations %>%
  dplyr::inner_join(
    driver_gene_reference,
    by = "hgnc_symbol"
  ) %>%
  dplyr::distinct()

# Add Cancer Gene Census annotations where available
loh_cancer_genes <- loh_driver_genes %>%
  dplyr::left_join(
    cancer_gene_census_clean,
    by = "hgnc_symbol"
  ) %>%
  dplyr::mutate(
    present_in_cancer_gene_census =
      !is.na(census_tier)
  ) %>%
  dplyr::distinct()

# -----------------------------------------------------------------------------
# Create a one-row-per-gene summary
# -----------------------------------------------------------------------------

loh_cancer_gene_summary <- loh_cancer_genes %>%
  dplyr::group_by(hgnc_symbol) %>%
  dplyr::summarise(
    affected_region_count =
      dplyr::n_distinct(region_id),
    
    chromosome_count =
      dplyr::n_distinct(ensembl_chr),
    
    maximum_sample_count =
      max(sample_count, na.rm = TRUE),
    
    total_annotation_records =
      dplyr::n(),
    
    present_in_cancer_gene_census =
      any(present_in_cancer_gene_census),
    
    census_tier =
      paste(
        sort(unique(stats::na.omit(census_tier))),
        collapse = "; "
      ),
    
    role_in_cancer =
      paste(
        sort(unique(stats::na.omit(role_in_cancer))),
        collapse = "; "
      ),
    
    tumour_types_somatic =
      paste(
        sort(unique(stats::na.omit(tumour_types_somatic))),
        collapse = "; "
      ),
    
    tumour_types_germline =
      paste(
        sort(unique(stats::na.omit(tumour_types_germline))),
        collapse = "; "
      ),
    
    mutation_types =
      paste(
        sort(unique(stats::na.omit(mutation_types))),
        collapse = "; "
      ),
    
    .groups = "drop"
  ) %>%
  dplyr::arrange(
    dplyr::desc(affected_region_count),
    hgnc_symbol
  )

# Create a one-row-per-gene reference table with selected census details
loh_cancer_genes_unique <- loh_cancer_genes %>%
  dplyr::arrange(
    hgnc_symbol,
    dplyr::desc(present_in_cancer_gene_census)
  ) %>%
  dplyr::distinct(
    hgnc_symbol,
    .keep_all = TRUE
  )

# -----------------------------------------------------------------------------
# Validation output
# -----------------------------------------------------------------------------

cat(
  "Annotated region-gene records:",
  nrow(gene_annotations),
  "\n"
)

cat(
  "Unique annotated HGNC genes:",
  dplyr::n_distinct(gene_annotations$hgnc_symbol),
  "\n"
)

cat(
  "Genes in curated driver-gene list:",
  nrow(driver_gene_reference),
  "\n"
)

cat(
  "Genes in Cancer Gene Census:",
  dplyr::n_distinct(
    cancer_gene_census_clean$hgnc_symbol
  ),
  "\n"
)

cat(
  "Matched LOH driver-gene records:",
  nrow(loh_cancer_genes),
  "\n"
)

cat(
  "Unique matched LOH driver genes:",
  dplyr::n_distinct(
    loh_cancer_genes$hgnc_symbol
  ),
  "\n"
)

cat(
  "Matched genes present in Cancer Gene Census:",
  sum(
    loh_cancer_gene_summary$
      present_in_cancer_gene_census
  ),
  "\n"
)

print(
  head(
    loh_cancer_gene_summary,
    20
  )
)

# -----------------------------------------------------------------------------
# Create output directories
# -----------------------------------------------------------------------------

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

# -----------------------------------------------------------------------------
# Save outputs
# -----------------------------------------------------------------------------

saveRDS(
  loh_cancer_genes,
  "data/processed/loh_cancer_genes.rds"
)

saveRDS(
  loh_cancer_gene_summary,
  "data/processed/loh_cancer_gene_summary.rds"
)

write_csv(
  loh_cancer_genes,
  "results/loh_cancer_genes.csv"
)

write_csv(
  loh_cancer_genes_unique,
  "results/loh_cancer_genes_unique.csv"
)

write_csv(
  loh_cancer_gene_summary,
  "results/loh_cancer_gene_summary.csv"
)

# -----------------------------------------------------------------------------
# Validation note:
# LOH-affected genes are first filtered against the curated
# driverGenes_120516_unique reference list. Cancer Gene Census information is
# then added where available. The complete region-gene table is retained,
# together with separate one-row-per-gene and gene-level summary outputs.
# -----------------------------------------------------------------------------