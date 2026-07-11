# -----------------------------------------------------------------------------
# Original MSc Cancer Dissertation analysis (2021)
# Refactored for reproducibility and GitHub portfolio (2026)
# -----------------------------------------------------------------------------

# 15_visualisation_LOH_patterns.R
# Purpose: generate plots showing LOH category, clonality, chromosome,
# sample and fragment-size patterns

library(dplyr)
library(readr)
library(ggplot2)

# -----------------------------------------------------------------------------
# Load validated processed data
# -----------------------------------------------------------------------------

loh_segments <- readRDS(
  "data/processed/loh_segments_fragment_sizes.rds"
)

category_clonality_summary <- readRDS(
  "data/processed/category_clonality_summary.rds"
)

loh_chromosome_clonality_summary <- readRDS(
  "data/processed/loh_chromosome_clonality_summary.rds"
)

sample_clonality_summary <- readRDS(
  "data/processed/sample_clonality_summary.rds"
)

# -----------------------------------------------------------------------------
# Confirm required columns exist
# -----------------------------------------------------------------------------

required_loh_columns <- c(
  "category",
  "clonality",
  "fragment_size_mb"
)

if (!all(required_loh_columns %in% names(loh_segments))) {
  stop(
    "The LOH dataset is missing one or more required columns: ",
    paste(required_loh_columns, collapse = ", ")
  )
}

required_category_columns <- c(
  "category",
  "clonality",
  "count",
  "proportion"
)

if (
  !all(
    required_category_columns %in%
    names(category_clonality_summary)
  )
) {
  stop(
    "The category-clonality summary is missing one or more ",
    "required columns."
  )
}

required_chromosome_columns <- c(
  "chr",
  "clonality",
  "count",
  "proportion"
)

if (
  !all(
    required_chromosome_columns %in%
    names(loh_chromosome_clonality_summary)
  )
) {
  stop(
    "The chromosome-clonality summary is missing one or more ",
    "required columns."
  )
}

required_sample_columns <- c(
  "sample",
  "clonality",
  "count",
  "proportion"
)

if (
  !all(
    required_sample_columns %in%
    names(sample_clonality_summary)
  )
) {
  stop(
    "The sample-clonality summary is missing one or more ",
    "required columns."
  )
}

# -----------------------------------------------------------------------------
# Prepare variables for plotting
# -----------------------------------------------------------------------------

loh_segments <- loh_segments %>%
  dplyr::mutate(
    category = factor(
      category,
      levels = as.character(1:8)
    )
  )

category_clonality_summary <-
  category_clonality_summary %>%
  dplyr::mutate(
    category = factor(
      category,
      levels = as.character(1:8)
    ),
    clonality = factor(
      clonality,
      levels = c(
        "clonal",
        "subclonal",
        "none"
      )
    )
  )

loh_chromosome_clonality_summary <-
  loh_chromosome_clonality_summary %>%
  dplyr::mutate(
    chr = factor(
      chr,
      levels = 1:24,
      labels = c(
        as.character(1:22),
        "X",
        "Y"
      )
    ),
    clonality = factor(
      clonality,
      levels = c(
        "clonal",
        "subclonal",
        "none"
      )
    )
  )

sample_clonality_summary <-
  sample_clonality_summary %>%
  dplyr::mutate(
    clonality = factor(
      clonality,
      levels = c(
        "clonal",
        "subclonal",
        "none"
      )
    )
  )

# Create figures directory if needed
dir.create(
  "figures",
  recursive = TRUE,
  showWarnings = FALSE
)

# -----------------------------------------------------------------------------
# Plot 1: LOH category distribution
# -----------------------------------------------------------------------------

plot_category_distribution <- loh_segments %>%
  ggplot(
    aes(x = category)
  ) +
  geom_bar() +
  labs(
    title = "Distribution of LOH categories",
    x = "LOH category",
    y = "Number of LOH fragments"
  ) +
  theme_bw()

ggsave(
  filename = "figures/LOH_category_distribution.png",
  plot = plot_category_distribution,
  width = 8,
  height = 5,
  dpi = 300
)

# -----------------------------------------------------------------------------
# Plot 2: LOH category by clonality proportion
# -----------------------------------------------------------------------------

plot_category_clonality <-
  category_clonality_summary %>%
  ggplot(
    aes(
      x = category,
      y = proportion,
      fill = clonality
    )
  ) +
  geom_col() +
  scale_y_continuous(
    labels = scales::percent
  ) +
  labs(
    title = "Clonality proportions by LOH category",
    x = "LOH category",
    y = "Proportion",
    fill = "Clonality"
  ) +
  theme_bw()

ggsave(
  filename =
    "figures/LOH_category_clonality_proportion.png",
  plot = plot_category_clonality,
  width = 8,
  height = 5,
  dpi = 300
)

# -----------------------------------------------------------------------------
# Plot 3: Chromosome by clonality proportion
# -----------------------------------------------------------------------------

plot_chromosome_clonality <-
  loh_chromosome_clonality_summary %>%
  ggplot(
    aes(
      x = chr,
      y = proportion,
      fill = clonality
    )
  ) +
  geom_col() +
  scale_y_continuous(
    labels = scales::percent
  ) +
  labs(
    title = "Clonality proportions by chromosome",
    x = "Chromosome",
    y = "Proportion",
    fill = "Clonality"
  ) +
  theme_bw()

ggsave(
  filename =
    "figures/chromosome_clonality_proportion.png",
  plot = plot_chromosome_clonality,
  width = 10,
  height = 5,
  dpi = 300
)

# -----------------------------------------------------------------------------
# Plot 4: Fragment size by LOH category
# -----------------------------------------------------------------------------

plot_fragment_size <- loh_segments %>%
  ggplot(
    aes(
      x = category,
      y = fragment_size_mb
    )
  ) +
  geom_boxplot(
    outlier.alpha = 0.2
  ) +
  labs(
    title = "Fragment-size distribution by LOH category",
    x = "LOH category",
    y = "Fragment size (Mb)"
  ) +
  theme_bw()

ggsave(
  filename =
    "figures/fragment_size_by_LOH_category.png",
  plot = plot_fragment_size,
  width = 8,
  height = 5,
  dpi = 300
)

# -----------------------------------------------------------------------------
# Plot 5: Sample by clonality proportion
# -----------------------------------------------------------------------------

plot_sample_clonality <-
  sample_clonality_summary %>%
  ggplot(
    aes(
      x = sample,
      y = proportion,
      fill = clonality
    )
  ) +
  geom_col() +
  scale_y_continuous(
    labels = scales::percent
  ) +
  labs(
    title = "Clonality proportions across tumour-region samples",
    x = "Tumour-region sample",
    y = "Proportion",
    fill = "Clonality"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )

ggsave(
  filename =
    "figures/sample_clonality_proportion.png",
  plot = plot_sample_clonality,
  width = 12,
  height = 6,
  dpi = 300
)

# -----------------------------------------------------------------------------
# Validation output
# -----------------------------------------------------------------------------

expected_figure_files <- c(
  "figures/LOH_category_distribution.png",
  "figures/LOH_category_clonality_proportion.png",
  "figures/chromosome_clonality_proportion.png",
  "figures/fragment_size_by_LOH_category.png",
  "figures/sample_clonality_proportion.png"
)

figure_validation <- tibble::tibble(
  figure_file = expected_figure_files,
  file_exists = file.exists(expected_figure_files),
  file_size_bytes = file.info(expected_figure_files)$size
)

print(
  figure_validation,
  n = Inf
)

if (!all(figure_validation$file_exists)) {
  warning(
    "One or more expected figure files were not created."
  )
}

# -----------------------------------------------------------------------------
# Validation note:
# Figures are generated from the validated outputs of Scripts 5, 7, 8 and 9.
# LOH categories are ordered from 1 to 8, chromosomes 23 and 24 are displayed
# as X and Y, and rare clonality values classified as "none" are retained.
# -----------------------------------------------------------------------------
