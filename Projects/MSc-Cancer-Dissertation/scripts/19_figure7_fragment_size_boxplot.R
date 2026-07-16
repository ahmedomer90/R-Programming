# -----------------------------------------------------------------------------
# Original MSc Cancer Dissertation analysis (2021)
# Refactored for reproducibility and GitHub portfolio (2026)
# -----------------------------------------------------------------------------

# 19_figure7_fragment_size_boxplot.R
# Purpose: reproduce Dissertation Figure 7, showing total fragment sizes
# across LOH categories 1–8 and non-LOH fragments ("none")

library(dplyr)
library(ggplot2)
library(readr)

# -----------------------------------------------------------------------------
# Input files
# -----------------------------------------------------------------------------

loh_summary_file <-
  "data/processed/sample_category_fragment_size_summary.rds"

non_loh_file <-
  "data/processed/non_loh_segments.rds"

required_files <- c(
  loh_summary_file,
  non_loh_file
)

missing_files <- required_files[
  !file.exists(required_files)
]

if (length(missing_files) > 0) {
  stop(
    "The following required input file(s) were not found:\n",
    paste(missing_files, collapse = "\n")
  )
}

# -----------------------------------------------------------------------------
# Load validated data
# -----------------------------------------------------------------------------

sample_category_fragment_size <- readRDS(
  loh_summary_file
)

non_loh_segments <- readRDS(
  non_loh_file
)

# -----------------------------------------------------------------------------
# Confirm required columns exist
# -----------------------------------------------------------------------------

required_loh_columns <- c(
  "sample",
  "category",
  "total_fragment_size_bp"
)

if (
  !all(
    required_loh_columns %in%
    names(sample_category_fragment_size)
  )
) {
  stop(
    "The LOH summary is missing one or more required columns: ",
    paste(required_loh_columns, collapse = ", ")
  )
}

required_non_loh_columns <- c(
  "sample",
  "startpos",
  "endpos"
)

if (
  !all(
    required_non_loh_columns %in%
    names(non_loh_segments)
  )
) {
  stop(
    "The non-LOH dataset is missing one or more required columns: ",
    paste(required_non_loh_columns, collapse = ", ")
  )
}

# -----------------------------------------------------------------------------
# Prepare LOH sample-category totals
# -----------------------------------------------------------------------------

loh_plot_data <- sample_category_fragment_size %>%
  dplyr::transmute(
    sample = as.character(sample),
    category = as.character(category),
    total_fragment_size_bp =
      as.numeric(total_fragment_size_bp)
  )

# -----------------------------------------------------------------------------
# Calculate total non-LOH fragment size for each sample
# -----------------------------------------------------------------------------

non_loh_sample_totals <- non_loh_segments %>%
  dplyr::mutate(
    fragment_size_bp = endpos - startpos
  ) %>%
  dplyr::group_by(sample) %>%
  dplyr::summarise(
    total_fragment_size_bp = sum(
      fragment_size_bp,
      na.rm = TRUE
    ),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    sample = as.character(sample),
    category = "none"
  ) %>%
  dplyr::select(
    sample,
    category,
    total_fragment_size_bp
  )

# -----------------------------------------------------------------------------
# Combine LOH categories and non-LOH sample totals
# -----------------------------------------------------------------------------

figure_7_plot_data <- dplyr::bind_rows(
  loh_plot_data,
  non_loh_sample_totals
) %>%
  dplyr::mutate(
    category = factor(
      category,
      levels = c(
        as.character(1:8),
        "none"
      )
    )
  ) %>%
  dplyr::arrange(
    category,
    sample
  )

# Check for invalid totals
if (
  any(
    figure_7_plot_data$total_fragment_size_bp < 0,
    na.rm = TRUE
  )
) {
  stop(
    "One or more negative total fragment sizes were detected."
  )
}

# -----------------------------------------------------------------------------
# Create dissertation-style jittered boxplot
# -----------------------------------------------------------------------------

set.seed(123)

figure_7 <- ggplot(
  figure_7_plot_data,
  aes(
    x = category,
    y = total_fragment_size_bp
  )
) +
  geom_boxplot(
    width = 0.7,
    fill = "white",
    colour = "black",
    linewidth = 0.6,
    outlier.shape = 16,
    outlier.size = 1.5
  ) +
  geom_jitter(
    width = 0.13,
    height = 0,
    shape = 16,
    size = 1.5,
    alpha = 0.25,
    colour = "black"
  ) +
  scale_y_continuous(
    labels = scales::label_scientific(),
    expand = expansion(
      mult = c(0.02, 0.05)
    )
  ) +
  labs(
    title = "Fragment sizes of different categories",
    x = "category",
    y = "Total Fragment Size"
  ) +
  theme_gray(
    base_size = 12
  ) +
  theme(
    plot.title = element_text(
      size = 16,
      face = "plain",
      hjust = 0
    ),
    axis.title.x = element_text(
      size = 12
    ),
    axis.title.y = element_text(
      size = 12
    ),
    axis.text.x = element_text(
      size = 10
    ),
    axis.text.y = element_text(
      size = 10
    ),
    plot.margin = margin(
      t = 10,
      r = 12,
      b = 10,
      l = 10
    )
  )

# -----------------------------------------------------------------------------
# Save figure
# -----------------------------------------------------------------------------

dir.create(
  "figures",
  recursive = TRUE,
  showWarnings = FALSE
)

figure_png <-
  "figures/Figure7_fragment_size_boxplot.png"

figure_pdf <-
  "figures/Figure7_fragment_size_boxplot.pdf"

ggsave(
  filename = figure_png,
  plot = figure_7,
  width = 10,
  height = 7,
  units = "in",
  dpi = 300,
  bg = "white"
)

ggsave(
  filename = figure_pdf,
  plot = figure_7,
  width = 10,
  height = 7,
  units = "in",
  bg = "white"
)

# -----------------------------------------------------------------------------
# Validation
# -----------------------------------------------------------------------------

figure_7_group_summary <- figure_7_plot_data %>%
  dplyr::group_by(category) %>%
  dplyr::summarise(
    observation_count = dplyr::n(),
    zero_total_count = sum(
      total_fragment_size_bp == 0,
      na.rm = TRUE
    ),
    minimum_total_bp = min(
      total_fragment_size_bp,
      na.rm = TRUE
    ),
    median_total_bp = median(
      total_fragment_size_bp,
      na.rm = TRUE
    ),
    maximum_total_bp = max(
      total_fragment_size_bp,
      na.rm = TRUE
    ),
    .groups = "drop"
  )

print(
  figure_7_group_summary,
  n = Inf
)

figure_7_validation <- tibble::tibble(
  LOH_sample_category_rows =
    nrow(loh_plot_data),
  
  non_LOH_sample_rows =
    nrow(non_loh_sample_totals),
  
  total_plot_rows =
    nrow(figure_7_plot_data),
  
  categories_displayed =
    dplyr::n_distinct(
      figure_7_plot_data$category
    ),
  
  negative_total_sizes =
    sum(
      figure_7_plot_data$total_fragment_size_bp < 0,
      na.rm = TRUE
    )
)

print(
  figure_7_validation,
  n = Inf
)

figure_file_validation <- tibble::tibble(
  figure_file = c(
    figure_png,
    figure_pdf
  ),
  file_exists = file.exists(
    c(
      figure_png,
      figure_pdf
    )
  ),
  file_size_bytes = file.info(
    c(
      figure_png,
      figure_pdf
    )
  )$size
)

print(
  figure_file_validation,
  n = Inf
)

if (!all(figure_file_validation$file_exists)) {
  warning(
    "One or more Figure 7 output files were not created."
  )
}

# -----------------------------------------------------------------------------
# Figure legend:
#
# Figure 7. Jittered boxplot showing total fragment sizes across the eight
# categories of LOH and the non-LOH group.
#
# "none" represents the total size of non-LOH fragments within each
# tumour-region sample. Each point represents the total fragment size for
# one sample within the corresponding LOH category or non-LOH group.
# -----------------------------------------------------------------------------
