# -----------------------------------------------------------------------------
# Original MSc Cancer Dissertation analysis (2021)
# Refactored for reproducibility and GitHub portfolio (2026)
# -----------------------------------------------------------------------------

# 30_figure19_sample_clonality_counts_area.R
# Purpose: reproduce Dissertation Figure 19 as a stacked area chart showing
# the number of clonal and subclonal LOH fragments within each tumour-region
# sample

library(dplyr)
library(tidyr)
library(ggplot2)

# -----------------------------------------------------------------------------
# Input file
# -----------------------------------------------------------------------------

sample_summary_file <-
  "data/processed/sample_clonality_summary.rds"

if (!file.exists(sample_summary_file)) {
  stop(
    "Required input file was not found: ",
    sample_summary_file
  )
}

# -----------------------------------------------------------------------------
# Load validated sample-clonality summary
# -----------------------------------------------------------------------------

sample_clonality_summary <- readRDS(
  sample_summary_file
)

required_columns <- c(
  "sample",
  "clonality",
  "count"
)

missing_columns <- setdiff(
  required_columns,
  names(sample_clonality_summary)
)

if (length(missing_columns) > 0) {
  stop(
    "The sample-clonality summary is missing the following column(s): ",
    paste(
      missing_columns,
      collapse = ", "
    )
  )
}

# -----------------------------------------------------------------------------
# Validate source data
# -----------------------------------------------------------------------------

if (any(is.na(sample_clonality_summary$sample))) {
  stop(
    "Missing sample identifiers were detected."
  )
}

if (any(is.na(sample_clonality_summary$count))) {
  stop(
    "Missing fragment counts were detected."
  )
}

if (any(sample_clonality_summary$count < 0)) {
  stop(
    "Negative fragment counts were detected."
  )
}

unexpected_clonality <- setdiff(
  unique(
    as.character(
      sample_clonality_summary$clonality
    )
  ),
  c(
    "clonal",
    "subclonal",
    "none"
  )
)

if (length(unexpected_clonality) > 0) {
  stop(
    "Unexpected clonality value(s) detected: ",
    paste(
      unexpected_clonality,
      collapse = ", "
    )
  )
}

# -----------------------------------------------------------------------------
# Preserve validated sample order
# -----------------------------------------------------------------------------

sample_order <- sample_clonality_summary %>%
  dplyr::distinct(sample) %>%
  dplyr::pull(sample) %>%
  as.character()

# -----------------------------------------------------------------------------
# Prepare plotting data
# -----------------------------------------------------------------------------

# Figure 19 compares the biological LOH clonality classes only.
# Rare records classified as "none" are excluded because they arose from
# strict interval-overlap edge cases.

figure_19_plot_data <- sample_clonality_summary %>%
  dplyr::filter(
    clonality %in% c(
      "clonal",
      "subclonal"
    )
  ) %>%
  dplyr::select(
    sample,
    clonality,
    count
  ) %>%
  tidyr::complete(
    sample = sample_order,
    clonality = c(
      "subclonal",
      "clonal"
    ),
    fill = list(
      count = 0
    )
  ) %>%
  dplyr::mutate(
    sample = factor(
      sample,
      levels = sample_order
    ),
    
    sample_index = match(
      as.character(sample),
      sample_order
    ),
    
    # Combined with reverse stacking below, this produces:
    # bottom: subclonal
    # top: clonal
    clonality = factor(
      clonality,
      levels = c(
        "subclonal",
        "clonal"
      )
    )
  ) %>%
  dplyr::arrange(
    sample_index,
    clonality
  )

# -----------------------------------------------------------------------------
# Validate prepared plotting data
# -----------------------------------------------------------------------------

if (any(is.na(figure_19_plot_data$count))) {
  stop(
    "Missing counts were produced in the plotting data."
  )
}

if (any(figure_19_plot_data$count < 0)) {
  stop(
    "Negative counts were detected in the plotting data."
  )
}

# -----------------------------------------------------------------------------
# Dissertation-style colours
# -----------------------------------------------------------------------------

clonality_colours <- c(
  subclonal = "#56B937",
  clonal = "#F06F67"
)

# -----------------------------------------------------------------------------
# Create stacked area chart
# -----------------------------------------------------------------------------

figure_19 <- ggplot(
  figure_19_plot_data,
  aes(
    x = sample_index,
    y = count,
    fill = clonality
  )
) +
  geom_area(
    position = position_stack(
      reverse = TRUE
    ),
    colour = NA
  ) +
  scale_fill_manual(
    values = clonality_colours,
    breaks = c(
      "clonal",
      "subclonal"
    ),
    labels = c(
      "clonal",
      "subclonal"
    ),
    name = "clonality",
    drop = FALSE
  ) +
  scale_x_continuous(
    breaks = NULL,
    expand = expansion(
      mult = c(0, 0)
    )
  ) +
  scale_y_continuous(
    breaks = c(
      0,
      50,
      100,
      150
    ),
    limits = c(
      0,
      150
    ),
    expand = expansion(
      mult = c(0, 0)
    )
  ) +
  labs(
    title = paste(
      "number of LOH fragments and proportion of clonality",
      "in each sample"
    ),
    x = "sample",
    y = "number of fragments"
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
    
    axis.text.x = element_blank(),
    
    axis.ticks.x = element_blank(),
    
    axis.text.y = element_text(
      size = 10
    ),
    
    legend.title = element_text(
      size = 11
    ),
    
    legend.text = element_text(
      size = 10
    ),
    
    legend.position = "right",
    
    plot.margin = margin(
      t = 10,
      r = 15,
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
  "figures/Figure19_sample_clonality_counts_area.png"

figure_pdf <-
  "figures/Figure19_sample_clonality_counts_area.pdf"

ggsave(
  filename = figure_png,
  plot = figure_19,
  width = 11,
  height = 7,
  units = "in",
  dpi = 300,
  bg = "white"
)

ggsave(
  filename = figure_pdf,
  plot = figure_19,
  width = 11,
  height = 7,
  units = "in",
  bg = "white"
)

# -----------------------------------------------------------------------------
# Validation
# -----------------------------------------------------------------------------

sample_count_summary <- figure_19_plot_data %>%
  dplyr::group_by(sample) %>%
  dplyr::summarise(
    subclonal_LOH_count = sum(
      count[
        clonality == "subclonal"
      ]
    ),
    
    clonal_LOH_count = sum(
      count[
        clonality == "clonal"
      ]
    ),
    
    total_LOH_count = sum(count),
    
    .groups = "drop"
  )

excluded_none_summary <- sample_clonality_summary %>%
  dplyr::filter(
    clonality == "none"
  ) %>%
  dplyr::summarise(
    excluded_none_fragments = sum(
      count,
      na.rm = TRUE
    )
  )

figure_19_validation <- tibble::tibble(
  number_of_samples =
    dplyr::n_distinct(
      figure_19_plot_data$sample
    ),
  
  number_of_clonality_classes =
    dplyr::n_distinct(
      figure_19_plot_data$clonality
    ),
  
  clonal_LOH_fragment_count =
    sum(
      figure_19_plot_data$count[
        figure_19_plot_data$clonality ==
          "clonal"
      ]
    ),
  
  subclonal_LOH_fragment_count =
    sum(
      figure_19_plot_data$count[
        figure_19_plot_data$clonality ==
          "subclonal"
      ]
    ),
  
  total_plotted_LOH_fragments =
    sum(
      figure_19_plot_data$count
    ),
  
  excluded_none_fragments =
    excluded_none_summary$
    excluded_none_fragments,
  
  minimum_sample_LOH_total =
    min(
      sample_count_summary$total_LOH_count
    ),
  
  maximum_sample_LOH_total =
    max(
      sample_count_summary$total_LOH_count
    ),
  
  negative_fragment_counts =
    sum(
      figure_19_plot_data$count < 0
    )
)

print(
  figure_19_validation,
  n = Inf,
  width = Inf
)

print(
  sample_count_summary,
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
    "One or more Figure 19 output files were not created."
  )
}

# -----------------------------------------------------------------------------
# Figure legend:
#
# Figure 19. Area chart showing the total number of clonal and subclonal LOH
# fragments within each tumour-region sample from the combined LUAD and LUSC
# cohort.
#
# Green represents subclonal LOH and salmon-red represents clonal LOH. Rare
# records classified as "none" because of strict interval-overlap edge cases
# are excluded from this biological comparison.
# -----------------------------------------------------------------------------