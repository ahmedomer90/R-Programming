# -----------------------------------------------------------------------------
# Original MSc Cancer Dissertation analysis (2021)
# Refactored for reproducibility and GitHub portfolio (2026)
# -----------------------------------------------------------------------------

# 29_figure18_sample_clonality_proportions_area.R
# Purpose: reproduce Dissertation Figure 18 as a proportional stacked area
# chart showing the proportions of clonal and subclonal LOH fragments within
# each tumour-region sample

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
    paste(missing_columns, collapse = ", ")
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
# Preserve the sample order stored in the validated summary
# -----------------------------------------------------------------------------

sample_order <- sample_clonality_summary %>%
  dplyr::distinct(sample) %>%
  dplyr::pull(sample) %>%
  as.character()

# -----------------------------------------------------------------------------
# Prepare plotting data
# -----------------------------------------------------------------------------

# Figure 18 compares the biological LOH clonality classes only.
# Rare records classified as "none" are excluded because they arose from
# strict interval-overlap edge cases rather than a biological clonality class.

figure_18_plot_data <- sample_clonality_summary %>%
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
  dplyr::group_by(sample) %>%
  dplyr::mutate(
    sample_total_LOH = sum(count),
    
    clonality_proportion =
      count / sample_total_LOH
  ) %>%
  dplyr::ungroup() %>%
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
# Confirm that each sample has at least one plotted LOH fragment
# -----------------------------------------------------------------------------

if (any(figure_18_plot_data$sample_total_LOH <= 0)) {
  stop(
    "One or more samples have no clonal or subclonal LOH fragments."
  )
}

if (any(is.na(figure_18_plot_data$clonality_proportion))) {
  stop(
    "Missing clonality proportions were produced."
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
# Create proportional stacked area chart
# -----------------------------------------------------------------------------

figure_18 <- ggplot(
  figure_18_plot_data,
  aes(
    x = sample_index,
    y = clonality_proportion,
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
      0.25,
      0.50,
      0.75,
      1
    ),
    limits = c(0, 1),
    expand = expansion(
      mult = c(0, 0)
    ),
    labels = scales::label_number(
      accuracy = 0.01
    )
  ) +
  labs(
    title = "clonality of LOH in the samples with LOH",
    x = "sample",
    y = "proportion of each LOH sample"
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
  "figures/Figure18_sample_clonality_proportions_area.png"

figure_pdf <-
  "figures/Figure18_sample_clonality_proportions_area.pdf"

ggsave(
  filename = figure_png,
  plot = figure_18,
  width = 11,
  height = 7,
  units = "in",
  dpi = 300,
  bg = "white"
)

ggsave(
  filename = figure_pdf,
  plot = figure_18,
  width = 11,
  height = 7,
  units = "in",
  bg = "white"
)

# -----------------------------------------------------------------------------
# Validation
# -----------------------------------------------------------------------------

sample_proportion_check <- figure_18_plot_data %>%
  dplyr::group_by(sample) %>%
  dplyr::summarise(
    total_proportion = sum(
      clonality_proportion
    ),
    
    sample_total_LOH =
      dplyr::first(
        sample_total_LOH
      ),
    
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

figure_18_validation <- tibble::tibble(
  number_of_samples =
    dplyr::n_distinct(
      figure_18_plot_data$sample
    ),
  
  number_of_clonality_classes =
    dplyr::n_distinct(
      figure_18_plot_data$clonality
    ),
  
  clonal_LOH_fragment_count =
    sum(
      figure_18_plot_data$count[
        figure_18_plot_data$clonality ==
          "clonal"
      ]
    ),
  
  subclonal_LOH_fragment_count =
    sum(
      figure_18_plot_data$count[
        figure_18_plot_data$clonality ==
          "subclonal"
      ]
    ),
  
  total_plotted_LOH_fragments =
    sum(
      figure_18_plot_data$count
    ),
  
  excluded_none_fragments =
    excluded_none_summary$
    excluded_none_fragments,
  
  minimum_sample_LOH_total =
    min(
      sample_proportion_check$sample_total_LOH
    ),
  
  maximum_sample_LOH_total =
    max(
      sample_proportion_check$sample_total_LOH
    ),
  
  all_sample_proportions_equal_one =
    all(
      abs(
        sample_proportion_check$
          total_proportion - 1
      ) < 1e-12
    ),
  
  negative_fragment_counts =
    sum(
      figure_18_plot_data$count < 0
    )
)

print(
  figure_18_validation,
  n = Inf,
  width = Inf
)

print(
  sample_proportion_check,
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
    "One or more Figure 18 output files were not created."
  )
}

# -----------------------------------------------------------------------------
# Figure legend:
#
# Figure 18. Area chart showing the proportions of clonal and subclonal LOH
# fragments within each tumour-region sample from the combined LUAD and LUSC
# cohort.
#
# Green represents subclonal LOH and salmon-red represents clonal LOH. Each
# sample is scaled so that its clonal and subclonal proportions sum to 1.
# Rare records classified as "none" because of strict interval-overlap edge
# cases are excluded from this biological comparison.
# -----------------------------------------------------------------------------