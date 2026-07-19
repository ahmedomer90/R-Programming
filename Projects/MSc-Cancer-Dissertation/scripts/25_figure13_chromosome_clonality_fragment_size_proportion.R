# -----------------------------------------------------------------------------
# Original MSc Cancer Dissertation analysis (2021)
# Refactored for reproducibility and GitHub portfolio (2026)
# -----------------------------------------------------------------------------

# 25_figure13_chromosome_clonality_fragment_size_proportion.R
# Purpose: reproduce Dissertation Figure 13 as a proportional stacked area
# chart showing the total LOH fragment-size contribution of clonal and
# subclonal LOH on each chromosome

library(dplyr)
library(tidyr)
library(ggplot2)

# -----------------------------------------------------------------------------
# Input file
# -----------------------------------------------------------------------------

loh_fragment_size_file <-
  "data/processed/loh_segments_fragment_sizes.rds"

if (!file.exists(loh_fragment_size_file)) {
  stop(
    "Required input file was not found: ",
    loh_fragment_size_file
  )
}

# -----------------------------------------------------------------------------
# Load validated LOH fragment-size data
# -----------------------------------------------------------------------------

loh_segments <- readRDS(
  loh_fragment_size_file
)

required_columns <- c(
  "chr",
  "clonality",
  "fragment_size_bp"
)

missing_columns <- setdiff(
  required_columns,
  names(loh_segments)
)

if (length(missing_columns) > 0) {
  stop(
    "The LOH fragment-size dataset is missing the following column(s): ",
    paste(missing_columns, collapse = ", ")
  )
}

# -----------------------------------------------------------------------------
# Validate source data
# -----------------------------------------------------------------------------

if (any(is.na(loh_segments$chr))) {
  stop(
    "Missing chromosome values were detected."
  )
}

if (any(is.na(loh_segments$fragment_size_bp))) {
  stop(
    "Missing fragment-size values were detected."
  )
}

if (any(loh_segments$fragment_size_bp < 0)) {
  stop(
    "Negative fragment sizes were detected."
  )
}

# -----------------------------------------------------------------------------
# Summarise total fragment size by chromosome and clonality
# -----------------------------------------------------------------------------

chromosome_clonality_fragment_size <-
  loh_segments %>%
  dplyr::filter(
    clonality %in% c(
      "clonal",
      "subclonal"
    )
  ) %>%
  dplyr::group_by(
    chr,
    clonality
  ) %>%
  dplyr::summarise(
    total_fragment_size_bp = sum(
      fragment_size_bp,
      na.rm = TRUE
    ),
    fragment_count = dplyr::n(),
    .groups = "drop"
  ) %>%
  tidyr::complete(
    chr = 1:24,
    clonality = c(
      "subclonal",
      "clonal"
    ),
    fill = list(
      total_fragment_size_bp = 0,
      fragment_count = 0
    )
  )

# -----------------------------------------------------------------------------
# Calculate clonality proportions within each chromosome
# -----------------------------------------------------------------------------

figure_13_plot_data <-
  chromosome_clonality_fragment_size %>%
  dplyr::group_by(chr) %>%
  dplyr::mutate(
    chromosome_total_fragment_size_bp =
      sum(total_fragment_size_bp),
    
    fragment_size_proportion =
      total_fragment_size_bp /
      chromosome_total_fragment_size_bp
  ) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(
    chr = as.integer(chr),
    
    # Subclonal is drawn at the bottom in blue;
    # clonal is stacked above it in salmon-red.
    clonality = factor(
      clonality,
      levels = c(
        "subclonal",
        "clonal"
      )
    )
  ) %>%
  dplyr::arrange(
    chr,
    clonality
  )

# Confirm that each chromosome has a positive total fragment size
if (
  any(
    figure_13_plot_data$
    chromosome_total_fragment_size_bp <= 0
  )
) {
  stop(
    "One or more chromosomes have no clonal or subclonal ",
    "LOH fragment size."
  )
}

# -----------------------------------------------------------------------------
# Dissertation-style colours
# -----------------------------------------------------------------------------

clonality_colours <- c(
  subclonal = "#619CFF",
  clonal = "#F8766D"
)

# -----------------------------------------------------------------------------
# Create proportional stacked area chart
# -----------------------------------------------------------------------------

figure_13 <- ggplot(
  figure_13_plot_data,
  aes(
    x = chr,
    y = fragment_size_proportion,
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
    breaks = c(0, 5, 10, 15, 20, 25),
    limits = c(0, 25),
    expand = expansion(mult = c(0, 0))
  ) +
  scale_y_continuous(
    breaks = c(0, 0.25, 0.50, 0.75, 1),
    limits = c(0, 1),
    expand = expansion(mult = c(0, 0)),
    labels = scales::label_number(
      accuracy = 0.01
    )
  ) +
  labs(
    title = paste(
      "clonality proportion of the total fragment sizes",
      "with LOH in chromosomes"
    ),
    x = "chromosome",
    y = "proportion of total fragments size with LOH"
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
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12),
    axis.text.x = element_text(size = 10),
    axis.text.y = element_text(size = 10),
    legend.title = element_text(size = 11),
    legend.text = element_text(size = 10),
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
  paste0(
    "figures/",
    "Figure13_chromosome_clonality_",
    "fragment_size_proportion.png"
  )

figure_pdf <-
  paste0(
    "figures/",
    "Figure13_chromosome_clonality_",
    "fragment_size_proportion.pdf"
  )

ggsave(
  filename = figure_png,
  plot = figure_13,
  width = 10,
  height = 7,
  units = "in",
  dpi = 300,
  bg = "white"
)

ggsave(
  filename = figure_pdf,
  plot = figure_13,
  width = 10,
  height = 7,
  units = "in",
  bg = "white"
)

# -----------------------------------------------------------------------------
# Validation
# -----------------------------------------------------------------------------

chromosome_proportion_check <-
  figure_13_plot_data %>%
  dplyr::group_by(chr) %>%
  dplyr::summarise(
    total_proportion = sum(
      fragment_size_proportion
    ),
    
    chromosome_total_fragment_size_bp =
      dplyr::first(
        chromosome_total_fragment_size_bp
      ),
    
    .groups = "drop"
  )

excluded_none_summary <-
  loh_segments %>%
  dplyr::filter(
    clonality == "none"
  ) %>%
  dplyr::summarise(
    excluded_none_fragments =
      dplyr::n(),
    
    excluded_none_fragment_size_bp =
      sum(
        fragment_size_bp,
        na.rm = TRUE
      )
  )

figure_13_validation <- tibble::tibble(
  number_of_chromosomes =
    dplyr::n_distinct(
      figure_13_plot_data$chr
    ),
  
  number_of_clonality_classes =
    dplyr::n_distinct(
      figure_13_plot_data$clonality
    ),
  
  clonal_fragment_count =
    sum(
      figure_13_plot_data$fragment_count[
        figure_13_plot_data$clonality ==
          "clonal"
      ]
    ),
  
  subclonal_fragment_count =
    sum(
      figure_13_plot_data$fragment_count[
        figure_13_plot_data$clonality ==
          "subclonal"
      ]
    ),
  
  clonal_total_fragment_size_bp =
    sum(
      figure_13_plot_data$
        total_fragment_size_bp[
          figure_13_plot_data$clonality ==
            "clonal"
        ]
    ),
  
  subclonal_total_fragment_size_bp =
    sum(
      figure_13_plot_data$
        total_fragment_size_bp[
          figure_13_plot_data$clonality ==
            "subclonal"
        ]
    ),
  
  total_plotted_fragment_size_bp =
    sum(
      figure_13_plot_data$
        total_fragment_size_bp
    ),
  
  excluded_none_fragments =
    excluded_none_summary$
    excluded_none_fragments,
  
  excluded_none_fragment_size_bp =
    excluded_none_summary$
    excluded_none_fragment_size_bp,
  
  all_chromosome_proportions_equal_one =
    all(
      abs(
        chromosome_proportion_check$
          total_proportion - 1
      ) < 1e-12
    )
)

print(
  figure_13_validation,
  n = Inf,
  width = Inf
)

print(
  chromosome_proportion_check,
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
    "One or more Figure 13 output files were not created."
  )
}

# -----------------------------------------------------------------------------
# Figure legend:
#
# Figure 13. Area chart showing the proportion of total LOH fragment size
# contributed by clonal and subclonal LOH on each chromosome across the
# combined LUAD and LUSC cohort.
#
# Blue represents the proportion of total LOH fragment size classified as
# subclonal, while salmon-red represents the proportion classified as clonal.
# Each chromosome is scaled so that the two proportions sum to 1. Rare records
# classified as "none" because of strict interval-overlap edge cases are
# excluded from the biological comparison.
# -----------------------------------------------------------------------------