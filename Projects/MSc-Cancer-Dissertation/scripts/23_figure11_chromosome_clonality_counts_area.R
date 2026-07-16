# -----------------------------------------------------------------------------
# Original MSc Cancer Dissertation analysis (2021)
# Refactored for reproducibility and GitHub portfolio (2026)
# -----------------------------------------------------------------------------

# 23_figure11_chromosome_clonality_counts_area.R
# Purpose: reproduce Dissertation Figure 11, showing the number of clonal
# and subclonal LOH fragments across chromosomes

library(dplyr)
library(tidyr)
library(ggplot2)

# -----------------------------------------------------------------------------
# Input file
# -----------------------------------------------------------------------------

chromosome_summary_file <-
  "data/processed/loh_chromosome_clonality_summary.rds"

if (!file.exists(chromosome_summary_file)) {
  stop(
    "Required input file was not found: ",
    chromosome_summary_file
  )
}

# -----------------------------------------------------------------------------
# Load validated chromosome-clonality summary
# -----------------------------------------------------------------------------

chromosome_clonality_summary <- readRDS(
  chromosome_summary_file
)

required_columns <- c(
  "chr",
  "clonality",
  "count"
)

if (!all(required_columns %in% names(chromosome_clonality_summary))) {
  stop(
    "The chromosome-clonality summary is missing one or more ",
    "required columns: ",
    paste(required_columns, collapse = ", ")
  )
}

# -----------------------------------------------------------------------------
# Prepare plotting data
# -----------------------------------------------------------------------------

# Figure 11 compares the biological clonality classes only.
# Rare "none" records are excluded because they arise from strict
# interval-overlap edge cases rather than a biological clonality category.
figure_11_plot_data <- chromosome_clonality_summary %>%
  dplyr::filter(
    clonality %in% c(
      "clonal",
      "subclonal"
    )
  ) %>%
  dplyr::select(
    chr,
    clonality,
    count
  ) %>%
  tidyr::complete(
    chr = 1:24,
    clonality = c(
      "subclonal",
      "clonal"
    ),
    fill = list(
      count = 0
    )
  ) %>%
  dplyr::mutate(
    chr = as.integer(chr),
    
    # This factor order places subclonal in green at the bottom
    # and clonal in red above it.
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

# Confirm valid values
if (any(is.na(figure_11_plot_data$count))) {
  stop("Missing fragment counts were detected.")
}

if (any(figure_11_plot_data$count < 0)) {
  stop("Negative fragment counts were detected.")
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

figure_11 <- ggplot(
  figure_11_plot_data,
  aes(
    x = chr,
    y = count,
    fill = clonality
  )
) +
  geom_area(
    position = "stack",
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
    name = "clonality"
  ) +
  scale_x_continuous(
    breaks = c(
      0,
      5,
      10,
      15,
      20,
      25
    ),
    limits = c(0, 25),
    expand = expansion(
      mult = c(0, 0)
    )
  ) +
  scale_y_continuous(
    breaks = c(
      0,
      400,
      800,
      1200
    ),
    expand = expansion(
      mult = c(0, 0.04)
    )
  ) +
  labs(
    title = "number of LOH fragments with a particular clonality in each chromosome",
    x = "chr",
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
    
    axis.text.x = element_text(
      size = 10
    ),
    
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
  "figures/Figure11_chromosome_clonality_counts_area.png"

figure_pdf <-
  "figures/Figure11_chromosome_clonality_counts_area.pdf"

ggsave(
  filename = figure_png,
  plot = figure_11,
  width = 10,
  height = 7,
  units = "in",
  dpi = 300,
  bg = "white"
)

ggsave(
  filename = figure_pdf,
  plot = figure_11,
  width = 10,
  height = 7,
  units = "in",
  bg = "white"
)

# -----------------------------------------------------------------------------
# Validation
# -----------------------------------------------------------------------------

chromosome_count_summary <- figure_11_plot_data %>%
  dplyr::group_by(chr) %>%
  dplyr::summarise(
    clonal_count = sum(
      count[clonality == "clonal"]
    ),
    
    subclonal_count = sum(
      count[clonality == "subclonal"]
    ),
    
    total_count = sum(count),
    
    .groups = "drop"
  )

figure_11_validation <- tibble::tibble(
  number_of_chromosomes =
    dplyr::n_distinct(
      figure_11_plot_data$chr
    ),
  
  number_of_clonality_classes =
    dplyr::n_distinct(
      figure_11_plot_data$clonality
    ),
  
  clonal_fragment_count =
    sum(
      figure_11_plot_data$count[
        figure_11_plot_data$clonality == "clonal"
      ]
    ),
  
  subclonal_fragment_count =
    sum(
      figure_11_plot_data$count[
        figure_11_plot_data$clonality == "subclonal"
      ]
    ),
  
  total_plotted_LOH_fragments =
    sum(
      figure_11_plot_data$count
    ),
  
  excluded_none_fragments =
    sum(
      chromosome_clonality_summary$count[
        chromosome_clonality_summary$clonality == "none"
      ],
      na.rm = TRUE
    ),
  
  maximum_chromosome_total =
    max(
      chromosome_count_summary$total_count
    ),
  
  minimum_chromosome_total =
    min(
      chromosome_count_summary$total_count
    )
)

print(
  figure_11_validation,
  n = Inf,
  width = Inf
)

print(
  chromosome_count_summary,
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
    "One or more Figure 11 output files were not created."
  )
}

# -----------------------------------------------------------------------------
# Figure legend:
#
# Figure 11. Area chart showing the total number of LOH fragments in each
# chromosome across the combined LUAD and LUSC cohort, together with the
# contribution of each clonality class.
#
# The green area represents subclonal LOH fragments and the red area represents
# clonal LOH fragments. Rare records classified as "none" because of strict
# interval-overlap edge cases are excluded from the biological comparison.
# -----------------------------------------------------------------------------