# -----------------------------------------------------------------------------
# Original MSc Cancer Dissertation analysis (2021)
# Refactored for reproducibility and GitHub portfolio (2026)
# -----------------------------------------------------------------------------

# 24_figure12_chromosome_clonality_counts_bar.R
# Purpose: reproduce Dissertation Figure 12 as a stacked bar chart showing
# the numbers of clonal and subclonal LOH fragments on each chromosome

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

missing_columns <- setdiff(
  required_columns,
  names(chromosome_clonality_summary)
)

if (length(missing_columns) > 0) {
  stop(
    "The chromosome-clonality summary is missing the following column(s): ",
    paste(missing_columns, collapse = ", ")
  )
}

# -----------------------------------------------------------------------------
# Prepare plotting data
# -----------------------------------------------------------------------------

# Figure 12 compares the biological clonal and subclonal LOH classes.
# The 10 rare records classified as "none" are excluded because they arose
# from strict interval-overlap edge cases in the reconstructed MSc workflow.

figure_12_plot_data <- chromosome_clonality_summary %>%
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
    
    # Subclonal is placed at the bottom of each bar,
    # with clonal stacked above it.
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

# -----------------------------------------------------------------------------
# Validate plotting values
# -----------------------------------------------------------------------------

if (any(is.na(figure_12_plot_data$count))) {
  stop(
    "Missing fragment counts were detected in the plotting data."
  )
}

if (any(figure_12_plot_data$count < 0)) {
  stop(
    "Negative fragment counts were detected in the plotting data."
  )
}

# -----------------------------------------------------------------------------
# Dissertation-style colours
# -----------------------------------------------------------------------------

clonality_colours <- c(
  subclonal = "#F0442D",
  clonal = "#1749F5"
)

# -----------------------------------------------------------------------------
# Create stacked bar chart
# -----------------------------------------------------------------------------

figure_12 <- ggplot(
  figure_12_plot_data,
  aes(
    x = chr,
    y = count,
    fill = clonality
  )
) +
  geom_col(
    width = 0.88,
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
    title = paste(
      "number of LOH fragments in each chromosome",
      "and clonality proportion"
    ),
    x = "chromosome",
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
  "figures/Figure12_chromosome_clonality_counts_bar.png"

figure_pdf <-
  "figures/Figure12_chromosome_clonality_counts_bar.pdf"

ggsave(
  filename = figure_png,
  plot = figure_12,
  width = 10,
  height = 7,
  units = "in",
  dpi = 300,
  bg = "white"
)

ggsave(
  filename = figure_pdf,
  plot = figure_12,
  width = 10,
  height = 7,
  units = "in",
  bg = "white"
)

# -----------------------------------------------------------------------------
# Validation
# -----------------------------------------------------------------------------

chromosome_count_summary <- figure_12_plot_data %>%
  dplyr::group_by(chr) %>%
  dplyr::summarise(
    subclonal_count = sum(
      count[clonality == "subclonal"]
    ),
    
    clonal_count = sum(
      count[clonality == "clonal"]
    ),
    
    total_count = sum(count),
    
    .groups = "drop"
  )

figure_12_validation <- tibble::tibble(
  number_of_chromosomes =
    dplyr::n_distinct(
      figure_12_plot_data$chr
    ),
  
  number_of_clonality_classes =
    dplyr::n_distinct(
      figure_12_plot_data$clonality
    ),
  
  clonal_fragment_count =
    sum(
      figure_12_plot_data$count[
        figure_12_plot_data$clonality == "clonal"
      ]
    ),
  
  subclonal_fragment_count =
    sum(
      figure_12_plot_data$count[
        figure_12_plot_data$clonality == "subclonal"
      ]
    ),
  
  total_plotted_LOH_fragments =
    sum(
      figure_12_plot_data$count
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
  figure_12_validation,
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
    "One or more Figure 12 output files were not created."
  )
}

# -----------------------------------------------------------------------------
# Figure legend:
#
# Figure 12. Stacked bar chart showing the total number of LOH fragments on
# each chromosome across the combined LUAD and LUSC cohort, including the
# contribution of each clonality class to every chromosome.
#
# Red represents subclonal LOH fragments and blue represents clonal LOH
# fragments. Rare records classified as "none" because of strict
# interval-overlap edge cases are excluded from this biological comparison.
# -----------------------------------------------------------------------------