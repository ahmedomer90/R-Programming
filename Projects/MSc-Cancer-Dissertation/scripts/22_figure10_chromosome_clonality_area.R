# -----------------------------------------------------------------------------
# Original MSc Cancer Dissertation analysis (2021)
# Refactored for reproducibility and GitHub portfolio (2026)
# -----------------------------------------------------------------------------

# 22_figure10_chromosome_clonality_area.R
# Purpose: reproduce Dissertation Figure 10, showing the proportion of
# clonal and subclonal LOH fragments across chromosomes

library(dplyr)
library(tidyr)
library(ggplot2)
library(readr)

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
  "count",
  "proportion"
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

# The rare "none" records arose from strict interval-overlap edge cases.
# They are excluded here because Figure 10 compares biological clonality:
# clonal versus subclonal LOH.
figure_10_plot_data <- chromosome_clonality_summary %>%
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
  dplyr::group_by(chr) %>%
  dplyr::mutate(
    chromosome_total =
      sum(count),
    
    proportion =
      count / chromosome_total
  ) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(
    chr = as.integer(chr),
    
    # Order is important:
    # subclonal is drawn at the bottom and clonal at the top.
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

# Confirm that no chromosome has a zero total
if (any(figure_10_plot_data$chromosome_total <= 0)) {
  stop(
    "One or more chromosomes have no clonal or subclonal LOH fragments."
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
# Create stacked proportional area chart
# -----------------------------------------------------------------------------

figure_10 <- ggplot(
  figure_10_plot_data,
  aes(
    x = chr,
    y = proportion,
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
    title = "clonality of LOH in the genome and chromosome",
    x = "chromosome",
    y = "proportion of the genome with LOH"
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
  "figures/Figure10_chromosome_clonality_area.png"

figure_pdf <-
  "figures/Figure10_chromosome_clonality_area.pdf"

ggsave(
  filename = figure_png,
  plot = figure_10,
  width = 10,
  height = 7,
  units = "in",
  dpi = 300,
  bg = "white"
)

ggsave(
  filename = figure_pdf,
  plot = figure_10,
  width = 10,
  height = 7,
  units = "in",
  bg = "white"
)

# -----------------------------------------------------------------------------
# Validation
# -----------------------------------------------------------------------------

chromosome_proportion_check <- figure_10_plot_data %>%
  dplyr::group_by(chr) %>%
  dplyr::summarise(
    total_proportion =
      sum(proportion),
    
    chromosome_total =
      dplyr::first(chromosome_total),
    
    .groups = "drop"
  )

figure_10_validation <- tibble::tibble(
  number_of_chromosomes =
    dplyr::n_distinct(
      figure_10_plot_data$chr
    ),
  
  number_of_clonality_classes =
    dplyr::n_distinct(
      figure_10_plot_data$clonality
    ),
  
  clonal_fragment_count =
    sum(
      figure_10_plot_data$count[
        figure_10_plot_data$clonality == "clonal"
      ]
    ),
  
  subclonal_fragment_count =
    sum(
      figure_10_plot_data$count[
        figure_10_plot_data$clonality == "subclonal"
      ]
    ),
  
  total_plotted_LOH_fragments =
    sum(
      figure_10_plot_data$count
    ),
  
  excluded_none_fragments =
    sum(
      chromosome_clonality_summary$count[
        chromosome_clonality_summary$clonality == "none"
      ],
      na.rm = TRUE
    ),
  
  all_chromosome_proportions_equal_one =
    all(
      abs(
        chromosome_proportion_check$total_proportion - 1
      ) < 1e-12
    )
)

print(
  figure_10_validation,
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
    "One or more Figure 10 output files were not created."
  )
}

# -----------------------------------------------------------------------------
# Figure legend:
#
# Figure 10. Area chart representing the proportions of clonal and subclonal
# LOH fragments across chromosomes in the combined LUAD and LUSC dataset.
#
# For each chromosome, the total area is scaled to 1. The green section
# represents the proportion of subclonal LOH fragments, while the red section
# represents the proportion of clonal LOH fragments. Rare records classified
# as "none" because of strict interval-overlap edge cases are excluded from
# this biological comparison.
# -----------------------------------------------------------------------------