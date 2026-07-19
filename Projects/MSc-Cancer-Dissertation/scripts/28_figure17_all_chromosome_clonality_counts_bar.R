# -----------------------------------------------------------------------------
# Original MSc Cancer Dissertation analysis (2021)
# Refactored for reproducibility and GitHub portfolio (2026)
# -----------------------------------------------------------------------------

# 28_figure17_all_chromosome_clonality_counts_bar.R
# Purpose: reproduce Dissertation Figure 17 as a stacked bar chart showing
# the number of clonal LOH, subclonal LOH, and non-LOH fragments on each
# chromosome in the complete ASCAT segment dataset

library(dplyr)
library(tidyr)
library(ggplot2)

# -----------------------------------------------------------------------------
# Input file
# -----------------------------------------------------------------------------

chromosome_summary_file <-
  "data/processed/all_chromosome_clonality_summary.rds"

if (!file.exists(chromosome_summary_file)) {
  stop(
    "Required input file was not found: ",
    chromosome_summary_file
  )
}

# -----------------------------------------------------------------------------
# Load validated chromosome-clonality summary
# -----------------------------------------------------------------------------

all_chromosome_clonality_summary <- readRDS(
  chromosome_summary_file
)

required_columns <- c(
  "chr",
  "clonality",
  "count"
)

missing_columns <- setdiff(
  required_columns,
  names(all_chromosome_clonality_summary)
)

if (length(missing_columns) > 0) {
  stop(
    "The chromosome-clonality summary is missing the following column(s): ",
    paste(missing_columns, collapse = ", ")
  )
}

# -----------------------------------------------------------------------------
# Validate source data
# -----------------------------------------------------------------------------

if (any(is.na(all_chromosome_clonality_summary$chr))) {
  stop(
    "Missing chromosome values were detected."
  )
}

if (any(is.na(all_chromosome_clonality_summary$count))) {
  stop(
    "Missing fragment counts were detected."
  )
}

if (any(all_chromosome_clonality_summary$count < 0)) {
  stop(
    "Negative fragment counts were detected."
  )
}

unexpected_clonality <- setdiff(
  unique(
    as.character(
      all_chromosome_clonality_summary$clonality
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
# Prepare plotting data
# -----------------------------------------------------------------------------

figure_17_plot_data <-
  all_chromosome_clonality_summary %>%
  dplyr::select(
    chr,
    clonality,
    count
  ) %>%
  tidyr::complete(
    chr = 1:24,
    clonality = c(
      "subclonal",
      "none",
      "clonal"
    ),
    fill = list(
      count = 0
    )
  ) %>%
  dplyr::mutate(
    chr = as.integer(chr),
    
    # The factor order is combined with reverse stacking below to produce:
    # bottom: Subclonal LOH
    # middle: No LOH
    # top: Clonal LOH
    clonality = factor(
      clonality,
      levels = c(
        "subclonal",
        "none",
        "clonal"
      ),
      labels = c(
        "Subclonal LOH",
        "No LOH",
        "Clonal LOH"
      )
    )
  ) %>%
  dplyr::arrange(
    chr,
    clonality
  )

# -----------------------------------------------------------------------------
# Dissertation-style colours
# -----------------------------------------------------------------------------

clonality_colours <- c(
  "Subclonal LOH" = "#FFF200",
  "No LOH" = "#FF2A00",
  "Clonal LOH" = "#003CFF"
)

# -----------------------------------------------------------------------------
# Create stacked bar chart
# -----------------------------------------------------------------------------

figure_17 <- ggplot(
  figure_17_plot_data,
  aes(
    x = chr,
    y = count,
    fill = clonality
  )
) +
  geom_col(
    position = position_stack(
      reverse = TRUE
    ),
    width = 0.88,
    colour = "white",
    linewidth = 0.25
  ) +
  scale_fill_manual(
    values = clonality_colours,
    breaks = c(
      "Clonal LOH",
      "No LOH",
      "Subclonal LOH"
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
      1000,
      2000,
      3000,
      4000
    ),
    expand = expansion(
      mult = c(0, 0.05)
    )
  ) +
  labs(
    title = paste(
      "number of fragments in each chromosome",
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
  "figures/Figure17_all_chromosome_clonality_counts_bar.png"

figure_pdf <-
  "figures/Figure17_all_chromosome_clonality_counts_bar.pdf"

ggsave(
  filename = figure_png,
  plot = figure_17,
  width = 10,
  height = 7,
  units = "in",
  dpi = 300,
  bg = "white"
)

ggsave(
  filename = figure_pdf,
  plot = figure_17,
  width = 10,
  height = 7,
  units = "in",
  bg = "white"
)

# -----------------------------------------------------------------------------
# Validation
# -----------------------------------------------------------------------------

chromosome_count_summary <-
  figure_17_plot_data %>%
  dplyr::group_by(chr) %>%
  dplyr::summarise(
    subclonal_LOH_count = sum(
      count[
        clonality == "Subclonal LOH"
      ]
    ),
    
    no_LOH_count = sum(
      count[
        clonality == "No LOH"
      ]
    ),
    
    clonal_LOH_count = sum(
      count[
        clonality == "Clonal LOH"
      ]
    ),
    
    total_count = sum(count),
    
    .groups = "drop"
  )

clonality_count_summary <-
  figure_17_plot_data %>%
  dplyr::group_by(clonality) %>%
  dplyr::summarise(
    fragment_count = sum(count),
    .groups = "drop"
  )

figure_17_validation <- tibble::tibble(
  number_of_chromosomes =
    dplyr::n_distinct(
      figure_17_plot_data$chr
    ),
  
  number_of_clonality_classes =
    dplyr::n_distinct(
      figure_17_plot_data$clonality
    ),
  
  clonal_LOH_fragment_count =
    sum(
      figure_17_plot_data$count[
        figure_17_plot_data$clonality ==
          "Clonal LOH"
      ]
    ),
  
  subclonal_LOH_fragment_count =
    sum(
      figure_17_plot_data$count[
        figure_17_plot_data$clonality ==
          "Subclonal LOH"
      ]
    ),
  
  no_LOH_fragment_count =
    sum(
      figure_17_plot_data$count[
        figure_17_plot_data$clonality ==
          "No LOH"
      ]
    ),
  
  total_ASCAT_fragments =
    sum(
      figure_17_plot_data$count
    ),
  
  minimum_chromosome_total =
    min(
      chromosome_count_summary$total_count
    ),
  
  maximum_chromosome_total =
    max(
      chromosome_count_summary$total_count
    ),
  
  negative_fragment_counts =
    sum(
      figure_17_plot_data$count < 0
    )
)

print(
  figure_17_validation,
  n = Inf,
  width = Inf
)

print(
  chromosome_count_summary,
  n = Inf
)

print(
  clonality_count_summary,
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
    "One or more Figure 17 output files were not created."
  )
}

# -----------------------------------------------------------------------------
# Figure legend:
#
# Figure 17. Stacked bar chart showing the number of fragments on each
# chromosome classified as clonal LOH, subclonal LOH, or no LOH across the
# complete tracerx.ascat.seg dataset for the combined LUAD and LUSC cohort.
#
# Yellow represents subclonal LOH, red represents fragments classified as
# having no LOH, and blue represents clonal LOH.
#
# The original MSc analysis used the clonality value "none" for fragments
# without classified clonal or subclonal LOH. This category is displayed as
# "No LOH" for clarity.
# -----------------------------------------------------------------------------