# -----------------------------------------------------------------------------
# Original MSc Cancer Dissertation analysis (2021)
# Refactored for reproducibility and GitHub portfolio (2026)
# -----------------------------------------------------------------------------

# 26_figure15_all_chromosome_clonality_proportions.R
# Purpose: reproduce Dissertation Figure 15 as a proportional stacked area
# chart showing clonal LOH, subclonal LOH, and non-LOH fragments across
# chromosomes in the complete ASCAT segment dataset

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

figure_15_plot_data <-
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
  dplyr::group_by(chr) %>%
  dplyr::mutate(
    chromosome_total_fragments =
      sum(count),
    
    chromosome_proportion =
      count / chromosome_total_fragments
  ) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(
    chr = as.integer(chr),
    
    # This ordering, together with reverse stacking below, produces:
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

if (
  any(
    figure_15_plot_data$
    chromosome_total_fragments <= 0
  )
) {
  stop(
    "One or more chromosomes have no fragments."
  )
}

# -----------------------------------------------------------------------------
# Dissertation-style colours
# -----------------------------------------------------------------------------

clonality_colours <- c(
  "Subclonal LOH" = "#619CFF",
  "No LOH" = "#00BA38",
  "Clonal LOH" = "#F8766D"
)

# -----------------------------------------------------------------------------
# Create proportional stacked area chart
# -----------------------------------------------------------------------------

figure_15 <- ggplot(
  figure_15_plot_data,
  aes(
    x = chr,
    y = chromosome_proportion,
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
    title = paste(
      "proportion of each clonality type in each chromosome",
      "in whole data studied"
    ),
    x = "chromosome",
    y = "proportion of the chromosome"
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
  paste0(
    "figures/",
    "Figure15_all_chromosome_",
    "clonality_proportions.png"
  )

figure_pdf <-
  paste0(
    "figures/",
    "Figure15_all_chromosome_",
    "clonality_proportions.pdf"
  )

ggsave(
  filename = figure_png,
  plot = figure_15,
  width = 10,
  height = 7,
  units = "in",
  dpi = 300,
  bg = "white"
)

ggsave(
  filename = figure_pdf,
  plot = figure_15,
  width = 10,
  height = 7,
  units = "in",
  bg = "white"
)

# -----------------------------------------------------------------------------
# Validation
# -----------------------------------------------------------------------------

chromosome_proportion_check <-
  figure_15_plot_data %>%
  dplyr::group_by(chr) %>%
  dplyr::summarise(
    total_proportion = sum(
      chromosome_proportion
    ),
    
    chromosome_total_fragments =
      dplyr::first(
        chromosome_total_fragments
      ),
    
    .groups = "drop"
  )

clonality_count_summary <-
  figure_15_plot_data %>%
  dplyr::group_by(clonality) %>%
  dplyr::summarise(
    fragment_count = sum(count),
    .groups = "drop"
  )

figure_15_validation <- tibble::tibble(
  number_of_chromosomes =
    dplyr::n_distinct(
      figure_15_plot_data$chr
    ),
  
  number_of_clonality_classes =
    dplyr::n_distinct(
      figure_15_plot_data$clonality
    ),
  
  clonal_LOH_fragment_count =
    sum(
      figure_15_plot_data$count[
        figure_15_plot_data$clonality ==
          "Clonal LOH"
      ]
    ),
  
  subclonal_LOH_fragment_count =
    sum(
      figure_15_plot_data$count[
        figure_15_plot_data$clonality ==
          "Subclonal LOH"
      ]
    ),
  
  no_LOH_fragment_count =
    sum(
      figure_15_plot_data$count[
        figure_15_plot_data$clonality ==
          "No LOH"
      ]
    ),
  
  total_ASCAt_fragments =
    sum(
      figure_15_plot_data$count
    ),
  
  all_chromosome_proportions_equal_one =
    all(
      abs(
        chromosome_proportion_check$
          total_proportion - 1
      ) < 1e-12
    ),
  
  minimum_chromosome_total =
    min(
      chromosome_proportion_check$
        chromosome_total_fragments
    ),
  
  maximum_chromosome_total =
    max(
      chromosome_proportion_check$
        chromosome_total_fragments
    )
)

print(
  figure_15_validation,
  n = Inf,
  width = Inf
)

print(
  chromosome_proportion_check,
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
    "One or more Figure 15 output files were not created."
  )
}

# -----------------------------------------------------------------------------
# Figure legend:
#
# Figure 15. Area chart showing the proportion of fragments on each chromosome
# classified as clonal LOH, subclonal LOH, or no LOH across the complete
# tracerx.ascat.seg dataset for the combined LUAD and LUSC cohort.
#
# Blue represents subclonal LOH, green represents fragments classified as
# having no LOH, and salmon-red represents clonal LOH. For each chromosome,
# the three proportions sum to 1.
#
# The original MSc analysis used the clonality value "none" to represent the
# absence of classified clonal or subclonal LOH. This category is displayed
# here as "No LOH" to make the legend easier to interpret.
# -----------------------------------------------------------------------------