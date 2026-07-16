# -----------------------------------------------------------------------------
# Original MSc Cancer Dissertation analysis (2021)
# Refactored for reproducibility and GitHub portfolio (2026)
# -----------------------------------------------------------------------------

# 20_figure8_sample_LOH_percentage.R
# Purpose: reproduce Dissertation Figure 8, showing tumour-region samples
# ordered by their percentage contribution to all LOH fragments, with
# category composition displayed as stacked bars

library(dplyr)
library(ggplot2)

# -----------------------------------------------------------------------------
# Input file
# -----------------------------------------------------------------------------

summary_file <-
  "data/processed/sample_category_summary.rds"

if (!file.exists(summary_file)) {
  stop(
    "Required input file was not found: ",
    summary_file
  )
}

# -----------------------------------------------------------------------------
# Load validated sample-category summary
# -----------------------------------------------------------------------------

sample_category_summary <- readRDS(
  summary_file
)

required_columns <- c(
  "sample",
  "category",
  "count"
)

if (!all(required_columns %in% names(sample_category_summary))) {
  stop(
    "The input dataset is missing one or more required columns: ",
    paste(required_columns, collapse = ", ")
  )
}

# -----------------------------------------------------------------------------
# Prepare data
# -----------------------------------------------------------------------------

total_LOH_fragments <- sum(
  sample_category_summary$count,
  na.rm = TRUE
)

if (total_LOH_fragments <= 0) {
  stop(
    "The total number of LOH fragments must be greater than zero."
  )
}

# Calculate total LOH fragments per sample
sample_totals <- sample_category_summary %>%
  dplyr::group_by(sample) %>%
  dplyr::summarise(
    sample_total_LOH = sum(
      count,
      na.rm = TRUE
    ),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    sample_percentage_of_all_LOH =
      sample_total_LOH / total_LOH_fragments
  ) %>%
  dplyr::arrange(
    dplyr::desc(sample_percentage_of_all_LOH),
    sample
  )

# Preserve descending sample order in the plot
ordered_samples <- sample_totals$sample

figure_8_plot_data <- sample_category_summary %>%
  dplyr::left_join(
    sample_totals,
    by = "sample"
  ) %>%
  dplyr::mutate(
    category = factor(
      category,
      levels = as.character(1:8),
      labels = paste0(
        "category",
        1:8
      )
    ),
    
    sample = factor(
      sample,
      levels = ordered_samples
    ),
    
    # Contribution of this sample-category combination
    # to the complete set of LOH fragments
    percentage_of_all_LOH =
      count / total_LOH_fragments
  ) %>%
  dplyr::arrange(
    sample,
    category
  )

# -----------------------------------------------------------------------------
# Original dissertation-style category colours
# -----------------------------------------------------------------------------

category_colours <- c(
  category1 = "#00E600",
  category2 = "#FF1E00",
  category3 = "#FFF200",
  category4 = "#FFA500",
  category5 = "#0047FF",
  category6 = "#F4A6B8",
  category7 = "#A52A2A",
  category8 = "#000000"
)

# -----------------------------------------------------------------------------
# Create stacked bar chart
# -----------------------------------------------------------------------------

figure_8 <- ggplot(
  figure_8_plot_data,
  aes(
    x = sample,
    y = percentage_of_all_LOH,
    fill = category
  )
) +
  geom_col(
    width = 1,
    colour = NA
  ) +
  scale_fill_manual(
    values = category_colours,
    drop = FALSE,
    name = "category"
  ) +
  scale_y_continuous(
    labels = scales::label_percent(
      accuracy = 0.01
    ),
    breaks = c(
      0,
      0.0025,
      0.0050,
      0.0075
    ),
    expand = expansion(
      mult = c(0, 0.04)
    )
  ) +
  labs(
    title = "descending order of number of categories of LOH samples",
    x = "sample",
    y = "percentage of the total LOH samples"
  ) +
  theme_minimal(
    base_size = 12
  ) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    
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
    
    # There are 434 samples, so individual labels would be unreadable
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
  "figures/Figure8_sample_LOH_percentage.png"

figure_pdf <-
  "figures/Figure8_sample_LOH_percentage.pdf"

ggsave(
  filename = figure_png,
  plot = figure_8,
  width = 11,
  height = 7,
  units = "in",
  dpi = 300,
  bg = "white"
)

ggsave(
  filename = figure_pdf,
  plot = figure_8,
  width = 11,
  height = 7,
  units = "in",
  bg = "white"
)

# -----------------------------------------------------------------------------
# Validation
# -----------------------------------------------------------------------------

sample_order_validation <- sample_totals %>%
  dplyr::mutate(
    descending_order_valid =
      sample_percentage_of_all_LOH ==
      cummin(sample_percentage_of_all_LOH)
  )

figure_8_validation <- tibble::tibble(
  number_of_samples =
    dplyr::n_distinct(
      figure_8_plot_data$sample
    ),
  
  number_of_categories =
    dplyr::n_distinct(
      figure_8_plot_data$category
    ),
  
  sample_category_combinations =
    nrow(figure_8_plot_data),
  
  total_LOH_fragments =
    total_LOH_fragments,
  
  sum_of_sample_percentages =
    sum(
      sample_totals$sample_percentage_of_all_LOH
    ),
  
  highest_sample_percentage =
    max(
      sample_totals$sample_percentage_of_all_LOH
    ),
  
  lowest_sample_percentage =
    min(
      sample_totals$sample_percentage_of_all_LOH
    ),
  
  descending_order_confirmed =
    all(
      diff(
        sample_totals$sample_percentage_of_all_LOH
      ) <= 0
    )
)

print(
  figure_8_validation,
  n = Inf
)

print(
  head(
    sample_totals,
    10
  )
)

print(
  tail(
    sample_totals,
    10
  )
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
    "One or more Figure 8 output files were not created."
  )
}

# -----------------------------------------------------------------------------
# Figure legend:
#
# Figure 8. Tumour-region samples arranged in descending order according to
# their percentage contribution to all LOH fragments.
#
# Each bar represents one tumour-region sample. The total bar height represents
# the proportion of all LOH fragments contributed by that sample, and the
# coloured sections represent the eight LOH categories. The sample with the
# greatest percentage contribution appears at the far left, while the sample
# with the smallest contribution appears at the far right.
# -----------------------------------------------------------------------------
