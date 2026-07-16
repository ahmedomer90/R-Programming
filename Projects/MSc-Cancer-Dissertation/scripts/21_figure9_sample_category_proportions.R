# -----------------------------------------------------------------------------
# Original MSc Cancer Dissertation analysis (2021)
# Refactored for reproducibility and GitHub portfolio (2026)
# -----------------------------------------------------------------------------

# 21_figure9_sample_category_proportions.R
# Purpose: reproduce Dissertation Figure 9, showing the proportional
# composition of LOH categories within each tumour-region sample

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
# Calculate total LOH fragments for each sample
# -----------------------------------------------------------------------------

sample_totals <- sample_category_summary %>%
  dplyr::group_by(sample) %>%
  dplyr::summarise(
    sample_total_LOH = sum(
      count,
      na.rm = TRUE
    ),
    .groups = "drop"
  )

if (any(sample_totals$sample_total_LOH <= 0)) {
  stop(
    "One or more samples have no LOH fragments and cannot be ",
    "displayed in a proportional stacked bar chart."
  )
}

# Order samples from the smallest total LOH count on the left
# to the largest total LOH count on the right.
ordered_samples <- sample_totals %>%
  dplyr::arrange(
    sample_total_LOH,
    sample
  ) %>%
  dplyr::pull(sample)

# -----------------------------------------------------------------------------
# Calculate the proportion of each category within each sample
# -----------------------------------------------------------------------------

figure_9_plot_data <- sample_category_summary %>%
  dplyr::left_join(
    sample_totals,
    by = "sample"
  ) %>%
  dplyr::mutate(
    category_proportion =
      count / sample_total_LOH,
    
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
    )
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
# Create the 100% stacked bar chart
# -----------------------------------------------------------------------------

figure_9 <- ggplot(
  figure_9_plot_data,
  aes(
    x = sample,
    y = category_proportion,
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
      accuracy = 1
    ),
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
    )
  ) +
  labs(
    title = paste(
      "percentage proportion of each category in each sample",
      "(descending order)"
    ),
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
    
    # Individual sample labels would not be readable for 434 bars.
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
  "figures/Figure9_sample_category_proportions.png"

figure_pdf <-
  "figures/Figure9_sample_category_proportions.pdf"

ggsave(
  filename = figure_png,
  plot = figure_9,
  width = 11,
  height = 7,
  units = "in",
  dpi = 300,
  bg = "white"
)

ggsave(
  filename = figure_pdf,
  plot = figure_9,
  width = 11,
  height = 7,
  units = "in",
  bg = "white"
)

# -----------------------------------------------------------------------------
# Validation
# -----------------------------------------------------------------------------

sample_proportion_check <- figure_9_plot_data %>%
  dplyr::group_by(sample) %>%
  dplyr::summarise(
    total_proportion = sum(
      category_proportion,
      na.rm = TRUE
    ),
    .groups = "drop"
  )

category_leaders <- figure_9_plot_data %>%
  dplyr::group_by(category) %>%
  dplyr::slice_max(
    order_by = count,
    n = 1,
    with_ties = TRUE
  ) %>%
  dplyr::select(
    category,
    sample,
    count,
    sample_total_LOH,
    category_proportion
  ) %>%
  dplyr::arrange(category)

figure_9_validation <- tibble::tibble(
  number_of_samples =
    dplyr::n_distinct(
      figure_9_plot_data$sample
    ),
  
  number_of_categories =
    dplyr::n_distinct(
      figure_9_plot_data$category
    ),
  
  sample_category_combinations =
    nrow(figure_9_plot_data),
  
  total_LOH_fragments =
    sum(
      sample_category_summary$count
    ),
  
  minimum_sample_total =
    min(
      sample_totals$sample_total_LOH
    ),
  
  maximum_sample_total =
    max(
      sample_totals$sample_total_LOH
    ),
  
  all_sample_proportions_equal_one =
    all(
      abs(
        sample_proportion_check$total_proportion - 1
      ) < 1e-12
    ),
  
  ascending_sample_order_confirmed =
    all(
      diff(
        sample_totals %>%
          dplyr::arrange(
            sample_total_LOH,
            sample
          ) %>%
          dplyr::pull(sample_total_LOH)
      ) >= 0
    )
)

print(
  figure_9_validation,
  n = Inf,
  width = Inf
)

print(
  category_leaders,
  n = Inf,
  width = Inf
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
    "One or more Figure 9 output files were not created."
  )
}

# -----------------------------------------------------------------------------
# Figure legend:
#
# Figure 9. Proportional composition of the eight LOH categories within each
# tumour-region sample.
#
# Each bar represents one sample and is scaled to 100%. The coloured sections
# show the proportion of that sample's LOH fragments assigned to categories
# 1–8. Samples are ordered by total LOH fragment count, with the smallest
# sample total on the left and the largest sample total on the right.
# -----------------------------------------------------------------------------