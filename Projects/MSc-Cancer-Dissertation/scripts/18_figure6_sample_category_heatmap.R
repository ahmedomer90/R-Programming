# -----------------------------------------------------------------------------
# Original MSc Cancer Dissertation analysis (2021)
# Refactored for reproducibility and GitHub portfolio (2026)
# -----------------------------------------------------------------------------

# 18_figure6_sample_category_heatmap.R
# Purpose: reproduce Dissertation Figure 6, showing the number of LOH
# fragments in each category for every tumour-region sample

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

sample_category_plot_data <-
  sample_category_summary %>%
  dplyr::mutate(
    category = factor(
      category,
      levels = as.character(1:8),
      labels = paste0(
        "category",
        1:8
      )
    ),
    
    # Preserve the sample ordering present in the validated dataset
    sample = factor(
      sample,
      levels = rev(
        unique(sample)
      )
    )
  )

# -----------------------------------------------------------------------------
# Create dissertation-style heatmap
# -----------------------------------------------------------------------------

figure_6 <- ggplot(
  sample_category_plot_data,
  aes(
    x = category,
    y = sample,
    fill = count
  )
) +
  geom_tile() +
  scale_fill_gradient(
    low = "lightblue",
    high = "darkblue",
    name = "number of fragments"
  ) +
  labs(
    title = "number of each category of LOH in whole data",
    x = "category",
    y = "sample"
  ) +
  theme_minimal(
    base_size = 12
  ) +
  theme(
    panel.grid = element_blank(),
    
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
      size = 9,
      angle = 0,
      hjust = 0.5
    ),
    
    # The original figure did not display hundreds of sample labels
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    
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
  "figures/Figure6_sample_category_heatmap.png"

figure_pdf <-
  "figures/Figure6_sample_category_heatmap.pdf"

ggsave(
  filename = figure_png,
  plot = figure_6,
  width = 10,
  height = 7,
  units = "in",
  dpi = 300,
  bg = "white"
)

ggsave(
  filename = figure_pdf,
  plot = figure_6,
  width = 10,
  height = 7,
  units = "in",
  bg = "white"
)

# -----------------------------------------------------------------------------
# Validation
# -----------------------------------------------------------------------------

figure_6_validation <- tibble::tibble(
  number_of_samples =
    dplyr::n_distinct(
      sample_category_plot_data$sample
    ),
  
  number_of_categories =
    dplyr::n_distinct(
      sample_category_plot_data$category
    ),
  
  sample_category_combinations =
    nrow(sample_category_plot_data),
  
  total_LOH_fragments =
    sum(sample_category_plot_data$count),
  
  maximum_fragments_in_one_cell =
    max(sample_category_plot_data$count)
)

print(
  figure_6_validation,
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
    "One or more Figure 6 output files were not created."
  )
}

# -----------------------------------------------------------------------------
# Figure legend:
#
# Figure 6. Number of each category of LOH in all LOH samples in the complete
# tracerx.ascat.seg dataset.
#
# The sample variable represents individual tumour-region samples with LOH.
# The number of fragments represents the number of fragments in each sample
# belonging to a particular LOH category. The category variable represents
# the eight distinct LOH categories.
# -----------------------------------------------------------------------------
