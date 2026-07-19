# -----------------------------------------------------------------------------
# Original MSc Cancer Dissertation analysis (2021)
# Refactored for reproducibility and GitHub portfolio (2026)
# -----------------------------------------------------------------------------

# 33_figure29_category1_fragment_position_heatmap.R
#
# Purpose:
# Reproduce Dissertation Figure 29 as a chromosome-position heatmap showing
# the number of category 1 LOH fragments overlapping genomic intervals.
#
# Category 1 represents homozygous deletion in the reconstructed LOH
# classification workflow.

library(dplyr)
library(tidyr)
library(ggplot2)
library(tibble)

# -----------------------------------------------------------------------------
# Input file
# -----------------------------------------------------------------------------

loh_segments_file <-
  "data/processed/loh_segments_fragment_sizes.rds"

if (!file.exists(loh_segments_file)) {
  stop(
    "Required input file was not found: ",
    loh_segments_file
  )
}

# -----------------------------------------------------------------------------
# Load processed LOH fragments
# -----------------------------------------------------------------------------

loh_segments <- readRDS(
  loh_segments_file
)

required_columns <- c(
  "chr",
  "startpos",
  "endpos",
  "category"
)

missing_columns <- setdiff(
  required_columns,
  names(loh_segments)
)

if (length(missing_columns) > 0) {
  stop(
    "The LOH segment dataset is missing the following column(s): ",
    paste(
      missing_columns,
      collapse = ", "
    )
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

if (
  any(
    is.na(loh_segments$startpos) |
    is.na(loh_segments$endpos)
  )
) {
  stop(
    "Missing genomic start or end coordinates were detected."
  )
}

if (
  any(
    loh_segments$startpos < 0 |
    loh_segments$endpos < 0
  )
) {
  stop(
    "Negative genomic coordinates were detected."
  )
}

if (
  any(
    loh_segments$endpos <
    loh_segments$startpos
  )
) {
  stop(
    "One or more fragments have an end position below its start position."
  )
}

if (
  any(
    !loh_segments$chr %in% 1:24
  )
) {
  stop(
    "Unexpected chromosome codes were found. Expected chromosome codes 1–24."
  )
}

# -----------------------------------------------------------------------------
# Inspect available LOH categories
# -----------------------------------------------------------------------------

available_categories <- sort(
  unique(
    as.character(
      loh_segments$category
    )
  )
)

if (!"1" %in% available_categories) {
  stop(
    "The expected LOH category '1' was not found. ",
    "Available categories are: ",
    paste(
      available_categories,
      collapse = ", "
    )
  )
}

# -----------------------------------------------------------------------------
# Filter category 1 LOH fragments
# -----------------------------------------------------------------------------

category1_segments <- loh_segments %>%
  dplyr::filter(
    category == "1"
  ) %>%
  dplyr::transmute(
    fragment_id = dplyr::row_number(),
    chr = as.integer(chr),
    startpos = as.numeric(startpos),
    endpos = as.numeric(endpos)
  )

if (nrow(category1_segments) == 0) {
  stop(
    "No category 1 LOH fragments were found."
  )
}

# -----------------------------------------------------------------------------
# Human chromosome lengths: GRCh37 / hg19
# -----------------------------------------------------------------------------
#
# Chromosome coding:
#   1–22 = autosomes
#   23   = chromosome X
#   24   = chromosome Y

chromosome_lengths <- tibble::tibble(
  chr = 1:24,
  
  chromosome_length_bp = c(
    249250621,  # 1
    243199373,  # 2
    198022430,  # 3
    191154276,  # 4
    180915260,  # 5
    171115067,  # 6
    159138663,  # 7
    146364022,  # 8
    141213431,  # 9
    135534747,  # 10
    135006516,  # 11
    133851895,  # 12
    115169878,  # 13
    107349540,  # 14
    102531392,  # 15
    90354753,   # 16
    81195210,   # 17
    78077248,   # 18
    59128983,   # 19
    63025520,   # 20
    48129895,   # 21
    51304566,   # 22
    155270560,  # X
    59373566    # Y
  )
)

# -----------------------------------------------------------------------------
# Define genomic bin width
# -----------------------------------------------------------------------------

bin_width_bp <- 5000000

# -----------------------------------------------------------------------------
# Create a complete chromosome-bin framework
# -----------------------------------------------------------------------------

chromosome_bins <- chromosome_lengths %>%
  dplyr::mutate(
    number_of_bins = ceiling(
      chromosome_length_bp /
        bin_width_bp
    )
  ) %>%
  tidyr::uncount(
    weights = number_of_bins,
    .id = "bin_index"
  ) %>%
  dplyr::mutate(
    bin_start_bp =
      (bin_index - 1) *
      bin_width_bp,
    
    bin_end_bp = pmin(
      bin_index *
        bin_width_bp,
      chromosome_length_bp
    ),
    
    bin_midpoint_bp =
      (
        bin_start_bp +
          bin_end_bp
      ) / 2
  )

# -----------------------------------------------------------------------------
# Assign each category 1 fragment to every genomic bin it overlaps
# -----------------------------------------------------------------------------

category1_fragment_bins <- category1_segments %>%
  dplyr::mutate(
    start_bin = floor(
      pmax(
        startpos - 1,
        0
      ) /
        bin_width_bp
    ) + 1,
    
    end_bin = floor(
      pmax(
        endpos - 1,
        0
      ) /
        bin_width_bp
    ) + 1
  ) %>%
  dplyr::mutate(
    bin_index = Map(
      function(first_bin, last_bin) {
        seq.int(
          from = first_bin,
          to = last_bin
        )
      },
      start_bin,
      end_bin
    )
  ) %>%
  tidyr::unnest_longer(
    bin_index,
    values_to = "bin_index"
  ) %>%
  dplyr::mutate(
    bin_index = as.integer(
      bin_index
    )
  )

# -----------------------------------------------------------------------------
# Count category 1 fragments overlapping each chromosome-position bin
# -----------------------------------------------------------------------------

category1_counts_by_bin <- category1_fragment_bins %>%
  dplyr::count(
    chr,
    bin_index,
    name = "category1_fragment_count"
  )

# -----------------------------------------------------------------------------
# Combine counts with chromosome positions
# -----------------------------------------------------------------------------

figure_29_all_bins <- chromosome_bins %>%
  dplyr::left_join(
    category1_counts_by_bin,
    by = c(
      "chr",
      "bin_index"
    )
  ) %>%
  dplyr::mutate(
    category1_fragment_count =
      tidyr::replace_na(
        category1_fragment_count,
        0L
      )
  ) %>%
  dplyr::arrange(
    chr,
    bin_index
  )

# Only non-zero bins are drawn to reproduce the sparse appearance
# of the dissertation figure.
figure_29_plot_data <- figure_29_all_bins %>%
  dplyr::filter(
    category1_fragment_count > 0
  )

if (nrow(figure_29_plot_data) == 0) {
  stop(
    "No chromosome-position bins contained category 1 fragments."
  )
}

# -----------------------------------------------------------------------------
# Create the category 1 chromosome-position heatmap
# -----------------------------------------------------------------------------

figure_29 <- ggplot(
  figure_29_plot_data,
  aes(
    x = chr,
    y = bin_midpoint_bp,
    fill = category1_fragment_count
  )
) +
  geom_tile(
    width = 0.96,
    height = bin_width_bp,
    colour = NA
  ) +
  scale_fill_gradient(
    low = "#ADD8E6",
    high = "#191970",
    name = "number of category 1 fragments",
    breaks = c(
      10,
      20,
      30,
      40
    ),
    limits = c(
      0,
      max(
        figure_29_plot_data$
          category1_fragment_count
      )
    ),
    oob = scales::squish
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
    limits = c(
      0,
      25
    ),
    expand = expansion(
      mult = c(
        0,
        0
      )
    )
  ) +
  scale_y_continuous(
    breaks = c(
      0,
      50000000,
      100000000,
      150000000,
      200000000,
      250000000
    ),
    limits = c(
      0,
      260000000
    ),
    expand = expansion(
      mult = c(
        0,
        0
      )
    ),
    labels = scales::label_scientific(
      digits = 1
    )
  ) +
  labs(
    title = paste(
      "number of category 1 fragments in the genome",
      "of patients with LOH"
    ),
    x = "chromosome",
    y = "position on chromosome"
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
    
    panel.grid.major = element_line(
      colour = "white",
      linewidth = 0.6
    ),
    
    panel.grid.minor = element_line(
      colour = "white",
      linewidth = 0.3
    ),
    
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
  "figures/Figure29_category1_fragment_position_heatmap.png"

figure_pdf <-
  "figures/Figure29_category1_fragment_position_heatmap.pdf"

ggsave(
  filename = figure_png,
  plot = figure_29,
  width = 10,
  height = 7,
  units = "in",
  dpi = 300,
  bg = "white"
)

ggsave(
  filename = figure_pdf,
  plot = figure_29,
  width = 10,
  height = 7,
  units = "in",
  bg = "white"
)

# -----------------------------------------------------------------------------
# Validation
# -----------------------------------------------------------------------------

chromosome_category1_summary <- figure_29_all_bins %>%
  dplyr::group_by(chr) %>%
  dplyr::summarise(
    chromosome_length_bp =
      dplyr::first(
        chromosome_length_bp
      ),
    
    number_of_bins =
      dplyr::n(),
    
    bins_with_category1 =
      sum(
        category1_fragment_count > 0
      ),
    
    maximum_category1_fragments_in_bin =
      max(
        category1_fragment_count
      ),
    
    sum_of_category1_bin_overlap_counts =
      sum(
        category1_fragment_count
      ),
    
    .groups = "drop"
  )

figure_29_validation <- tibble::tibble(
  number_of_input_LOH_fragments =
    nrow(
      loh_segments
    ),
  
  number_of_category1_fragments =
    nrow(
      category1_segments
    ),
  
  number_of_chromosomes_in_source =
    dplyr::n_distinct(
      category1_segments$chr
    ),
  
  genomic_bin_width_bp =
    bin_width_bp,
  
  total_number_of_bins =
    nrow(
      figure_29_all_bins
    ),
  
  bins_with_at_least_one_category1_fragment =
    nrow(
      figure_29_plot_data
    ),
  
  maximum_category1_fragments_in_one_bin =
    max(
      figure_29_plot_data$
        category1_fragment_count
    ),
  
  minimum_positive_category1_fragments_in_one_bin =
    min(
      figure_29_plot_data$
        category1_fragment_count
    ),
  
  chromosomes_with_plotted_category1 =
    dplyr::n_distinct(
      figure_29_plot_data$chr
    ),
  
  negative_fragment_counts =
    sum(
      figure_29_all_bins$
        category1_fragment_count < 0
    )
)

print(
  figure_29_validation,
  n = Inf,
  width = Inf
)

print(
  chromosome_category1_summary,
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
    "One or more Figure 29 output files were not created."
  )
}

# -----------------------------------------------------------------------------
# Figure legend:
#
# Figure 29. Chromosome-position heatmap showing the number of category 1 LOH
# fragments overlapping genomic intervals across all tumour-region samples.
#
# Category 1 represents homozygous deletion. Each chromosome was divided into
# 5-Mb intervals. Lighter blue regions indicate fewer overlapping category 1
# fragments, while darker blue regions indicate positions affected by more
# category 1 fragments. Empty regions indicate bins with no category 1 LOH.
# -----------------------------------------------------------------------------