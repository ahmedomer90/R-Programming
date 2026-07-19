# -----------------------------------------------------------------------------
# Original MSc Cancer Dissertation analysis (2021)
# Refactored for reproducibility and GitHub portfolio (2026)
# -----------------------------------------------------------------------------

# 36_figure32_category4_fragment_position_heatmap.R
#
# Purpose:
# Reproduce Dissertation Figure 32 as a chromosome-position heatmap showing
# the number of category 4 LOH fragments overlapping genomic intervals.
#
# Each chromosome is divided into 5-Mb genomic intervals. The script counts
# the number of category 4 fragments overlapping each interval across all
# tumour-region samples.

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

if (nrow(loh_segments) == 0) {
  stop(
    "The LOH segment dataset contains no rows."
  )
}

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
    paste(
      "One or more fragments have an end position",
      "below their start position."
    )
  )
}

if (
  any(
    !as.integer(loh_segments$chr) %in% 1:24
  )
) {
  stop(
    paste(
      "Unexpected chromosome codes were found.",
      "Expected chromosome codes 1–24."
    )
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

if (!"4" %in% available_categories) {
  stop(
    "The expected LOH category '4' was not found. ",
    "Available categories are: ",
    paste(
      available_categories,
      collapse = ", "
    )
  )
}

# -----------------------------------------------------------------------------
# Filter category 4 LOH fragments
# -----------------------------------------------------------------------------

category4_segments <- loh_segments %>%
  dplyr::filter(
    as.character(category) == "4"
  ) %>%
  dplyr::transmute(
    fragment_id = dplyr::row_number(),
    chr = as.integer(chr),
    startpos = as.numeric(startpos),
    endpos = as.numeric(endpos)
  )

if (nrow(category4_segments) == 0) {
  stop(
    "No category 4 LOH fragments were found."
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
    249250621,  # chromosome 1
    243199373,  # chromosome 2
    198022430,  # chromosome 3
    191154276,  # chromosome 4
    180915260,  # chromosome 5
    171115067,  # chromosome 6
    159138663,  # chromosome 7
    146364022,  # chromosome 8
    141213431,  # chromosome 9
    135534747,  # chromosome 10
    135006516,  # chromosome 11
    133851895,  # chromosome 12
    115169878,  # chromosome 13
    107349540,  # chromosome 14
    102531392,  # chromosome 15
    90354753,   # chromosome 16
    81195210,   # chromosome 17
    78077248,   # chromosome 18
    59128983,   # chromosome 19
    63025520,   # chromosome 20
    48129895,   # chromosome 21
    51304566,   # chromosome 22
    155270560,  # chromosome X
    59373566    # chromosome Y
  )
)

# -----------------------------------------------------------------------------
# Validate fragment coordinates against chromosome lengths
# -----------------------------------------------------------------------------

coordinate_validation <- category4_segments %>%
  dplyr::left_join(
    chromosome_lengths,
    by = "chr"
  )

if (any(is.na(coordinate_validation$chromosome_length_bp))) {
  stop(
    paste(
      "One or more category 4 fragments could not",
      "be matched to a chromosome."
    )
  )
}

fragments_beyond_chromosome <- coordinate_validation %>%
  dplyr::filter(
    endpos > chromosome_length_bp
  )

if (nrow(fragments_beyond_chromosome) > 0) {
  warning(
    nrow(fragments_beyond_chromosome),
    " category 4 fragment(s) extend beyond the listed GRCh37 chromosome ",
    "length. Their end coordinates will be limited to the chromosome length."
  )
}

category4_segments <- coordinate_validation %>%
  dplyr::mutate(
    startpos = pmax(
      startpos,
      0
    ),
    
    endpos = pmin(
      endpos,
      chromosome_length_bp
    )
  ) %>%
  dplyr::select(
    fragment_id,
    chr,
    startpos,
    endpos
  )

# -----------------------------------------------------------------------------
# Define genomic bin width
# -----------------------------------------------------------------------------

bin_width_bp <- 5000000

# -----------------------------------------------------------------------------
# Create complete chromosome-bin framework
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
  ) %>%
  dplyr::arrange(
    chr,
    bin_index
  )

# -----------------------------------------------------------------------------
# Assign every category 4 fragment to all genomic bins it overlaps
# -----------------------------------------------------------------------------

category4_fragment_bins <- category4_segments %>%
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
# Count category 4 fragments overlapping each chromosome-position bin
# -----------------------------------------------------------------------------

category4_counts_by_bin <- category4_fragment_bins %>%
  dplyr::count(
    chr,
    bin_index,
    name = "category4_fragment_count"
  )

# -----------------------------------------------------------------------------
# Combine overlap counts with all chromosome-position bins
# -----------------------------------------------------------------------------

figure_32_all_bins <- chromosome_bins %>%
  dplyr::left_join(
    category4_counts_by_bin,
    by = c(
      "chr",
      "bin_index"
    )
  ) %>%
  dplyr::mutate(
    category4_fragment_count =
      tidyr::replace_na(
        category4_fragment_count,
        0L
      )
  ) %>%
  dplyr::arrange(
    chr,
    bin_index
  )

# Only bins containing at least one category 4 fragment are drawn.
# Bins with zero fragments therefore appear blank.

figure_32_plot_data <- figure_32_all_bins %>%
  dplyr::filter(
    category4_fragment_count > 0
  )

if (nrow(figure_32_plot_data) == 0) {
  stop(
    "No chromosome-position bins contained category 4 fragments."
  )
}

# -----------------------------------------------------------------------------
# Determine colour-scale limits and legend breaks
# -----------------------------------------------------------------------------

maximum_category4_count <- max(
  figure_32_plot_data$category4_fragment_count
)

# The original dissertation Figure 32 displays legend labels around
# 20, 40 and 60 fragments.

legend_breaks <- c(
  20,
  40,
  60
)

legend_breaks <- legend_breaks[
  legend_breaks <= maximum_category4_count
]

if (length(legend_breaks) == 0) {
  legend_breaks <- pretty(
    c(
      0,
      maximum_category4_count
    ),
    n = 3
  )
  
  legend_breaks <- legend_breaks[
    legend_breaks > 0 &
      legend_breaks <= maximum_category4_count
  ]
}

# -----------------------------------------------------------------------------
# Create Figure 32
# -----------------------------------------------------------------------------

figure_32 <- ggplot(
  figure_32_plot_data,
  aes(
    x = chr,
    y = bin_midpoint_bp,
    fill = category4_fragment_count
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
    name = "number of category 4 fragments",
    breaks = legend_breaks,
    limits = c(
      0,
      maximum_category4_count
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
      "number of category 4 fragments in the genome",
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

# Display Figure 32 in the RStudio Plots pane.

print(
  figure_32
)

# -----------------------------------------------------------------------------
# Save Figure 32
# -----------------------------------------------------------------------------

dir.create(
  "figures",
  recursive = TRUE,
  showWarnings = FALSE
)

figure_png <-
  "figures/Figure32_category4_fragment_position_heatmap.png"

figure_pdf <-
  "figures/Figure32_category4_fragment_position_heatmap.pdf"

ggsave(
  filename = figure_png,
  plot = figure_32,
  width = 10,
  height = 7,
  units = "in",
  dpi = 300,
  bg = "white"
)

ggsave(
  filename = figure_pdf,
  plot = figure_32,
  width = 10,
  height = 7,
  units = "in",
  bg = "white"
)

# -----------------------------------------------------------------------------
# Count original category 4 fragments on each chromosome
# -----------------------------------------------------------------------------

category4_fragments_by_chromosome <- category4_segments %>%
  dplyr::count(
    chr,
    name = "number_of_category4_fragments"
  ) %>%
  tibble::as_tibble()

# -----------------------------------------------------------------------------
# Chromosome-level validation table
# -----------------------------------------------------------------------------

chromosome_category4_summary <- figure_32_all_bins %>%
  dplyr::group_by(
    chr
  ) %>%
  dplyr::summarise(
    chromosome_length_bp =
      dplyr::first(
        chromosome_length_bp
      ),
    
    number_of_bins =
      dplyr::n(),
    
    bins_with_category4 =
      sum(
        category4_fragment_count > 0
      ),
    
    maximum_category4_fragments_in_bin =
      max(
        category4_fragment_count
      ),
    
    sum_of_category4_bin_overlap_counts =
      sum(
        category4_fragment_count
      ),
    
    .groups = "drop"
  ) %>%
  dplyr::left_join(
    category4_fragments_by_chromosome,
    by = "chr"
  ) %>%
  dplyr::mutate(
    number_of_category4_fragments =
      tidyr::replace_na(
        number_of_category4_fragments,
        0L
      )
  ) %>%
  dplyr::select(
    chr,
    chromosome_length_bp,
    number_of_bins,
    number_of_category4_fragments,
    bins_with_category4,
    maximum_category4_fragments_in_bin,
    sum_of_category4_bin_overlap_counts
  ) %>%
  tibble::as_tibble()

# -----------------------------------------------------------------------------
# Overall Figure 32 validation
# -----------------------------------------------------------------------------

figure_32_validation <- tibble::tibble(
  number_of_input_LOH_fragments =
    nrow(
      loh_segments
    ),
  
  number_of_category4_fragments =
    nrow(
      category4_segments
    ),
  
  number_of_chromosomes_in_category4_source =
    dplyr::n_distinct(
      category4_segments$chr
    ),
  
  genomic_bin_width_bp =
    bin_width_bp,
  
  total_number_of_bins =
    nrow(
      figure_32_all_bins
    ),
  
  bins_with_at_least_one_category4_fragment =
    nrow(
      figure_32_plot_data
    ),
  
  maximum_category4_fragments_in_one_bin =
    max(
      figure_32_plot_data$
        category4_fragment_count
    ),
  
  minimum_positive_category4_fragments_in_one_bin =
    min(
      figure_32_plot_data$
        category4_fragment_count
    ),
  
  chromosomes_with_plotted_category4 =
    dplyr::n_distinct(
      figure_32_plot_data$chr
    ),
  
  negative_fragment_counts =
    sum(
      figure_32_all_bins$
        category4_fragment_count < 0
    )
)

print(
  figure_32_validation,
  n = Inf,
  width = Inf
)

print(
  chromosome_category4_summary,
  n = Inf,
  width = Inf
)

# -----------------------------------------------------------------------------
# Validate output files
# -----------------------------------------------------------------------------

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
  n = Inf,
  width = Inf
)

if (!all(figure_file_validation$file_exists)) {
  warning(
    "One or more Figure 32 output files were not created."
  )
}

if (
  any(
    figure_file_validation$file_size_bytes <= 0,
    na.rm = TRUE
  )
) {
  warning(
    "One or more Figure 32 output files have a zero file size."
  )
}

# -----------------------------------------------------------------------------
# Optional validation: counts for all LOH categories
# -----------------------------------------------------------------------------

category_count_validation <- loh_segments %>%
  dplyr::count(
    category,
    name = "number_of_fragments"
  ) %>%
  dplyr::arrange(
    as.numeric(
      as.character(category)
    )
  ) %>%
  tibble::as_tibble()

print(
  category_count_validation,
  n = Inf,
  width = Inf
)

# -----------------------------------------------------------------------------
# Figure legend
# -----------------------------------------------------------------------------
#
# Figure 32. Chromosome-position heatmap conveying the total number of
# category 4 LOH fragments overlapping each genomic region of each chromosome
# across all patients and tumour-region samples.
#
# Each chromosome was divided into 5-Mb genomic intervals. Lighter blue
# intervals indicate fewer overlapping category 4 fragments, while darker
# blue intervals indicate genomic positions affected by a larger number of
# category 4 fragments. Blank regions indicate genomic bins containing no
# category 4 LOH fragments.
# -----------------------------------------------------------------------------