# -----------------------------------------------------------------------------
# Original MSc Cancer Dissertation analysis (2021)
# Refactored for reproducibility and GitHub portfolio (2026)
# -----------------------------------------------------------------------------

# 32_figure28_LOH_fragment_position_heatmap.R
#
# Purpose:
# Reproduce Dissertation Figure 28 as a chromosome-position heatmap showing
# how many LOH fragments overlap each genomic interval, regardless of LOH
# category or clonality.

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
  "endpos"
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
# Validate source coordinates
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
    "One or more fragments have an end position below their start position."
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
# Human chromosome lengths
# -----------------------------------------------------------------------------
#
# TRACERx ASCAT data used human genome build GRCh37/hg19.
#
# Chromosome coding used in this reconstructed dataset:
#   1–22 = autosomes
#   23   = chromosome X
#   24   = chromosome Y
#
# Chromosome lengths allow the plot to preserve the recognisable descending
# chromosome shapes visible in the original dissertation figure.

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
#
# A 5-Mb bin width gives a block-like heatmap close to the dissertation
# figure while preserving positional detail.

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
  )

# -----------------------------------------------------------------------------
# Assign every LOH fragment to all genomic bins that it overlaps
# -----------------------------------------------------------------------------

loh_fragment_bins <- loh_segments %>%
  dplyr::transmute(
    fragment_id = dplyr::row_number(),
    
    chr = as.integer(chr),
    
    startpos = as.numeric(startpos),
    
    endpos = as.numeric(endpos),
    
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
# Count the number of LOH fragments overlapping every chromosome-position bin
# -----------------------------------------------------------------------------

fragment_counts_by_bin <- loh_fragment_bins %>%
  dplyr::count(
    chr,
    bin_index,
    name = "fragment_count"
  )

# -----------------------------------------------------------------------------
# Combine counts with the complete chromosome-bin framework
# -----------------------------------------------------------------------------

figure_28_plot_data <- chromosome_bins %>%
  dplyr::left_join(
    fragment_counts_by_bin,
    by = c(
      "chr",
      "bin_index"
    )
  ) %>%
  dplyr::mutate(
    fragment_count =
      tidyr::replace_na(
        fragment_count,
        0L
      )
  ) %>%
  dplyr::arrange(
    chr,
    bin_index
  )

# -----------------------------------------------------------------------------
# Validate prepared plotting data
# -----------------------------------------------------------------------------

if (
  any(
    figure_28_plot_data$
    fragment_count < 0
  )
) {
  stop(
    "Negative fragment counts were detected in the plotting data."
  )
}

if (
  any(
    figure_28_plot_data$
    bin_midpoint_bp >
    figure_28_plot_data$
    chromosome_length_bp
  )
) {
  stop(
    "One or more plotted bins extend beyond the chromosome length."
  )
}

# -----------------------------------------------------------------------------
# Create chromosome-position heatmap
# -----------------------------------------------------------------------------

figure_28 <- ggplot(
  figure_28_plot_data,
  aes(
    x = chr,
    y = bin_midpoint_bp,
    fill = fragment_count
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
    name = "number of fragments",
    breaks = c(
      50,
      100,
      150,
      200,
      250
    ),
    limits = c(
      0,
      max(
        figure_28_plot_data$
          fragment_count
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
      "number of fragments on each position",
      "in the chromosome with LOH"
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
  "figures/Figure28_LOH_fragment_position_heatmap.png"

figure_pdf <-
  "figures/Figure28_LOH_fragment_position_heatmap.pdf"

ggsave(
  filename = figure_png,
  plot = figure_28,
  width = 10,
  height = 7,
  units = "in",
  dpi = 300,
  bg = "white"
)

ggsave(
  filename = figure_pdf,
  plot = figure_28,
  width = 10,
  height = 7,
  units = "in",
  bg = "white"
)

# -----------------------------------------------------------------------------
# Validation
# -----------------------------------------------------------------------------

chromosome_heatmap_summary <- figure_28_plot_data %>%
  dplyr::group_by(chr) %>%
  dplyr::summarise(
    chromosome_length_bp =
      dplyr::first(
        chromosome_length_bp
      ),
    
    number_of_bins =
      dplyr::n(),
    
    bins_with_LOH =
      sum(
        fragment_count > 0
      ),
    
    maximum_bin_fragment_count =
      max(
        fragment_count
      ),
    
    sum_of_bin_overlap_counts =
      sum(
        fragment_count
      ),
    
    .groups = "drop"
  )

figure_28_validation <- tibble::tibble(
  number_of_input_LOH_fragments =
    nrow(
      loh_segments
    ),
  
  number_of_chromosomes =
    dplyr::n_distinct(
      figure_28_plot_data$chr
    ),
  
  genomic_bin_width_bp =
    bin_width_bp,
  
  total_number_of_bins =
    nrow(
      figure_28_plot_data
    ),
  
  bins_with_at_least_one_LOH_fragment =
    sum(
      figure_28_plot_data$
        fragment_count > 0
    ),
  
  maximum_fragments_in_one_bin =
    max(
      figure_28_plot_data$
        fragment_count
    ),
  
  minimum_fragments_in_one_bin =
    min(
      figure_28_plot_data$
        fragment_count
    ),
  
  chromosomes_with_plotted_LOH =
    dplyr::n_distinct(
      figure_28_plot_data$chr[
        figure_28_plot_data$
          fragment_count > 0
      ]
    ),
  
  negative_fragment_counts =
    sum(
      figure_28_plot_data$
        fragment_count < 0
    )
)

print(
  figure_28_validation,
  n = Inf,
  width = Inf
)

print(
  chromosome_heatmap_summary,
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

if (
  !all(
    figure_file_validation$
    file_exists
  )
) {
  warning(
    "One or more Figure 28 output files were not created."
  )
}

# -----------------------------------------------------------------------------
# Figure legend:
#
# Figure 28. Chromosome-position heatmap showing the number of LOH fragments
# overlapping genomic intervals on each chromosome, regardless of LOH category
# or clonality.
#
# Each chromosome is divided into 5-Mb intervals. Lighter blue intervals
# indicate fewer overlapping LOH fragments, while darker blue intervals
# indicate genomic positions affected by larger numbers of LOH fragments.
# -----------------------------------------------------------------------------