# -----------------------------------------------------------------------------
# Original MSc Cancer Dissertation analysis (2021)
# Refactored for reproducibility and GitHub portfolio (2026)
# -----------------------------------------------------------------------------

# 17_figure4_major_allele_distribution.R
# Purpose: reproduce Dissertation Figure 4, showing the distribution of
# major-allele copy numbers in LOH, complete, and non-LOH fragment datasets

library(dplyr)
library(ggplot2)
library(patchwork)

# -----------------------------------------------------------------------------
# Input files
# -----------------------------------------------------------------------------

all_fragments_file <-
  "data/processed/tracerx_ascat_seg_categories.rds"

loh_fragments_file <-
  "data/processed/loh_segments_fragment_sizes.rds"

non_loh_fragments_file <-
  "data/processed/non_loh_segments.rds"

required_files <- c(
  all_fragments_file,
  loh_fragments_file,
  non_loh_fragments_file
)

missing_files <- required_files[
  !file.exists(required_files)
]

if (length(missing_files) > 0) {
  stop(
    "The following required input file(s) were not found:\n",
    paste(missing_files, collapse = "\n")
  )
}

# -----------------------------------------------------------------------------
# Load validated data
# -----------------------------------------------------------------------------

all_fragments <- readRDS(
  all_fragments_file
)

loh_fragments <- readRDS(
  loh_fragments_file
)

non_loh_fragments <- readRDS(
  non_loh_fragments_file
)

# Confirm that nMajor exists in all three datasets
datasets <- list(
  all_fragments = all_fragments,
  loh_fragments = loh_fragments,
  non_loh_fragments = non_loh_fragments
)

missing_nmajor <- names(datasets)[
  !vapply(
    datasets,
    function(x) "nMajor" %in% names(x),
    logical(1)
  )
]

if (length(missing_nmajor) > 0) {
  stop(
    "The column 'nMajor' was not found in: ",
    paste(missing_nmajor, collapse = ", ")
  )
}

# Remove missing nMajor observations if any are present
all_fragments_plot_data <- all_fragments %>%
  dplyr::filter(!is.na(nMajor))

loh_fragments_plot_data <- loh_fragments %>%
  dplyr::filter(!is.na(nMajor))

non_loh_fragments_plot_data <- non_loh_fragments %>%
  dplyr::filter(!is.na(nMajor))

# -----------------------------------------------------------------------------
# Shared dissertation-style theme
# -----------------------------------------------------------------------------

dissertation_theme <- theme_gray(
  base_size = 11
) +
  theme(
    plot.title = element_text(
      size = 12,
      face = "plain",
      hjust = 0
    ),
    axis.title = element_text(
      size = 10
    ),
    axis.text = element_text(
      size = 9
    ),
    plot.tag = element_text(
      size = 22,
      face = "bold",
      hjust = 0,
      vjust = 1
    ),
    plot.tag.position = "topleft",
    plot.margin = margin(
      t = 12,
      r = 10,
      b = 8,
      l = 10
    )
  )

# -----------------------------------------------------------------------------
# Panel A: nMajor distribution in LOH fragments
# -----------------------------------------------------------------------------

plot_A <- ggplot(
  loh_fragments_plot_data,
  aes(x = nMajor)
) +
  geom_histogram(
    binwidth = 1,
    boundary = -0.5,
    fill = "blue",
    colour = "black",
    linewidth = 0.35
  ) +
  scale_x_continuous(
    breaks = c(0, 20, 40, 60, 80),
    expand = expansion(mult = c(0, 0.02))
  ) +
  coord_cartesian(
    xlim = c(-1, 85)
  ) +
  labs(
    tag = "A",
    title = "Copy number of Major alleles in all LOH samples",
    x = "nMajor",
    y = "count"
  ) +
  dissertation_theme

# -----------------------------------------------------------------------------
# Panel B: nMajor distribution in all ASCAT fragments
# -----------------------------------------------------------------------------

plot_B <- ggplot(
  all_fragments_plot_data,
  aes(x = nMajor)
) +
  geom_histogram(
    binwidth = 1,
    boundary = -0.5,
    fill = "blue",
    colour = "black",
    linewidth = 0.35
  ) +
  scale_x_continuous(
    breaks = c(0, 25, 50, 75, 100),
    expand = expansion(mult = c(0, 0.02))
  ) +
  coord_cartesian(
    xlim = c(-1, 100)
  ) +
  labs(
    tag = "B",
    title = "Copy number of Major alleles",
    x = "nMajor",
    y = "count"
  ) +
  dissertation_theme

# -----------------------------------------------------------------------------
# Panel C: nMajor distribution in non-LOH fragments
# -----------------------------------------------------------------------------

plot_C <- ggplot(
  non_loh_fragments_plot_data,
  aes(x = nMajor)
) +
  geom_histogram(
    binwidth = 1,
    boundary = -0.5,
    fill = "blue",
    colour = "black",
    linewidth = 0.35
  ) +
  scale_x_continuous(
    breaks = c(0, 25, 50, 75, 100),
    expand = expansion(mult = c(0, 0.02))
  ) +
  coord_cartesian(
    xlim = c(-1, 100)
  ) +
  labs(
    tag = "C",
    title = "Copy number of Major alleles in non-LOH fragments",
    x = "nMajor",
    y = "count"
  ) +
  dissertation_theme

# -----------------------------------------------------------------------------
# Assemble the dissertation-style layout
#
#               A
#
#          B         C
# -----------------------------------------------------------------------------

figure_design <- "
#AA#
BBCC
"

figure_4 <-
  plot_A +
  plot_B +
  plot_C +
  plot_layout(
    design = figure_design,
    heights = c(1, 1)
  )

# -----------------------------------------------------------------------------
# Save the combined figure
# -----------------------------------------------------------------------------

dir.create(
  "figures",
  recursive = TRUE,
  showWarnings = FALSE
)

figure_png <-
  "figures/Figure4_major_allele_distribution.png"

figure_pdf <-
  "figures/Figure4_major_allele_distribution.pdf"

ggsave(
  filename = figure_png,
  plot = figure_4,
  width = 13,
  height = 9,
  units = "in",
  dpi = 300,
  bg = "white"
)

ggsave(
  filename = figure_pdf,
  plot = figure_4,
  width = 13,
  height = 9,
  units = "in",
  bg = "white"
)

# -----------------------------------------------------------------------------
# Validation
# -----------------------------------------------------------------------------

figure_4_validation <- tibble::tibble(
  panel = c(
    "A: LOH fragments",
    "B: all fragments",
    "C: non-LOH fragments"
  ),
  fragment_count = c(
    nrow(loh_fragments_plot_data),
    nrow(all_fragments_plot_data),
    nrow(non_loh_fragments_plot_data)
  ),
  minimum_nMajor = c(
    min(loh_fragments_plot_data$nMajor),
    min(all_fragments_plot_data$nMajor),
    min(non_loh_fragments_plot_data$nMajor)
  ),
  maximum_nMajor = c(
    max(loh_fragments_plot_data$nMajor),
    max(all_fragments_plot_data$nMajor),
    max(non_loh_fragments_plot_data$nMajor)
  )
)

print(
  figure_4_validation,
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
    "One or more Figure 4 output files were not created."
  )
}

# -----------------------------------------------------------------------------
# Figure legend:
#
# Figure 4. Distribution of major-allele copy number.
#
# A: Copy number of major alleles in LOH fragments.
# B: Copy number of major alleles in all fragments, regardless of LOH status.
# C: Copy number of major alleles in non-LOH fragments.
#
# In all panels, count represents the number of genomic fragments with the
# corresponding number of major-allele copies, and nMajor represents the
# number of copies of the major allele.
# -----------------------------------------------------------------------------