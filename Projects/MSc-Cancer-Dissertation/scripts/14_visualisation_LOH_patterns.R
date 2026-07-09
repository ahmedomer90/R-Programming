# 14_visualisation_LOH_patterns.R
# MSc Cancer Dissertation reconstruction
# Purpose: generate plots showing LOH category, clonality, chromosome, sample, and fragment-size patterns

# Load packages
library(dplyr)
library(readr)
library(ggplot2)

# Load processed data
loh_segments <- readRDS("data/processed/loh_segments_fragment_sizes.rds")
category_clonality_summary <- readRDS("data/processed/category_clonality_summary.rds")
chromosome_clonality_summary <- readRDS("data/processed/chromosome_clonality_summary.rds")
sample_clonality_summary <- readRDS("data/processed/sample_clonality_summary.rds")

# Create figures directory if it does not exist
dir.create("figures", showWarnings = FALSE)

# Plot 1: LOH category distribution
plot_category_distribution <- loh_segments %>%
  ggplot(aes(x = LOH_category)) +
  geom_bar() +
  xlab("LOH category") +
  ylab("Number of LOH fragments") +
  ggtitle("Distribution of LOH categories")

ggsave(
  filename = "figures/LOH_category_distribution.png",
  plot = plot_category_distribution,
  width = 8,
  height = 5,
  dpi = 300
)

# Plot 2: LOH category by clonality
plot_category_clonality <- category_clonality_summary %>%
  ggplot(aes(x = LOH_category, y = proportion, fill = clonality)) +
  geom_col() +
  xlab("LOH category") +
  ylab("Proportion") +
  ggtitle("Proportion of clonal and subclonal LOH by category")

ggsave(
  filename = "figures/LOH_category_clonality_proportion.png",
  plot = plot_category_clonality,
  width = 8,
  height = 5,
  dpi = 300
)

# Plot 3: Chromosome by clonality
plot_chromosome_clonality <- chromosome_clonality_summary %>%
  ggplot(aes(x = factor(chr), y = proportion, fill = clonality)) +
  geom_col() +
  xlab("Chromosome") +
  ylab("Proportion") +
  ggtitle("Proportion of clonal and subclonal LOH by chromosome")

ggsave(
  filename = "figures/chromosome_clonality_proportion.png",
  plot = plot_chromosome_clonality,
  width = 10,
  height = 5,
  dpi = 300
)

# Plot 4: Fragment size by LOH category
plot_fragment_size <- loh_segments %>%
  ggplot(aes(x = LOH_category, y = fragment_size_mb)) +
  geom_boxplot() +
  xlab("LOH category") +
  ylab("Fragment size (Mb)") +
  ggtitle("Fragment size distribution by LOH category")

ggsave(
  filename = "figures/fragment_size_by_LOH_category.png",
  plot = plot_fragment_size,
  width = 8,
  height = 5,
  dpi = 300
)
