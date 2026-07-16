# -----------------------------------------------------------------------------
# Original MSc Cancer Dissertation analysis (2021)
# Refactored for reproducibility and GitHub portfolio (2026)
# -----------------------------------------------------------------------------

# 06_cochran_armitage_test.R
# Purpose: test for a trend in clonality across ordered LOH categories

library(dplyr)
library(tidyr)
library(readr)
library(DescTools)

# Load the validated category and clonality summary from Script 5
category_clonality_summary <- readRDS(
  "data/processed/category_clonality_summary.rds"
)

# Confirm required columns exist
required_columns <- c("category", "clonality", "count")

if (!all(required_columns %in% names(category_clonality_summary))) {
  stop(
    "Required columns 'category', 'clonality', and/or 'count' were not found."
  )
}

# Record rows with clonality classified as "none"
none_cases <- category_clonality_summary %>%
  filter(clonality == "none")

print(none_cases)

# Prepare the ordered category-by-clonality table
# The Cochran-Armitage test requires exactly two outcome columns.
trend_data <- category_clonality_summary %>%
  filter(clonality %in% c("clonal", "subclonal")) %>%
  select(category, clonality, count) %>%
  pivot_wider(
    names_from = clonality,
    values_from = count,
    values_fill = 0
  ) %>%
  mutate(category = as.integer(category)) %>%
  arrange(category)

print(trend_data)

# Convert clonal and subclonal counts into a numeric matrix
trend_matrix <- as.matrix(
  trend_data[, c("clonal", "subclonal")]
)

rownames(trend_matrix) <- trend_data$category

print(trend_matrix)

# Perform Cochran-Armitage test for trend
cochran_armitage_test <- CochranArmitageTest(
  trend_matrix,
  alternative = "two.sided"
)

# Display the complete test and p-value
print(cochran_armitage_test)
print(cochran_armitage_test$p.value)

# Create results directory if needed
dir.create("results", recursive = TRUE, showWarnings = FALSE)

# Save test inputs and result
write_csv(
  trend_data,
  "results/cochran_armitage_input_table.csv"
)

saveRDS(
  cochran_armitage_test,
  "results/cochran_armitage_test.rds"
)

# Save a readable text version of the test
capture.output(
  cochran_armitage_test,
  file = "results/cochran_armitage_test.txt"
)

# -----------------------------------------------------------------------------
# Validation note:
# Ten rare fragments classified as "none" by Script 3 are excluded from this
# test because the Cochran-Armitage procedure requires an ordered r x 2 table.
# The test therefore compares clonal versus subclonal fragments across the
# eight ordered LOH categories.
# -----------------------------------------------------------------------------