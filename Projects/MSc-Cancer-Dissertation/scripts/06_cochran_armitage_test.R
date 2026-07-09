# 06_cochran_armitage_test.R
# MSc Cancer Dissertation reconstruction
# Purpose: test for a trend between LOH category and clonality

# Load packages
library(dplyr)
library(readr)

# Load LOH category and clonality summary
loh_segments <- readRDS("data/processed/loh_segments_categories.rds")

# Prepare counts for trend test
trend_data <- loh_segments %>%
  filter(!is.na(LOH_category)) %>%
  count(LOH_category, clonality) %>%
  tidyr::pivot_wider(
    names_from = clonality,
    values_from = n,
    values_fill = 0
  ) %>%
  arrange(LOH_category)

# View prepared table
print(trend_data)

# Cochran-Armitage style test for trend
# Here, "Clonal" counts are tested across ordered LOH categories
trend_test <- prop.trend.test(
  x = trend_data$Clonal,
  n = trend_data$Clonal + trend_data$Subclonal,
  score = 1:nrow(trend_data)
)

# View result
print(trend_test)

# Save outputs
saveRDS(
  trend_test,
  "results/cochran_armitage_trend_test.rds"
)

write_csv(
  trend_data,
  "results/cochran_armitage_input_table.csv"
)
