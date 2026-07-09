# 15_export_session_info.R
# MSc Cancer Dissertation reconstruction
# Purpose: export R session information for reproducibility

# Create results directory if needed
dir.create("results", showWarnings = FALSE)

# Capture session information
session_info <- capture.output(sessionInfo())

# Save session information
writeLines(
  session_info,
  "results/sessionInfo.txt"
)
