# -----------------------------------------------------------------------------
# Original MSc Cancer Dissertation analysis (2021)
# Refactored for reproducibility and GitHub portfolio (2026)
# -----------------------------------------------------------------------------

# 16_export_session_info.R
# Purpose: export R session information for reproducibility

# Create results directory if needed
dir.create(
  "results",
  recursive = TRUE,
  showWarnings = FALSE
)

# Capture session information
session_info <- capture.output(
  sessionInfo()
)

# Export session information
writeLines(
  session_info,
  "results/sessionInfo.txt"
)

# Validation
cat(
  "Session information exported successfully.\n"
)

cat(
  "Output file exists:",
  file.exists("results/sessionInfo.txt"),
  "\n"
)

cat(
  "File size (bytes):",
  file.info("results/sessionInfo.txt")$size,
  "\n"
)

# -----------------------------------------------------------------------------
# Validation note:
# Session information records the R version, operating system,
# locale, attached packages and loaded namespaces to support
# reproducibility of the analysis.
# -----------------------------------------------------------------------------