library(conumee)
library(here)

# =============================================================
# Conumee validation analysis: first 5 tumours
# One multi-page PDF: one tumour genome plot per page
# =============================================================

# -----------------------------
# 1. Define input files
# -----------------------------

tumour_file <- here("Data", "Mset.GM.RDS")
control_file <- here("Data", "mSet_Cerebellum_GSE54880.rds")

stopifnot(
  file.exists(tumour_file),
  file.exists(control_file)
)


# -----------------------------
# 2. Load tumour and control data
# -----------------------------

Mset <- readRDS(tumour_file)
cont <- readRDS(control_file)

stopifnot(
  ncol(Mset) >= 1,
  ncol(cont) >= 1
)

# Change this value later if you want to analyse more tumours, for example 150L.

# Select the tumour columns to analyse.
tumour_indices <- 601:763

stopifnot(
  length(tumour_indices) > 0,
  min(tumour_indices) >= 1,
  max(tumour_indices) <= ncol(Mset)
)

# Label used for output filenames.
range_label <- paste0(
  min(tumour_indices),
  "_",
  max(tumour_indices)
)

cat("Tumours available: ", ncol(Mset), "\n", sep = "")
cat("Tumours selected:  ", length(tumour_indices), "\n", sep = "")
cat("Controls available: ", ncol(cont), "\n\n", sep = "")

# Load the reconstructed normal-cerebellum reference once.
data.c <- CNV.load(cont)


# -----------------------------
# 3. Verify probe compatibility
# -----------------------------

# One tumour is sufficient for this structural check because all columns
# of the same MethylSet share the same probe rows.
probe_test_mset <- Mset[, tumour_indices[1], drop = FALSE]
probe_test_query <- CNV.load(probe_test_mset)

query_probes <- rownames(probe_test_query@intensity)
control_probes <- rownames(data.c@intensity)

same_probe_set <- setequal(query_probes, control_probes)
same_probe_order <- identical(query_probes, control_probes)

cat("Query probes:   ", length(query_probes), "\n", sep = "")
cat("Control probes: ", length(control_probes), "\n", sep = "")
cat("Same probe set? ", same_probe_set, "\n", sep = "")
cat("Same order?     ", same_probe_order, "\n\n", sep = "")

# Probe membership must match. Probe order may differ because CNV.fit()
# has already been shown to reconcile the ordering successfully.
stopifnot(same_probe_set)

rm(probe_test_mset, probe_test_query, query_probes, control_probes)


# -----------------------------
# 4. Create the 450k annotation
# -----------------------------

data(exclude_regions, package = "conumee")
data(detail_regions, package = "conumee")

anno <- CNV.create_anno(
  array_type = "450k",
  exclude_regions = exclude_regions,
  detail_regions = detail_regions
)


# -----------------------------
# 5. Create output locations
# -----------------------------

output_directory <- here("Figures", "Conumee")
segment_directory <- file.path(
  output_directory,
  paste0("Segments_", range_label)
)

combined_plot_file <- file.path(
  output_directory,
  paste0("GenomePlots_", range_label, ".pdf")
)

run_summary_file <- file.path(
  output_directory,
  paste0("Conumee_tumours_", range_label, "_run_summary.csv")
)

dir.create(
  output_directory,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  segment_directory,
  recursive = TRUE,
  showWarnings = FALSE
)


# -----------------------------
# 6. Helper functions
# -----------------------------

safe_filename <- function(x) {
  x <- as.character(x)
  x <- gsub("[^A-Za-z0-9._-]+", "_", x)
  x <- gsub("_+", "_", x)
  x <- sub("^_+", "", x)
  x <- sub("_+$", "", x)

  if (!nzchar(x)) {
    x <- "unnamed_sample"
  }

  x
}

is_graphics_device_open <- function(device_number) {
  open_devices <- grDevices::dev.list()
  !is.null(open_devices) && device_number %in% open_devices
}

process_one_tumour <- function(i) {
  sample_name <- colnames(Mset)[i]

  if (is.null(sample_name) || is.na(sample_name) || !nzchar(sample_name)) {
    sample_name <- paste0("tumour_", i)
  }

  # Prefixing the column number prevents accidental filename collisions.
  sample_stem <- sprintf(
    "%03d_%s",
    i,
    safe_filename(sample_name)
  )

  segments_file <- file.path(
    segment_directory,
    paste0(sample_stem, "_segments.csv")
  )

  start_time <- Sys.time()

  result <- tryCatch(
    {
      cat(
        sprintf(
          "[%d/%d] Processing tumour column %d: %s\n",
          match(i, tumour_indices),
          length(tumour_indices),
          i,
          sample_name
        )
      )

      MSet_test <- Mset[, i, drop = FALSE]
      data.q <- CNV.load(MSet_test)
      names(data.q) <- colnames(MSet_test)

      x <- CNV.fit(data.q, data.c, anno)
      x <- CNV.bin(x)
      x <- CNV.detail(x)
      x <- CNV.segment(x)

      # Export this tumour's segment-level results.
      segments <- CNV.write(
        x,
        what = "segments"
      )

      write.csv(
        segments,
        segments_file,
        row.names = FALSE
      )

      stopifnot(file.exists(segments_file))

      # Add this tumour's genome plot as the next page of the shared PDF.
      CNV.genomeplot(x)

      elapsed_minutes <- as.numeric(
        difftime(Sys.time(), start_time, units = "mins")
      )

      data.frame(
        tumour_column = i,
        sample_id = sample_name,
        status = "success",
        number_of_segments = nrow(segments),
        elapsed_minutes = elapsed_minutes,
        plot_file = combined_plot_file,
        segments_file = segments_file,
        error_message = NA_character_,
        stringsAsFactors = FALSE
      )
    },
    error = function(e) {
      elapsed_minutes <- as.numeric(
        difftime(Sys.time(), start_time, units = "mins")
      )

      message(
        sprintf(
          "Tumour column %d (%s) failed: %s",
          i,
          sample_name,
          conditionMessage(e)
        )
      )

      data.frame(
        tumour_column = i,
        sample_id = sample_name,
        status = "failed",
        number_of_segments = NA_integer_,
        elapsed_minutes = elapsed_minutes,
        plot_file = NA_character_,
        segments_file = NA_character_,
        error_message = conditionMessage(e),
        stringsAsFactors = FALSE
      )
    }
  )

  result
}


# -----------------------------
# 7. Analyse selected tumours
# -----------------------------

# =============================================================
# Create the combined PDF locally first
# =============================================================

# Writing an actively changing multi-page PDF directly into iCloud can fail.
# Therefore, create it in R's local temporary directory first.
temporary_plot_file <- tempfile(
  pattern = paste0("GenomePlots_", range_label, "_"),
  fileext = ".pdf"
)

cat("Temporary PDF: ", temporary_plot_file, "\n", sep = "")
cat("Final PDF:     ", combined_plot_file, "\n\n", sep = "")

# Close any graphics devices left open by an earlier failed run.
graphics.off()

# Open one local multi-page PDF before processing the selected tumours.
grDevices::pdf(
  file = temporary_plot_file,
  width = 14,
  height = 8,
  onefile = TRUE
)

pdf_device <- grDevices::dev.cur()

run_results <- tryCatch(
  {
    do.call(
      rbind,
      lapply(tumour_indices, process_one_tumour)
    )
  },
  finally = {
    # Always attempt to close the shared PDF device.
    if (is_graphics_device_open(pdf_device)) {
      close_result <- try(
        grDevices::dev.off(pdf_device),
        silent = TRUE
      )
      
      if (inherits(close_result, "try-error")) {
        warning(
          "The temporary PDF device could not be closed normally: ",
          as.character(close_result)
        )
      }
    }
  }
)

# Confirm that R successfully created the temporary PDF.
if (!file.exists(temporary_plot_file)) {
  stop("The temporary multi-page PDF was not created.")
}

if (file.info(temporary_plot_file)$size == 0) {
  stop("The temporary multi-page PDF is empty.")
}

# Remove the previous final PDF, if present.
if (file.exists(combined_plot_file)) {
  removed <- file.remove(combined_plot_file)
  
  if (!removed) {
    stop(
      "The old final PDF could not be removed. ",
      "Close it in Preview or another application and try again."
    )
  }
}

# Copy the fully completed PDF into the iCloud project folder.
copied <- file.copy(
  from = temporary_plot_file,
  to = combined_plot_file,
  overwrite = TRUE
)

if (!copied) {
  stop(
    "The PDF was created locally, but it could not be copied to: ",
    combined_plot_file
  )
}

stopifnot(
  file.exists(combined_plot_file),
  file.info(combined_plot_file)$size > 0
)

cat(
  "Combined PDF copied successfully to:\n",
  combined_plot_file,
  "\n\n",
  sep = ""
)


# -----------------------------
# 8. Save run summary
# -----------------------------

write.csv(
  run_results,
  run_summary_file,
  row.names = FALSE
)

stopifnot(file.exists(run_summary_file))

cat("\nAnalysis complete.\n")
cat("Successful tumours: ", sum(run_results$status == "success"), "\n", sep = "")
cat("Failed tumours:     ", sum(run_results$status == "failed"), "\n", sep = "")
cat("Combined PDF:       ", combined_plot_file, "\n", sep = "")
cat("Run summary:        ", run_summary_file, "\n", sep = "")
cat("Segment files:      ", segment_directory, "\n", sep = "")

run_results
