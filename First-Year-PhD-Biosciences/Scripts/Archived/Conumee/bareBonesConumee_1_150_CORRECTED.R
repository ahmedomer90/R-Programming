library(conumee)
library(here)

# =============================================================
# Conumee validation analysis: tumours 1 to 150
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

# Analyse tumours 1 to 150, or all available tumours if fewer than 150 exist.
tumour_indices <- seq_len(min(150L, ncol(Mset)))

cat("Tumours available: ", ncol(Mset), "\n", sep = "")
cat("Tumours selected:  ", length(tumour_indices), "\n", sep = "")
cat("Controls available:", ncol(cont), "\n\n")

# Load the reconstructed normal-cerebellum reference once.
data.c <- CNV.load(cont)


# -----------------------------
# 3. Verify probe compatibility
# -----------------------------

# A single tumour is sufficient for this structural probe-set check because
# all columns of the same MethylSet share the same probe rows.
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

# Probe membership must match. Probe order is allowed to differ because
# CNV.fit() successfully reconciles the order internally.
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
# 5. Create output directories
# -----------------------------

output_directory <- here("Figures", "Conumee")
plot_directory <- file.path(output_directory, "GenomePlots_1_150")
segment_directory <- file.path(output_directory, "Segments_1_150")

for (directory in c(output_directory, plot_directory, segment_directory)) {
  dir.create(
    directory,
    recursive = TRUE,
    showWarnings = FALSE
  )
}


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

process_one_tumour <- function(i) {
  sample_name <- colnames(Mset)[i]

  if (is.null(sample_name) || is.na(sample_name) || !nzchar(sample_name)) {
    sample_name <- paste0("tumour_", i)
  }

  # Prefixing the column number prevents accidental file-name collisions.
  sample_stem <- sprintf(
    "%03d_%s",
    i,
    safe_filename(sample_name)
  )

  plot_file <- file.path(
    plot_directory,
    paste0(sample_stem, "_genomeplot.pdf")
  )

  segments_file <- file.path(
    segment_directory,
    paste0(sample_stem, "_segments.csv")
  )

  start_time <- Sys.time()
  pdf_device <- NULL

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

      # Export segment-level results.
      segments <- CNV.write(
        x,
        what = "segments"
      )

      write.csv(
        segments,
        segments_file,
        row.names = FALSE
      )

      # Save one genome-wide plot for this tumour.
      grDevices::pdf(
        plot_file,
        width = 14,
        height = 8
      )
      pdf_device <- grDevices::dev.cur()

      CNV.genomeplot(x)

      grDevices::dev.off(pdf_device)
      pdf_device <- NULL

      stopifnot(
        file.exists(plot_file),
        file.exists(segments_file)
      )

      elapsed_minutes <- as.numeric(
        difftime(Sys.time(), start_time, units = "mins")
      )

      data.frame(
        tumour_column = i,
        sample_id = sample_name,
        status = "success",
        number_of_segments = nrow(segments),
        elapsed_minutes = elapsed_minutes,
        plot_file = plot_file,
        segments_file = segments_file,
        error_message = NA_character_,
        stringsAsFactors = FALSE
      )
    },
    error = function(e) {
      # Close only the PDF device opened by this function, if it is still open.
      if (!is.null(pdf_device) && pdf_device %in% grDevices::dev.list()) {
        try(grDevices::dev.off(pdf_device), silent = TRUE)
      }

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
# 7. Analyse tumours 1 to 150
# -----------------------------

run_results <- do.call(
  rbind,
  lapply(tumour_indices, process_one_tumour)
)


# -----------------------------
# 8. Save run summary
# -----------------------------

run_summary_file <- file.path(
  output_directory,
  "Conumee_tumours_1_150_run_summary.csv"
)

write.csv(
  run_results,
  run_summary_file,
  row.names = FALSE
)

cat("\nAnalysis complete.\n")
cat("Successful tumours: ", sum(run_results$status == "success"), "\n", sep = "")
cat("Failed tumours:     ", sum(run_results$status == "failed"), "\n", sep = "")
cat("Run summary:        ", run_summary_file, "\n", sep = "")
cat("Genome plots:       ", plot_directory, "\n", sep = "")
cat("Segment files:      ", segment_directory, "\n", sep = "")

stopifnot(file.exists(run_summary_file))

run_results
