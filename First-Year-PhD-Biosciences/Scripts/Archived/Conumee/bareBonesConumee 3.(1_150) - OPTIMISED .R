library(conumee)
library(here)

# -----------------------------
# 1. Load one tumour sample
# -----------------------------

tumour_file <- here("Data", "Mset.GM.RDS")

stopifnot(file.exists(tumour_file))

Mset <- readRDS(tumour_file)

# Select one tumour for the initial validation run
MSet_test <- Mset[, 1, drop = FALSE]

data.q <- CNV.load(MSet_test)
names(data.q) <- colnames(MSet_test)

data.q


# -----------------------------
# 2. Load reconstructed controls
# -----------------------------

control_file <- here(
  "Data",
  "mSet_Cerebellum_GSE54880.rds"
)

stopifnot(file.exists(control_file))

cont <- readRDS(control_file)
data.c <- CNV.load(cont)

query_probes <- rownames(data.q@intensity)
control_probes <- rownames(data.c@intensity)

cat("Query probes:   ", length(query_probes), "\n")
cat("Control probes: ", length(control_probes), "\n")

cat("Same probe set? ", setequal(query_probes, control_probes), "\n")
cat("Same order?     ", identical(query_probes, control_probes), "\n")


data.c


# -----------------------------
# 3. Verify probe compatibility
# -----------------------------

### THIS PART IS DELETED! ###

# -----------------------------
# 4. Create 450k annotation
# -----------------------------

data(exclude_regions, package = "conumee")
data(detail_regions, package = "conumee")

anno <- CNV.create_anno(
  array_type = "450k",
  exclude_regions = exclude_regions,
  detail_regions = detail_regions
)

anno


# -----------------------------
# 5. Run Conumee analysis
# -----------------------------

x <- CNV.fit(data.q, data.c, anno)
x <- CNV.bin(x)
x <- CNV.detail(x)
x <- CNV.segment(x)

# -----------------------------
# Save segmentation results
# -----------------------------
segments <- CNV.write(x, what = "segments")

segments_file <- here(
  "Figures",
  "Conumee",
  "TestOutput_one_tumour_segments.csv"
)

write.csv(
  segments,
  segments_file,
  row.names = FALSE
)

pdf(plot_file, width = 14, height = 8)

CNV.genomeplot(x)

dev.off()

# -----------------------------
# 6. Save plot
# -----------------------------

output_directory <- here("Figures", "Conumee")
dir.create(
  output_directory,
  recursive = TRUE,
  showWarnings = FALSE
)

plot_file <- here(
  "Figures",
  "Conumee",
  "TestOutput_one_tumour.pdf"
)

pdf(plot_file, width = 14, height = 8)

CNV.genomeplot(x)

dev.off()

file.exists(plot_file)
file.info(plot_file)[, c("size", "mtime")]


# -----------------------------
# 7. Export segments
# -----------------------------

segments <- CNV.write(
  x,
  what = "segments"
)

segments_file <- here(
  "Figures",
  "Conumee",
  "TestOutput_one_tumour_segments.csv"
)

write.csv(
  segments,
  segments_file,
  row.names = FALSE
)


# -----------------------------
# 8. Confirm outputs
# -----------------------------

file.exists(plot_file)
file.exists(segments_file)

plot_file
segments_file