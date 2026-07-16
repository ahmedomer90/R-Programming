# Script Validation Checklist

This document records the validation status of each script in the MSc Cancer Dissertation analysis pipeline.

Validation criteria:

- **Compared with Appendix A** – Script reviewed against the original MSc dissertation code.
- **Runs Successfully** – Script executed successfully from a clean R session.
- **Outputs Verified** – Expected output files and summary checks confirmed.

| Script | Compared with Appendix A | Runs Successfully | Outputs Verified | Status |
|:-------|:------------------------:|:-----------------:|:----------------:|:------:|
| 01_data_import.R | ✅ | ✅ | ✅ | ✅ Validated |
| 02_classify_LOH.R | ✅ | ✅ | ✅ | ✅ Validated |
| 03_LOH_clonality.R | ✅ | ✅ | ✅ | ✅ Validated |
| 04_LOH_categories.R | ✅ | ✅ | ✅ | ✅ Validated |
| 05_category_clonality_summary.R | ✅ | ✅ | ✅ | ✅ Validated |
| 06_cochran_armitage_test.R | ✅ | ✅ | ✅ | ✅ Validated |
| 07_chromosome_clonality_summary.R | ✅ | ✅ | ✅ | ✅ Validated |
| 08_sample_clonality_summary.R | ✅ | ✅ | ✅ | ✅ Validated |
| 09_fragment_size_analysis.R | ✅ | ✅ | ✅ | ✅ Validated |
| 10_sample_category_summary.R | ✅ | ✅ | ✅ | ✅ Validated |
| 11_sample_category_fragment_size_summary.R | ✅ | ✅ | ✅ | ✅ Validated |
| 12_high_prevalence_LOH_regions.R | ✅ | ✅ | ✅ | ✅ Validated |
| 13_gene_annotation.R | ✅ | ✅ | ✅ | ✅ Validated |
| 14_cancer_gene_filtering.R | ✅ | ✅ | ✅ | ✅ Validated |
| 15_visualisation_LOH_patterns.R | ✅ | ✅ | ✅ | ✅ Validated |
| 16_export_session_info.R | ✅ | ✅ | ✅ | ✅ Validated |
| 17_figure4_major_allele_distribution.R | ✅ | ✅ | ✅ | ✅ Validated |
| 18_figure6_sample_category_heatmap.R | ✅ | ✅ | ✅ | ✅ Validated |

---

## Validation Notes

### Script 01 – Data Import

- Imported the original TRACERx ASCAT segment data.
- Imported the original clinical metadata.
- Verified expected dimensions and column names.
- Successfully created:
  - `tracerx_ascat_seg_imported.rds`
  - `clinical_data_imported.rds`

---

### Script 02 – LOH Classification

- LOH classified using the original MSc criterion (`nMinor == 0`).
- Validation summary:
  - LOH fragments: **16,525**
  - Non-LOH fragments: **38,990**
  - Total fragments: **55,515**

---

### Script 03 – LOH Clonality

- Original overlap-based clonality algorithm reproduced from Appendix A.
- Validation summary:
  - Clonal: **8,073**
  - Subclonal: **8,442**
  - None: **10**
- Investigated the 10 `"none"` cases.
- Confirmed they correspond to rare single-position genomic intervals (`startpos == endpos`) retained from the original MSc methodology.

---

### Script 04 – LOH Categories

- Categories generated from `nMajor` using the original MSc classification:
  - Category 1 = nMajor 0
  - ...
  - Category 8 = nMajor ≥ 7
- Validation confirmed:
  - All **16,525** LOH fragments assigned to categories.
  - All **38,990** non-LOH fragments assigned `"none"`.

---

### Script 05 – Category Clonality Summary

- Summarised clonality counts and proportions within each LOH category.
- Validation confirmed:
  - Total LOH fragments = **16,525**
  - Category proportions sum to **1** for every category.
  - Output files successfully generated.

### Script 06 – Cochran-Armitage Trend Test

- Tested for a trend in clonality across the eight ordered LOH categories.
- Excluded the 10 `"none"` cases because the test requires a two-column outcome table.
- Validation result:
  - Z statistic: **28.885**
  - p-value: **1.853213 × 10⁻183**
- Confirmed that the test input table and result files were generated successfully.

### Script 07 – Chromosome Clonality Summary

- Summarised clonality counts and proportions for each chromosome.
- Analysed both the complete ASCAT dataset and the LOH-only subset.
- Validation confirmed:
  - Total complete fragments = **55,515**
  - Total LOH fragments = **16,525**
  - Proportions sum to **1** for all 24 chromosomes.
- Output summary files were generated successfully.

### Script 08 – Sample Clonality Summary

- Summarised LOH clonality counts and proportions for each tumour-region sample.
- Created separate datasets containing:
  - Clonal LOH fragments
  - Subclonal LOH fragments
  - Rare `"none"` clonality fragments
- Validation confirmed:
  - Total LOH fragments = **16,525**
  - Clonal fragments = **8,073**
  - Subclonal fragments = **8,442**
  - `"None"` fragments = **10**
  - Total tumour-region samples = **434**
  - Clonality proportions sum to **1** for every sample.
- Processed datasets and summary tables were generated successfully.

### Script 09 – Fragment Size Analysis

- Calculated LOH fragment sizes using the original MSc method:
  - `fragment_size_bp = endpos - startpos`
- Validation confirmed:
  - Total LOH fragments analysed = **16,525**
  - Negative fragment sizes = **0**
  - Zero-length fragments (`startpos == endpos`) = **339**
- Category and category-by-clonality summaries accounted for all **16,525** LOH fragments.
- Processed datasets and summary tables were generated successfully.

### Script 10 – Sample Category Summary

- Summarised the number of LOH fragments in each category for every tumour-region sample.
- Included zero-count sample-category combinations.
- Validation confirmed:
  - Distinct tumour-region samples = **434**
  - Sample-category combinations = **3,472**
  - Total LOH fragments = **16,525**
  - Category totals matched the validated results from Script 4.
- Processed and CSV output files were generated successfully.

### Script 11 – Sample-Category Fragment Size Summary

- Summarised the total LOH fragment size for every tumour-region sample and LOH category.
- Included missing sample-category combinations with zero fragment counts and zero total fragment size.
- Validation confirmed:
  - Distinct tumour-region samples = **434**
  - Sample-category combinations = **3,472**
  - Total LOH fragment count = **16,525**
  - Total LOH fragment size = **295,381,177,628 bp**
  - Zero-total sample-category combinations = **1,303**
- The completed summary reproduced the exact total fragment size from the validated LOH dataset.
- Processed RDS and CSV output files were generated successfully.

### Script 12 – High-Prevalence LOH Regions

- Grouped LOH fragments by chromosome, category, start position and end position.
- Counted the number of distinct tumour-region samples containing each identical genomic interval.
- Validation confirmed:
  - Unique LOH region-category combinations = **13,193**
  - Maximum sample prevalence for a region = **59**
  - Regions present in at least two samples = **1,610**
- The processed RDS and CSV output files were generated successfully.

### Script 13 – Gene Annotation

- Annotated high-prevalence LOH regions with overlapping genes from Ensembl BioMart.
- Used chromosome-level batched Ensembl queries and calculated interval overlaps locally with `IRanges`.
- Converted ASCAT chromosome numbers 23 and 24 to Ensembl chromosomes `X` and `Y`.
- Standardised chromosome data types and used explicit `dplyr::select()` to avoid namespace conflicts.
- Validation confirmed:
  - Input LOH regions = **13,193**
  - Regions with at least one HGNC gene = **12,694**
  - Regions without an HGNC gene = **499**
  - Annotated region–gene records = **3,171,328**
  - Unique HGNC gene symbols = **40,935**
  - Annotated and unannotated region counts together reproduced all **13,193** input regions.
- Processed RDS and CSV output files were generated successfully.

### Script 14 – Cancer Gene Filtering

- Filtered LOH-affected annotated genes using the curated `driverGenes_120516_unique.txt` reference list.
- Added Cancer Gene Census annotations where available.
- Retained the full region–gene table and created separate one-row-per-gene and gene-level summary outputs.
- Validation confirmed:
  - Curated driver genes = **641**
  - Unique Cancer Gene Census genes = **723**
  - Matched LOH driver-gene records = **47,940**
  - Unique matched LOH driver genes = **604**
  - Matched genes present in the Cancer Gene Census = **547**
  - Matched genes absent from the Cancer Gene Census = **57**
- All expected processed RDS and CSV outputs were generated successfully.

### Script 15 – Visualisation of LOH Patterns

- Generated publication-quality figures illustrating:
  - LOH category distribution
  - Clonality proportions by LOH category
  - Clonality proportions by chromosome
  - Fragment-size distributions by LOH category
  - Clonality proportions across tumour-region samples
- Validation confirmed:
  - All five expected PNG files were generated successfully.
  - All figure files had non-zero file sizes.
  - Figures were visually inspected and confirmed to display the expected LOH patterns.
- Figures were saved in the `figures/` directory for downstream reporting and publication.

### Script 16 – Export Session Information

- Exported the complete R session information to support computational reproducibility.
- Validation confirmed:
  - `results/sessionInfo.txt` was created successfully.
  - Output file size = **1,419 bytes**.
  - The file records the R version, operating system, platform, locale, time zone, attached packages and loaded namespaces.
- The exported environment included:
  - R **4.5.3**
  - Apple Silicon platform
  - macOS Tahoe 26.3
  - Europe/London time zone
- This script completes the reproducible analysis pipeline.

### Script 17 – Figure 4: Major-Allele Copy-Number Distributions

- Reproduced Dissertation Figure 4 as a three-panel histogram figure.
- Generated:
  - **Panel A:** Major-allele copy-number distribution in LOH fragments.
  - **Panel B:** Major-allele copy-number distribution in all ASCAT fragments.
  - **Panel C:** Major-allele copy-number distribution in non-LOH fragments.
- Used blue histogram bars with black outlines and the original grey `ggplot2` background to closely match the dissertation figure.
- Reproduced the dissertation layout:
  - Panel A centred at the top.
  - Panel B positioned at the bottom left.
  - Panel C positioned at the bottom right.
- Validation confirmed:
  - LOH fragments represented in Panel A = **16,525**
  - All fragments represented in Panel B = **55,515**
  - Non-LOH fragments represented in Panel C = **38,990**
  - `nMajor` range in LOH fragments = **0–81**
  - `nMajor` range in all fragments = **0–97**
  - `nMajor` range in non-LOH fragments = **0–97**
  - PNG and PDF figure files were generated successfully with non-zero file sizes.
  - Panel labels **A**, **B**, and **C** were visually inspected and confirmed to be positioned correctly.
- Outputs were saved as:
  - `figures/Figure4_major_allele_distribution.png`
  - `figures/Figure4_major_allele_distribution.pdf`
 
### Script 18 – Figure 6: Sample–Category LOH Heatmap

- Reproduced Dissertation Figure 6 as a heatmap showing the number of LOH fragments in each category for every tumour-region sample.
- Used the validated sample-category summary produced by Script 10.
- Validation confirmed:
  - Tumour-region samples = **434**
  - LOH categories = **8**
  - Sample-category combinations = **3,472**
  - Total LOH fragments represented = **16,525**
  - PNG and PDF figure files were generated successfully with non-zero file sizes.
- Low fragment counts are displayed in light blue and higher counts in dark blue.
