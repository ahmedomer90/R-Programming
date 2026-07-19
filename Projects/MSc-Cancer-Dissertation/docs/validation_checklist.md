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
| 19_figure7_fragment_size_boxplot.R | ✅ | ✅ | ✅ | ✅ Validated |
| 20_figure8_sample_LOH_percentage.R | ✅ | ✅ | ✅ | ✅ Validated |
| 21_figure9_sample_category_proportions.R | ✅ | ✅ | ✅ | ✅ Validated |
| 22_figure10_chromosome_clonality_area.R | ✅ | ✅ | ✅ | ✅ Validated |
| 23_figure11_chromosome_clonality_counts_area.R | ✅ | ✅ | ✅ | ✅ Validated |
| 24_figure12_chromosome_clonality_counts_bar.R | ✅ | ✅ | ✅ | ✅ Validated |
| 25_figure13_chromosome_clonality_fragment_size_proportion.R | ✅ | ✅ | ✅ | ✅ Validated |
| 26_figure15_all_chromosome_clonality_proportions.R | ✅ | ✅ | ✅ | ✅ Validated |
| 27_figure16_all_chromosome_clonality_counts.R | ✅ | ✅ | ✅ | ✅ Validated |
| 28_figure17_all_chromosome_clonality_counts_bar.R | ✅ | ✅ | ✅ | ✅ Validated |
| 29_figure18_sample_clonality_proportions_area.R | ✅ | ✅ | ✅ | ✅ Validated |
| 30_figure19_sample_clonality_counts_area.R | ✅ | ✅ | ✅ | ✅ Validated |
| 31_figure20_sample_clonality_counts_bar.R | ✅ | ✅ | ✅ | ✅ Validated |
| 32_figure28_LOH_fragment_position_heatmap.R | ✅ | ✅ | ✅ | ✅ Validated |
| 33_figure29_category1_fragment_position_heatmap.R | ✅ | ✅ | ✅ | ✅ Validated |
| 34_figure30_category2_fragment_position_heatmap.R | ✅ | ✅ | ✅ | ✅ Validated |
| 35_figure31_category3_fragment_position_heatmap.R | ✅ | ✅ | ✅ | ✅ Validated |
| 36_figure32_category4_fragment_position_heatmap.R | ✅ | ✅ | ✅ | ✅ Validated |
| 37_figure33_category5_fragment_position_heatmap.R | ✅ | ✅ | ✅ | ✅ Validated |
| 38_figure34_category6_fragment_position_heatmap.R | ✅ | ✅ | ✅ | ✅ Validated |
| 39_figure35_category7_fragment_position_heatmap.R | ✅ | ✅ | ✅ | ✅ Validated |
| 40_figure36_category8_fragment_position_heatmap.R | ✅ | ✅ | ✅ | ✅ Validated |

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

### Script 19 – Figure 7: Fragment-Size Boxplot

- Reproduced Dissertation Figure 7 as a jittered boxplot showing total fragment sizes across LOH categories 1–8 and the non-LOH group.
- Used the validated sample-category fragment-size summary for LOH categories.
- Calculated total non-LOH fragment size for each tumour-region sample and represented these observations as the `"none"` category.
- Displayed:
  - White boxplots with black outlines
  - Semi-transparent black jittered points
  - Scientific notation on the y-axis
  - Categories ordered as `1` through `8`, followed by `"none"`
- Validation confirmed:
  - Nine groups were displayed: eight LOH categories and one non-LOH group.
  - No negative total fragment sizes were detected.
  - PNG and PDF figure files were generated successfully with non-zero file sizes.
  - The figure was visually inspected and confirmed to reproduce the expected dissertation pattern, including the substantially larger total fragment sizes in the `"none"` group.
- Outputs were saved as:
  - `figures/Figure7_fragment_size_boxplot.png`
  - `figures/Figure7_fragment_size_boxplot.pdf`

### Script 20 – Figure 8: Sample LOH Percentage

- Reproduced Dissertation Figure 8 as a stacked bar chart showing tumour-region samples in descending order of their contribution to all LOH fragments.
- Each bar represents one tumour-region sample, with coloured sections showing LOH categories 1–8.
- Validation confirmed:
  - Tumour-region samples = **434**
  - LOH categories = **8**
  - Sample-category combinations = **3,472**
  - Total LOH fragments = **16,525**
  - Sum of sample percentages = **1**
  - Highest sample contribution = **0.89%**
  - Lowest sample contribution = **0.00605%**
  - Descending sample order was confirmed.
- PNG and PDF outputs were generated successfully with non-zero file sizes.
- Outputs were saved as:
  - `figures/Figure8_sample_LOH_percentage.png`
  - `figures/Figure8_sample_LOH_percentage.pdf`
 
### Script 21 – Figure 9: Sample Category Proportions

- Reproduced Dissertation Figure 9 as a 100% stacked bar chart showing the proportional composition of LOH categories 1–8 within each tumour-region sample.
- Each sample was scaled to 100%, with coloured sections representing the proportion of fragments assigned to each LOH category.
- Samples were ordered from the smallest total LOH fragment count on the left to the largest on the right, following the figure legend.
- Validation confirmed:
  - Tumour-region samples = **434**
  - LOH categories = **8**
  - Sample-category combinations = **3,472**
  - Total LOH fragments = **16,525**
  - Category proportions summed to **1** within every sample.
  - Ascending sample-total order was confirmed.
  - Category-leading samples were identified successfully.
  - PNG and PDF outputs were generated with non-zero file sizes.
- Outputs were saved as:
  - `figures/Figure9_sample_category_proportions.png`
  - `figures/Figure9_sample_category_proportions.pdf`

### Script 22 – Figure 10: Chromosome Clonality Area Chart

- Reproduced Dissertation Figure 10 as a stacked area chart showing the proportions of clonal and subclonal LOH fragments across chromosomes.
- Displayed subclonal LOH in green and clonal LOH in red.
- Excluded the 10 rare `"none"` classifications arising from strict interval-overlap edge cases.
- Validation confirmed:
  - Chromosomes represented = **24**
  - Clonality classes plotted = **2**
  - Clonal LOH fragments = **8,073**
  - Subclonal LOH fragments = **8,442**
  - Total plotted LOH fragments = **16,515**
  - Excluded `"none"` fragments = **10**
  - Clonal and subclonal proportions summed to **1** for every chromosome.
  - PNG and PDF outputs were generated successfully with non-zero file sizes.
- Outputs were saved as:
  - `figures/Figure10_chromosome_clonality_area.png`
  - `figures/Figure10_chromosome_clonality_area.pdf`

### Script 23 – Figure 11: Chromosome Clonality Counts Area Chart

- Reproduced Dissertation Figure 11 as a stacked area chart showing the number of clonal and subclonal LOH fragments across chromosomes.
- Displayed subclonal LOH in green and clonal LOH in red.
- Excluded the 10 rare `"none"` classifications caused by strict interval-overlap edge cases.
- Validation confirmed:
  - Chromosomes represented = **24**
  - Clonality classes plotted = **2**
  - Clonal LOH fragments = **8,073**
  - Subclonal LOH fragments = **8,442**
  - Total plotted LOH fragments = **16,515**
  - Excluded `"none"` fragments = **10**
  - Maximum chromosome total = **1,246**
  - Minimum chromosome total = **282**
  - PNG and PDF outputs were generated successfully with non-zero file sizes.
- Outputs were saved as:
  - `figures/Figure11_chromosome_clonality_counts_area.png`
  - `figures/Figure11_chromosome_clonality_counts_area.pdf`
  
  ### Script 24 – Figure 12: Chromosome Clonality Counts Bar Chart

- Reproduced Dissertation Figure 12 as a stacked bar chart showing the number of clonal and subclonal LOH fragments across chromosomes.
- Displayed subclonal LOH in red and clonal LOH in blue.
- Excluded the 10 rare `"none"` classifications caused by strict interval-overlap edge cases.
- Validation confirmed:
  - Chromosomes represented = **24**
  - Clonality classes plotted = **2**
  - Clonal LOH fragments = **8,073**
  - Subclonal LOH fragments = **8,442**
  - Total plotted LOH fragments = **16,515**
  - Excluded `"none"` fragments = **10**
  - Maximum chromosome total = **1,246**
  - Minimum chromosome total = **282**
  - PNG and PDF outputs were generated successfully with non-zero file sizes.
- Outputs were saved as:
  - `figures/Figure12_chromosome_clonality_counts_bar.png`
  - `figures/Figure12_chromosome_clonality_counts_bar.pdf`

### Script 25 – Figure 13: Chromosome Clonality Fragment-Size Proportions

- Reproduced Dissertation Figure 13 as a proportional stacked area chart showing the contribution of clonal and subclonal LOH fragment size across chromosomes.
- Displayed subclonal LOH in blue and clonal LOH in salmon-red.
- Validation confirmed:
  - Chromosomes represented = **24**
  - Clonality classes plotted = **2**
  - Clonal LOH fragments = **8,073**
  - Subclonal LOH fragments = **8,442**
  - Clonal total fragment size = **152,355,287,305 bp**
  - Subclonal total fragment size = **143,025,890,323 bp**
  - Total plotted fragment size = **295,381,177,628 bp**
  - Excluded `"none"` fragments = **10**
  - Excluded `"none"` fragment size = **0 bp**
  - Clonal and subclonal fragment-size proportions summed to **1** for every chromosome.
  - PNG and PDF outputs were generated successfully with non-zero file sizes.
- Outputs were saved as:
  - `figures/Figure13_chromosome_clonality_fragment_size_proportion.png`
  - `figures/Figure13_chromosome_clonality_fragment_size_proportion.pdf`
  - Visually confirmed that subclonal LOH is displayed in blue at the bottom and clonal LOH in salmon-red at the top, matching the dissertation figure.

### Script 26 – Figure 15: Whole Genome Chromosome Clonality Proportions

- Reproduced Dissertation Figure 15 as a proportional stacked area chart showing the proportion of clonal LOH, subclonal LOH and No LOH fragments across all chromosomes in the complete TRACERx ASCAT segmentation dataset.
- Validation confirmed:
  - Chromosomes represented = **24**
  - Clonality classes plotted = **3**
  - Clonal LOH fragments = **8,087**
  - Subclonal LOH fragments = **17,816**
  - No LOH fragments = **29,612**
  - Total ASCAT fragments = **55,515**
  - Minimum chromosome fragment count = **664**
  - Maximum chromosome fragment count = **4,471**
  - Proportions summed to **1** for every chromosome.
  - PNG and PDF outputs were generated successfully with non-zero file sizes.
- Outputs were saved as:
  - `figures/Figure15_all_chromosome_clonality_proportions.png`
  - `figures/Figure15_all_chromosome_clonality_proportions.pdf`
- Visually confirmed that the stacked area order matches the dissertation:
  - Clonal LOH (salmon-red) at the top
  - No LOH (green) in the middle
  - Subclonal LOH (blue) at the bottom

### Script 27 – Figure 16: Whole-Dataset Chromosome Clonality Counts

- Reproduced Dissertation Figure 16 as a stacked area chart showing the number of fragments on each chromosome classified as clonal LOH, subclonal LOH, or No LOH.
- Validation confirmed:
  - Chromosomes represented = **24**
  - Clonality classes plotted = **3**
  - Clonal LOH fragments = **8,087**
  - Subclonal LOH fragments = **17,816**
  - No LOH fragments = **29,612**
  - Total ASCAT fragments = **55,515**
  - Minimum chromosome fragment total = **664**
  - Maximum chromosome fragment total = **4,471**
  - Negative fragment counts = **0**
  - PNG and PDF outputs were generated successfully with non-zero file sizes.
- Visually confirmed:
  - Subclonal LOH is shown in blue at the bottom.
  - No LOH is shown in green in the middle.
  - Clonal LOH is shown in salmon-red at the top.
- Outputs were saved as:
  - `figures/Figure16_all_chromosome_clonality_counts.png`
  - `figures/Figure16_all_chromosome_clonality_counts.pdf`

### Script 28 – Figure 17: Whole-Dataset Chromosome Clonality Counts Bar Chart

- Reproduced Dissertation Figure 17 as a stacked bar chart showing the number of fragments on each chromosome classified as clonal LOH, subclonal LOH, or No LOH.
- Validation confirmed:
  - Chromosomes represented = **24**
  - Clonality classes plotted = **3**
  - Clonal LOH fragments = **8,087**
  - Subclonal LOH fragments = **17,816**
  - No LOH fragments = **29,612**
  - Total ASCAT fragments = **55,515**
  - Minimum chromosome fragment total = **664**
  - Maximum chromosome fragment total = **4,471**
  - Negative fragment counts = **0**
  - PNG and PDF outputs were generated successfully with non-zero file sizes.
- Visually confirmed:
  - Subclonal LOH is shown in yellow at the bottom.
  - No LOH is shown in red in the middle.
  - Clonal LOH is shown in blue at the top.
- Outputs were saved as:
  - `figures/Figure17_all_chromosome_clonality_counts_bar.png`
  - `figures/Figure17_all_chromosome_clonality_counts_bar.pdf`

### Script 29 – Figure 18: Sample Clonality Proportions Area Chart

- Reproduced Dissertation Figure 18 as a proportional stacked area chart showing clonal and subclonal LOH composition within each tumour-region sample.
- Validation confirmed:
  - Samples represented = **434**
  - Clonality classes plotted = **2**
  - Clonal LOH fragments = **8,073**
  - Subclonal LOH fragments = **8,442**
  - Total plotted LOH fragments = **16,515**
  - Excluded `"none"` fragments = **10**
  - Minimum sample LOH total = **1**
  - Maximum sample LOH total = **147**
  - Clonal and subclonal proportions summed to **1** for every sample.
  - Negative fragment counts = **0**
  - PNG and PDF outputs were generated successfully with non-zero file sizes.
- Visually confirmed:
  - Subclonal LOH is shown in green at the bottom.
  - Clonal LOH is shown in salmon-red at the top.
  - Every sample reaches a combined height of **1**.
- Outputs were saved as:
  - `figures/Figure18_sample_clonality_proportions_area.png`
  - `figures/Figure18_sample_clonality_proportions_area.pdf`

### Script 30 – Figure 19: Sample Clonality Counts Area Chart

- Reproduced Dissertation Figure 19 as a stacked area chart showing the number of clonal and subclonal LOH fragments within each tumour-region sample.
- Validation confirmed:
  - Samples represented = **434**
  - Clonality classes plotted = **2**
  - Clonal LOH fragments = **8,073**
  - Subclonal LOH fragments = **8,442**
  - Total plotted LOH fragments = **16,515**
  - Excluded `"none"` fragments = **10**
  - Minimum sample LOH total = **1**
  - Maximum sample LOH total = **147**
  - Negative fragment counts = **0**
  - PNG and PDF outputs were generated successfully with non-zero file sizes.
- Visually confirmed:
  - Subclonal LOH is shown in green at the bottom.
  - Clonal LOH is shown in salmon-red at the top.
  - The tallest sample reaches approximately **147 fragments**.
- Outputs were saved as:
  - `figures/Figure19_sample_clonality_counts_area.png`
  - `figures/Figure19_sample_clonality_counts_area.pdf`

### Script 31 – Figure 20: Sample Clonality Counts Bar Chart

- Reproduced Dissertation Figure 20 as a stacked bar chart showing the number of clonal and subclonal LOH fragments within each tumour-region sample.
- Validation confirmed:
  - Samples represented = **434**
  - Clonality classes plotted = **2**
  - Clonal LOH fragments = **8,073**
  - Subclonal LOH fragments = **8,442**
  - Total plotted LOH fragments = **16,515**
  - Excluded `"none"` fragments = **10**
  - Minimum sample LOH total = **1**
  - Maximum sample LOH total = **147**
  - Negative fragment counts = **0**
  - PNG and PDF outputs were generated successfully with non-zero file sizes.
- Visually confirmed:
  - Subclonal LOH is shown in red at the bottom.
  - Clonal LOH is shown in blue at the top.
  - The tallest sample reaches approximately **147 fragments**.
- Outputs were saved as:
  - `figures/Figure20_sample_clonality_counts_bar.png`
  - `figures/Figure20_sample_clonality_counts_bar.pdf`

### Script 32 – Figure 28: LOH Fragment Position Heatmap

- Reproduced Dissertation Figure 28 as a chromosome-position heatmap showing the number of LOH fragments overlapping genomic intervals, regardless of LOH category or clonality.
- Chromosomes were divided into **5 Mb bins**.
- Validation confirmed:
  - Input LOH fragments = **16,525**
  - Chromosomes represented = **24**
  - Genomic bin width = **5,000,000 bp**
  - Total genomic bins = **633**
  - Bins containing at least one LOH fragment = **601**
  - Maximum overlapping fragments in one bin = **361**
  - Minimum overlapping fragments in one bin = **0**
  - Chromosomes with plotted LOH = **24**
  - Negative fragment counts = **0**
  - PNG and PDF outputs were generated successfully with non-zero file sizes.
- Outputs were saved as:
  - `figures/Figure28_LOH_fragment_position_heatmap.png`
  - `figures/Figure28_LOH_fragment_position_heatmap.pdf`

### Script 33 – Figure 29: Category 1 Fragment Position Heatmap

- Reproduced Dissertation Figure 29 as a chromosome-position heatmap showing Category 1 LOH fragments across genomic intervals.
- Category 1 represents homozygous deletion in the reconstructed LOH classification workflow.
- Chromosomes were divided into **5 Mb bins**.
- Validation confirmed:
  - Input LOH fragments = **16,525**
  - Category 1 fragments = **544**
  - Chromosomes represented = **24**
  - Genomic bin width = **5,000,000 bp**
  - Total genomic bins = **633**
  - Bins containing at least one Category 1 fragment = **202**
  - Maximum Category 1 fragments in one bin = **50**
  - Minimum positive Category 1 count in one bin = **1**
  - Chromosomes containing plotted Category 1 fragments = **24**
  - Negative fragment counts = **0**
  - PNG and PDF outputs were generated successfully with non-zero file sizes.
- Outputs were saved as:
  - `figures/Figure29_category1_fragment_position_heatmap.png`
  - `figures/Figure29_category1_fragment_position_heatmap.pdf`

### Script 34 – Figure 30: Category 2 Fragment Position Heatmap

- Reproduced Dissertation Figure 30 as a chromosome-position heatmap showing Category 2 LOH fragments across genomic intervals.
- Category 2 fragments were counted within 5 Mb genomic bins across all tumour-region samples.
- Validation confirmed:
  - Input LOH fragments = **16,525**
  - Category 2 fragments = **3,883**
  - Chromosomes represented = **24**
  - Genomic bin width = **5,000,000 bp**
  - Total genomic bins = **633**
  - Bins containing at least one Category 2 fragment = **598**
  - Maximum Category 2 fragments in one bin = **112**
  - Minimum positive Category 2 fragments in one bin = **1**
  - Chromosomes with plotted Category 2 fragments = **24**
  - Negative fragment counts = **0**
  - PNG and PDF outputs were successfully generated with non-zero file sizes.
- Outputs:
  - `figures/Figure30_category2_fragment_position_heatmap.png`
  - `figures/Figure30_category2_fragment_position_heatmap.pdf`

### Script 35 – Figure 31: Category 3 Fragment Position Heatmap

- Reproduced Dissertation Figure 31 as a chromosome-position heatmap showing Category 3 LOH fragments across genomic intervals.
- Category 3 fragments were counted within 5 Mb genomic bins across all tumour-region samples.
- Validation confirmed:
  - Input LOH fragments = **16,525**
  - Category 3 fragments = **6,796**
  - Chromosomes represented = **24**
  - Genomic bin width = **5,000,000 bp**
  - Total genomic bins = **633**
  - Bins containing at least one Category 3 fragment = **601**
  - Maximum Category 3 fragments in one bin = **170**
  - Minimum positive Category 3 fragments in one bin = **4**
  - Chromosomes with plotted Category 3 fragments = **24**
  - Negative fragment counts = **0**
  - PNG and PDF outputs were generated successfully with non-zero file sizes.
- Outputs:
  - `figures/Figure31_category3_fragment_position_heatmap.png`
  - `figures/Figure31_category3_fragment_position_heatmap.pdf`

### Script 36 – Figure 32: Category 4 Fragment Position Heatmap

- Reproduced Dissertation Figure 32 as a chromosome-position heatmap showing Category 4 LOH fragments across genomic intervals.
- Category 4 fragments were counted within 5 Mb genomic bins across all tumour-region samples.
- Validation confirmed:
  - Input LOH fragments = **16,525**
  - Category 4 fragments = **2,379**
  - Chromosomes represented = **24**
  - Genomic bin width = **5,000,000 bp**
  - Total genomic bins = **633**
  - Bins containing at least one Category 4 fragment = **596**
  - Maximum Category 4 fragments in one bin = **77**
  - Minimum positive Category 4 fragments in one bin = **1**
  - Chromosomes with plotted Category 4 fragments = **24**
  - Negative fragment counts = **0**
  - PNG and PDF outputs were generated successfully with non-zero file sizes.
- Outputs:
  - `figures/Figure32_category4_fragment_position_heatmap.png`
  - `figures/Figure32_category4_fragment_position_heatmap.pdf`

### Script 37 – Figure 33: Category 5 Fragment Position Heatmap

- Reproduced Dissertation Figure 33 as a chromosome-position heatmap showing Category 5 LOH fragments across genomic intervals.
- Category 5 fragments were counted within 5 Mb genomic bins across all tumour-region samples.
- Validation confirmed:
  - Input LOH fragments = **16,525**
  - Category 5 fragments = **1,189**
  - Chromosomes represented = **24**
  - Genomic bin width = **5,000,000 bp**
  - Total genomic bins = **633**
  - Bins containing at least one Category 5 fragment = **599**
  - Maximum Category 5 fragments in one bin = **39**
  - Minimum positive Category 5 fragments in one bin = **1**
  - Chromosomes with plotted Category 5 fragments = **24**
  - Negative fragment counts = **0**
  - PNG and PDF outputs were generated successfully with non-zero file sizes.
- Outputs:
  - `figures/Figure33_category5_fragment_position_heatmap.png`
  - `figures/Figure33_category5_fragment_position_heatmap.pdf`

### Script 38 – Figure 34: Category 6 Fragment Position Heatmap

- Reproduced Dissertation Figure 34 as a chromosome-position heatmap showing Category 6 LOH fragments across genomic intervals.
- Category 6 fragments were counted within 5 Mb genomic bins across all tumour-region samples.
- The reconstructed figure retains all 24 chromosome codes:
  - chromosomes 1–22
  - chromosome 23 = X
  - chromosome 24 = Y
- The colour legend was corrected so that the labels `2.5`, `5.0`, `7.5`, and `10.0` are displayed separately rather than overlapping.
- Validation confirmed:
  - Input LOH fragments = **16,525**
  - Category 6 fragments = **631**
  - Chromosomes represented = **24**
  - Genomic bin width = **5,000,000 bp**
  - Total genomic bins = **633**
  - Bins containing at least one Category 6 fragment = **427**
  - Maximum Category 6 fragments in one bin = **67**
  - Minimum positive Category 6 fragments in one bin = **1**
  - Chromosomes with plotted Category 6 fragments = **24**
  - Negative fragment counts = **0**
  - PNG output exists and has a non-zero file size of **114,373 bytes**
  - PDF output exists and has a non-zero file size of **8,077 bytes**
- Outputs:
  - `figures/Figure34_category6_fragment_position_heatmap.png`
  - `figures/Figure34_category6_fragment_position_heatmap.pdf`

### Script 39 – Figure 35: Category 7 Fragment Position Heatmap

- Reproduced Dissertation Figure 35 as a chromosome-position heatmap showing Category 7 LOH fragments across genomic intervals.
- Category 7 fragments were counted within 5 Mb genomic bins across all tumour-region samples.
- The complete chromosome framework includes chromosome codes 1–24:
  - chromosomes 1–22
  - chromosome 23 = X
  - chromosome 24 = Y
- Category 7 fragments were observed on **22 chromosomes**; the remaining chromosomes contain no plotted Category 7 fragments and therefore appear blank.
- Validation confirmed:
  - Input LOH fragments = **16,525**
  - Category 7 fragments = **335**
  - Chromosomes containing Category 7 fragments = **22**
  - Genomic bin width = **5,000,000 bp**
  - Total genomic bins = **633**
  - Bins containing at least one Category 7 fragment = **242**
  - Maximum Category 7 fragments in one bin = **52**
  - Minimum positive Category 7 fragments in one bin = **1**
  - Chromosomes with plotted Category 7 fragments = **22**
  - Negative fragment counts = **0**
  - PNG output exists and has a non-zero file size of **108,798 bytes**
  - PDF output exists and has a non-zero file size of **7,383 bytes**
- Outputs:
  - `figures/Figure35_category7_fragment_position_heatmap.png`
  - `figures/Figure35_category7_fragment_position_heatmap.pdf`

### Script 40 – Figure 36: Category 8 Fragment Position Heatmap

- Reproduced Dissertation Figure 36 as a chromosome-position heatmap showing Category 8 LOH fragments across genomic intervals.
- Category 8 fragments were counted within 5 Mb genomic bins across all tumour-region samples.
- The complete plotting framework includes chromosome codes 1–24:
  - chromosomes 1–22
  - chromosome 23 = X
  - chromosome 24 = Y
- Category 8 fragments were observed on **20 chromosomes**; chromosomes without Category 8 fragments remain blank in the heatmap.
- Validation confirmed:
  - Input LOH fragments = **16,525**
  - Category 8 fragments = **768**
  - Chromosomes containing Category 8 fragments = **20**
  - Genomic bin width = **5,000,000 bp**
  - Total genomic bins = **633**
  - Bins containing at least one Category 8 fragment = **186**
  - Maximum Category 8 fragments in one bin = **123**
  - Minimum positive Category 8 fragments in one bin = **1**
  - Chromosomes with plotted Category 8 fragments = **20**
  - Negative fragment counts = **0**
  - PNG output exists and has a non-zero file size of **108,477 bytes**
  - PDF output exists and has a non-zero file size of **7,200 bytes**
- Outputs:
  - `figures/Figure36_category8_fragment_position_heatmap.png`
  - `figures/Figure36_category8_fragment_position_heatmap.pdf`



