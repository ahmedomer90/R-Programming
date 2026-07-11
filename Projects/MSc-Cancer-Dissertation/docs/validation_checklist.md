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
| 10_... | ⬜ | ⬜ | ⬜ | Pending |
| 11_... | ⬜ | ⬜ | ⬜ | Pending |
| 12_... | ⬜ | ⬜ | ⬜ | Pending |
| 13_... | ⬜ | ⬜ | ⬜ | Pending |
| 14_... | ⬜ | ⬜ | ⬜ | Pending |
| 15_... | ⬜ | ⬜ | ⬜ | Pending |

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

- 
