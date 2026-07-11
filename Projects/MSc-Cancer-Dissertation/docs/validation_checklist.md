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
| 06_... | ⬜ | ⬜ | ⬜ | Pending |
| 07_... | ⬜ | ⬜ | ⬜ | Pending |
| 08_... | ⬜ | ⬜ | ⬜ | Pending |
| 09_... | ⬜ | ⬜ | ⬜ | Pending |
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
