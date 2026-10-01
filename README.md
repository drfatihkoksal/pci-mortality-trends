# Temporal trends in all-cause mortality after PCI, 2015–2024

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23083746.svg)](https://doi.org/10.5281/zenodo.23083746)

Analysis code, aggregate results, and documentation for:

> Köksal F, Levent F, Koca F, Severgün K, Doğan T, Kapsız M, Melek M,
> Vatansever Ağca F, Tenekecioğlu E, Arı H. *Temporal Trends in One-Year
> All-Cause Mortality After Percutaneous Coronary Intervention: A Decade-Long
> Cohort Linked to Civil Registry Data.* (manuscript in preparation)

Department of Cardiology, University of Health Sciences, Bursa Yüksek İhtisas
Training and Research Hospital, Bursa, Türkiye. Ethics approval: Clinical
Research Ethics Committee of the same hospital, decision no. 2024-TBEK
2026/05-22 (20 May 2026).

## What is and is not in this repository

| Included | Not included |
| --- | --- |
| Full analysis pipeline (`analysis/`) | Source PCI report workbook |
| Validation tests (`tests/`) | Patient-level analysis files (`data/private/`) |
| Aggregate tables, figures, model results (`outputs/`) | Manuscript draft |
| Statistical analysis plan, data dictionary, STROBE/RECORD checklist (`docs/`) | |

Patient-level data are protected health information linked to the national
civil registry (MERNİS) and cannot be shared publicly. Every file in `outputs/`
is an aggregate summary; no file contains names, national identifiers, catheter
numbers, operator names, or study identifiers.

## Pipeline

```bash
bash analysis/run_pipeline.sh
```

1. `01_extract_deidentify.py` reads the source workbook, consolidates same-day
   reports, selects each patient's first observed PCI, and writes de-identified
   patient-level files to `data/private/` (not version-controlled). Direct
   identifiers are used only in memory.
2. `02_analyze.R` fits the fixed-horizon logistic models with natural-spline
   calendar year, marginal standardization with delta-method confidence
   intervals, Cox models, all sensitivity analyses (year specification,
   exclusion of 2015, report-source adjustment, fixed 2-year look-back episodes
   with patient-clustered robust variance), and the calendar-period Poisson
   analysis. It writes everything in `outputs/`.
3. `03_build_manuscript.py` generates the manuscript, supplement, and
   submission documents from `outputs/`.
4. `tests/test_pipeline.py` checks de-identification, cohort reconciliation, and
   agreement between tables and reported results.

Steps 1, 2, and 4 require the protected source data. The aggregate outputs can
be inspected without running anything.

Requirements: Python 3 (`requirements.txt`), R ≥ 4.3 with the `survival`
package, and Pandoc (LibreOffice optional for PDF). The R session used for the
published results is recorded in `outputs/audit/r_session_info.txt`.

## Key outputs

| File | Content |
| --- | --- |
| `outputs/tables/annual_adjusted_mortality.csv` | Standardized annual risks, all models |
| `outputs/tables/model_contrasts.csv` | First- vs final-year contrasts by adjustment tier |
| `outputs/tables/sensitivity_analyses.csv` | Peak and net contrasts for every sensitivity analysis |
| `outputs/tables/source_bridge_2022.csv` | PDF vs PARS outcomes in the 2022 transition year |
| `outputs/tables/period_mortality_2015_2018_cohorts.csv` | Calendar-period death rates among 1-year survivors |
| `outputs/results/cox_summary.csv` | Secondary Cox models and proportional-hazards tests |
| `outputs/audit/followup_maturity_by_year.csv` | Vital-status completeness by horizon |

## Licence

Code is released under the MIT licence (`LICENSE`). Aggregate data, figures,
and documentation are released under CC BY 4.0. Please cite the article and this
repository when reusing them: version 1.0.0, https://doi.org/10.5281/zenodo.23083747
(all versions: https://doi.org/10.5281/zenodo.23083746; see `CITATION.cff`).

## Contact

Fatih Köksal, MD — dr.fatihkoksal@hotmail.com — ORCID 0000-0002-4197-4683
