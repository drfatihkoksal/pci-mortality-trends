# Analysis data dictionary

## Persistent protected cohort

`data/private/analysis_cohort.csv` contains one row per patient's first observed
PCI. It is pseudonymized, not anonymous, and must remain in the protected
environment.

| Variable | Type | Definition |
|---|---|---|
| `study_id` | string | Sequential random-domain study identifier; no direct identifier is embedded. |
| `index_date` | date | Earliest observed therapeutic PCI report date, normalized to calendar day. |
| `index_year` | integer | Calendar year derived from `index_date`; the source `Yil` field is not used. |
| `age` | numeric | Age at PCI; values outside 18–110 are set missing pending review. |
| `sex` | category | Male/Female as recorded; missing values remain missing in models. |
| `dead` | binary | All-cause death recorded in MERNİS. |
| `death_date` | date | MERNİS death date, present only for decedents. |
| `last_vital_status_date` | date | Death date or latest verified individual MERNİS query date. |
| `followup_days` | integer | Calendar-day difference from index PCI to last vital-status date. |
| `death_before_index` | binary | Quality flag; unresolved records are excluded from outcome analyses. |
| `same_day_death` | binary | Death date equals index date; analysis time is set to 0.5 days. |
| `source_record_count` | integer | Number of source reports belonging to the patient. |
| `reports_same_day` | integer | Reports consolidated on the index calendar day. |
| `source` | category | PDF, PARS, or combined source label. |
| `clinical_presentation` | category | Preliminary rule-based report classification; not a validated EHR diagnosis. |
| `primary_pci` | binary | Report text contains a primary/primer PCI expression. |
| `left_main` | binary | Left-main expression found in the consolidated index report. |
| `graft_pci` | binary | Graft/CABG/saphenous/LIMA expression found. |
| `cto` | binary | CTO/chronic total occlusion expression found. |
| `restenosis` | binary | Restenosis/in-stent expression found. |
| `des`, `bms` | binary | Drug-eluting or bare-metal stent expression found. |
| `target_lad`, `target_rca`, `target_cx` | binary | Vessel expression found; fields are not mutually exclusive. |

Free-text fields themselves are not persisted in the protected cohort because
they may contain identifiers and are unnecessary after classification.

## Optional EHR enrichment input

`data/private/ehr_enrichment.csv` must contain one row per `patientID` and may
contain the following prespecified columns. Binary fields use 0/1; continuous
fields use the units below. Definitions and extraction code lists must be frozen
before outcome modeling.

| Variable | Unit/values | Prespecified source/window |
|---|---|---|
| `diabetes`, `hypertension`, `heart_failure` | 0/1 | Diagnosis/problem list documented before or during index admission. |
| `current_smoker` | 0/1 | Status at index admission. |
| `prior_mi`, `prior_pci`, `prior_cabg` | 0/1 | Any event/procedure before index PCI. |
| `prior_stroke`, `peripheral_artery_disease` | 0/1 | Documented before/index admission. |
| `dialysis` | 0/1 | Chronic dialysis before index PCI. |
| `egfr` | mL/min/1.73 m² | Derived consistently from the first pre-PCI admission creatinine. |
| `hemoglobin` | g/dL | First value within 24 h of admission and before PCI. |
| `lvef` | percent | Latest pre-PCI value within 180 days, or pre-PCI value during index admission. |
| `cardiogenic_shock`, `cardiac_arrest` | 0/1 | Present before PCI. |
| `ehr_presentation` | STEMI/NSTEMI/other ACS/chronic | Locked EHR hierarchy, preferred over text classification. |
| `radial_access` | 0/1 | Index procedure. |
| `multivessel_disease` | 0/1 | Pre-intervention angiographic anatomy. |
| `procedural_success` | 0/1 | Prespecified cath-lab definition. |
| `intravascular_imaging` | 0/1 | IVUS or OCT during index PCI. |

Outcome, exposure, direct identifiers, operator names, and post-discharge
medications are not imputed or included in the primary baseline-risk model.

