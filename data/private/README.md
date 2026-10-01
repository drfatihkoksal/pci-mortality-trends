# Protected data area

Files in this directory are excluded from version control.

Optional `ehr_enrichment.csv` must contain one row for each patient's first
observed PCI and must include `patientID` plus the prespecified fields documented
in `docs/data_dictionary.md`. A direct national identifier may be used only for
the in-memory linkage step and must not be copied to generated outputs.

The pipeline creates `analysis_cohort.csv` here. It contains no name, national
identifier, catheter number, or physician name, but it remains patient-level
health data and must be protected accordingly.

