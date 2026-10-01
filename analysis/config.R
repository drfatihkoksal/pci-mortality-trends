study_config <- list(
  input_csv = "data/private/analysis_cohort.csv",
  results_dir = "outputs/results",
  tables_dir = "outputs/tables",
  figures_dir = "outputs/figures",
  audit_dir = "outputs/audit",
  primary_start_year = 2015L,
  primary_end_year = 2024L,
  early_end_year = 2025L,
  long_term_end_year = 2020L,
  primary_horizon_days = 365L,
  early_horizon_days = 30L,
  long_term_horizon_days = 1825L,
  manuscript_title = paste(
    "Temporal Trends in One-Year All-Cause Mortality After",
    "Percutaneous Coronary Intervention: A Decade-Long Cohort",
    "Linked to Civil Registry Data"
  )
)

