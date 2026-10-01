# Statistical analysis plan

Version: 1.0  
Status: implemented for the available workbook; EHR enrichment pending  
Analysis platform: R with Python used only for secure Excel ingestion

## Objective and estimand

The primary objective is to estimate the association between index calendar year
and 1-year all-cause mortality after the first observed PCI, standardized to the
pooled 2015–2024 cohort. The principal contrast is the adjusted absolute risk
difference between 2024 and 2015; the joint nonlinear calendar-year test is also
reported. Neither contrast is interpreted causally.

## Cohorts and outcomes

- Primary: first observed PCI during 2015–2024; all-cause death by 365 days.
- Early: first observed PCI during 2015–2025; all-cause death by 30 days.
- Long-term: first observed PCI during 2015–2020; all-cause death by 1,825 days.
- Same-day deaths receive an analysis time of 0.5 days.
- A patient is mature for a horizon if vital status is known beyond the horizon
  or death occurred within it.
- Unresolved death-before-index records are excluded.

## Models

Fixed-horizon logistic regression models index year with a 3-df natural spline.
The main scientific output is the annual marginal standardized mortality risk,
not the conditional odds ratio. Confidence intervals use the model covariance
matrix and delta method. A joint likelihood-ratio test compares the spline model
with an otherwise identical model omitting calendar year.

Three model tiers are always produced:

1. Unadjusted calendar-year model.
2. Age/sex model; this is the selected preliminary model while EHR enrichment is
   absent.
3. Available case-mix model; it becomes selected only after the prespecified EHR
   fields are supplied and pass completeness/validation checks.

The EHR-enriched baseline model includes nonlinear age, presentation, shock,
cardiac arrest, diabetes, hypertension, smoking, heart failure, prior vascular
history, dialysis, nonlinear eGFR, hemoglobin, and LVEF when at least 60% observed.
Angiographic and procedural variables are treated as exploratory additions, not
as the primary baseline adjustment set.

A Cox model censored at each fixed horizon is a secondary time-to-event analysis;
Schoenfeld diagnostics are exported. Kaplan–Meier survival is shown for the
mature 2015–2017 and 2018–2020 cohorts.

## Missing data

The current limited model uses complete observations because age/sex missingness
is minimal. Once EHR enrichment is supplied, 40-fold multiple imputation is
required before the enriched model is designated final. The imputation model
will include calendar year, outcome, follow-up/Nelson–Aalen information, and all
analysis covariates. Exposure and outcome are never imputed. A core covariate
remaining more than 40% missing after data recovery is excluded from the primary
model and analyzed in an enriched-subset sensitivity analysis.

## Prespecified sensitivity and subgroup analyses

- Complete-case and age/sex-only estimates.
- Text-derived case-mix model, explicitly labeled preliminary.
- Separate PCI-admission analysis with patient-clustered standard errors after
  admission identifiers are supplied.
- Fixed-horizon modified Poisson model as a check of the logistic standardization.
- Sex, age ≥75 years, and ACS/chronic-presentation interaction analyses; all are
  exploratory and interpreted without claims of multiplicity-controlled proof.
- PDF/PARS bridge validation in 2022 and source-stratified text-classifier metrics.

## Reporting

All estimates use two-sided 95% confidence intervals. Baseline tables emphasize
standardized differences rather than significance tests. Exact analysis counts,
event counts, model formulas, session information, and diagnostics are generated
by code. Reporting follows STROBE and RECORD.

## Freeze conditions before submission

- Ethics approval entered (2024-TBEK 2026/05-22, 20 May 2026); confirm waiver wording.
- EHR enrichment either completed with multiple imputation or formally removed
  from the intended submission scope.
- The 53 death-before-index records are excluded from all analyses (decision 2026-10-01).
- Text extraction validated on the prespecified clinician-reviewed sample.
- Final cohort hash and key result table signed off by the clinical lead and an
  independent statistician.


## Amendment 1 (2026-10-01): post hoc sensitivity analyses

Added after review of the preliminary results; all are labelled exploratory in
the manuscript and do not replace the prespecified estimand.

- Descriptive peak contrasts (peak year vs first year; final year vs peak year),
  because the modelled calendar-year curve was non-monotonic.
- Calendar-year specification: natural splines with 2, 4, and 5 df, and year as
  a categorical variable.
- Exclusion of the first index year (2015).
- Report-source adjustment (PDF vs PARS, switch in September 2022), with risks
  standardized to PDF documentation, and a 2022 within-year bridge comparison.
- Fixed 2-year look-back episode cohort (2017 onward): every procedure day with
  no recorded PCI in the preceding 730 days; patient-clustered sandwich variance
  and robust Wald test. Addresses look-back that lengthens with calendar time
  when the first observed procedure is selected (left truncation).
- Calendar-period mortality among 1-year survivors of 2015–2018 procedures:
  Poisson model by calendar year (2016–2025), adjusted for attained age and sex.
- Cox models now use the same covariates as the selected logistic model, with a
  likelihood-ratio test for calendar year.
