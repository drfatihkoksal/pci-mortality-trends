#!/usr/bin/env python3
"""Generate the manuscript, supplement, and submission documents from outputs."""

from __future__ import annotations

import csv
import json
import re
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "outputs" / "results"
TABLES = ROOT / "outputs" / "tables"
AUDIT = ROOT / "outputs" / "audit"
MANUSCRIPT = ROOT / "manuscript"
SUBMISSION = ROOT / "submission"


def key_values() -> dict[str, str]:
    with (RESULTS / "key_results.csv").open(encoding="utf-8") as handle:
        return {row["key"]: row["value"] for row in csv.DictReader(handle)}


def fnum(value: str | float, digits: int = 1) -> str:
    return f"{float(value):,.{digits}f}"


def fint(value: str | float) -> str:
    return f"{int(round(float(value))):,}"


def fp(value: str | float) -> str:
    p = float(value)
    if p < 0.001:
        return "<0.001"
    return f"{p:.3f}"


def fcell(value: object, digits: int) -> str:
    if isinstance(value, str):
        return value
    if pd.isna(value):
        return ""
    if isinstance(value, (int,)) or (isinstance(value, float) and float(value).is_integer() and digits == 0):
        return f"{int(value):,}"
    return f"{float(value):,.{digits}f}"


def markdown_table(frame: pd.DataFrame, digits: int | dict[str, int] = 2) -> str:
    """Render a table with the column labels exactly as supplied."""
    body = []
    for _, row in frame.iterrows():
        cells = []
        for column in frame.columns:
            d = digits.get(column, 2) if isinstance(digits, dict) else digits
            value = row[column]
            if pd.api.types.is_integer_dtype(frame[column]):
                cells.append(f"{int(value):,}")
            else:
                cells.append(fcell(value, d))
        body.append(cells)
    headers = [str(x) for x in frame.columns]
    # Separator lengths set pandoc's relative column widths, so wide content
    # gets proportionally wide columns and the table spans the text width.
    widths = [
        4 + max(3, min(40, max([max(len(w) for w in h.split())] + [len(row[i]) for row in body])))
        for i, h in enumerate(headers)
    ]
    lines = [
        "| " + " | ".join(headers) + " |",
        "| " + " | ".join("-" * w for w in widths) + " |",
    ]
    lines += ["| " + " | ".join(cells) + " |" for cells in body]
    return "\n".join(lines)


def num(value: float, digits: int = 2) -> str:
    """Format a number with a typographic minus sign."""
    return f"{value:.{digits}f}".replace("-", "−")


def ci(row: pd.Series, prefix: str, digits: int = 2) -> str:
    """Table cell: risk difference in percentage points with 95% CI."""
    return (
        f"{num(100 * row[prefix], digits)} ({num(100 * row[prefix + '_lower'], digits)} "
        f"to {num(100 * row[prefix + '_upper'], digits)})"
    )


def cit(row: pd.Series, prefix: str, unit: str = "percentage points") -> str:
    """Running text: risk difference with 95% CI, without nested parentheses."""
    return (
        f"{num(100 * row[prefix])} {unit} (95% CI, {num(100 * row[prefix + '_lower'])} "
        f"to {num(100 * row[prefix + '_upper'])})"
    )


def cit_in(row: pd.Series, prefix: str) -> str:
    """Like cit(), for use inside parentheses."""
    return cit(row, prefix).replace(" (95% CI, ", "; 95% CI, ")[:-1]


def drop_in(row: pd.Series, prefix: str) -> str:
    return drop(row, prefix).replace(" (95% CI, ", "; 95% CI, ")[:-1]


def drop(row: pd.Series, prefix: str) -> str:
    """Running text for a decrease, expressed as a positive magnitude."""
    return (
        f"{-100 * row[prefix]:.2f} percentage points (95% CI, {-100 * row[prefix + '_upper']:.2f} "
        f"to {-100 * row[prefix + '_lower']:.2f})"
    )


def includes_zero(row: pd.Series, prefix: str) -> bool:
    return row[prefix + "_lower"] <= 0 <= row[prefix + "_upper"]


def selected_annual(outcome: str, model: str) -> pd.DataFrame:
    frame = pd.read_csv(TABLES / "annual_adjusted_mortality.csv")
    return frame[(frame["outcome"] == outcome) & (frame["model"] == model)].copy()


def selected_contrast(outcome: str, model: str) -> pd.Series:
    frame = pd.read_csv(TABLES / "model_contrasts.csv")
    row = frame[(frame["outcome"] == outcome) & (frame["model"] == model)]
    if len(row) != 1:
        raise RuntimeError(f"Expected one contrast for {outcome}/{model}; got {len(row)}")
    return row.iloc[0]


def pstr(value: str | float) -> str:
    p = float(value)
    return "*P* < 0.001" if p < 0.001 else f"*P* = {p:.3f}"


def build_documents() -> None:
    MANUSCRIPT.mkdir(parents=True, exist_ok=True)
    SUBMISSION.mkdir(parents=True, exist_ok=True)
    k = key_values()
    model = k["selected_model"]
    analysis_mode = k["analysis_mode"]
    enriched = analysis_mode == "EHR-enriched"
    adjustment = "prespecified baseline case mix" if enriched else "age and sex"
    internal_status = (
        "EHR enrichment detected; verify multiple-imputation and validation outputs before submission."
        if enriched
        else "INTERNAL WORKING DRAFT: EHR enrichment, clinician text validation, and confirmation of author contributions remain required before submission."
    )
    audit = json.loads((AUDIT / "data_quality_summary.json").read_text(encoding="utf-8"))

    one = selected_contrast("1-year all-cause mortality", model)
    early = selected_contrast("30-day all-cause mortality", model)
    long = selected_contrast("5-year all-cause mortality", model)
    case_mix_one = selected_contrast("1-year all-cause mortality", "available_case_mix")
    case_mix_early = selected_contrast("30-day all-cause mortality", "available_case_mix")

    # Sensitivity analyses: the first row per outcome is the selected model.
    sens = pd.read_csv(TABLES / "sensitivity_analyses.csv")
    sens_one = sens[sens["outcome"] == "1-year all-cause mortality"].reset_index(drop=True)
    sens_early = sens[sens["outcome"] == "30-day all-cause mortality"].reset_index(drop=True)
    s1 = sens_one.iloc[0]
    s30 = sens_early.iloc[0]
    alt_one = sens_one.iloc[1:]
    alt_early = sens_early.iloc[1:]
    n_alt = len(alt_one)
    one_net_null = int(sum(includes_zero(r, "rd_end_start") for _, r in alt_one.iterrows()))
    one_rise_robust = int(sum(r["rd_peak_start_lower"] > 0 for _, r in alt_one.iterrows()))
    one_fall_robust = int(sum(r["rd_end_peak_upper"] < 0 for _, r in alt_one.iterrows()))
    early_fall_robust = int(sum(r["rd_end_peak_upper"] < 0 for _, r in alt_early.iterrows()))
    peak_years_one = sorted(set(int(y) for y in sens_one["peak_year"]))
    lookback_one = sens_one[sens_one["analysis"].str.startswith("Fixed 2-year lookback")].iloc[0]
    lookback_early = sens_early[sens_early["analysis"].str.startswith("Fixed 2-year lookback")].iloc[0]
    source_one = sens_one[sens_one["analysis"].str.startswith("Adjusted for report source")].iloc[0]
    source_early = sens_early[sens_early["analysis"].str.startswith("Adjusted for report source")].iloc[0]
    categorical_one = sens_one[sens_one["analysis"] == "Year as categorical"].iloc[0]

    bridge = pd.read_csv(TABLES / "source_bridge_2022.csv")
    b1 = bridge[bridge["outcome"] == "1-year all-cause mortality"].set_index("source")
    b30 = bridge[bridge["outcome"] == "30-day all-cause mortality"].set_index("source")
    period = pd.read_csv(TABLES / "period_mortality_2015_2018_cohorts.csv").set_index("calendar_year")
    cox = pd.read_csv(RESULTS / "cox_summary.csv").set_index("outcome")
    cox1 = cox.loc["1-year all-cause mortality"]
    cox30 = cox.loc["30-day all-cause mortality"]
    cox5 = cox.loc["5-year all-cause mortality"]
    crude = pd.read_csv(TABLES / "annual_crude_mortality.csv")
    crude_one = crude[crude["outcome"] == "1-year all-cause mortality"]
    crude_min = crude_one.loc[crude_one["crude_risk"].idxmin()]
    crude_max = crude_one.loc[crude_one["crude_risk"].idxmax()]
    case_mix_annual = pd.read_csv(TABLES / "annual_case_mix.csv").set_index("index_year")

    def rr(year: int) -> str:
        row = period.loc[year]
        return (
            f"{row['adjusted_rate_ratio_vs_2019']:.2f}; "
            f"95% CI, {row['rr_lower']:.2f} to {row['rr_upper']:.2f}"
        )

    def case_mix_text(row: pd.Series) -> str:
        return (
            f"{num(100 * row['risk_difference'])} percentage points (95% CI, "
            f"{num(100 * row['rd_lower'])} to {num(100 * row['rd_upper'])})"
        )

    both_robust = one_rise_robust == n_alt and one_fall_robust == n_alt
    robust_phrase = (
        f"in all {n_alt} sensitivity analyses"
        if both_robust
        else f"in {one_rise_robust} and {one_fall_robust} of {n_alt} sensitivity analyses, respectively"
    )

    # Table 1.
    table1_full = pd.read_csv(TABLES / "table1_baseline_by_period.csv")
    table1_labels = {
        "period": "Index period",
        "n": "Patients, n",
        "age": "Age, y",
        "female": "Women",
        "primary_pci": "Primary PCI (text)",
        "stemi_primary": "STEMI or primary PCI (text)",
        "left_main": "Left main (text)",
        "graft_pci": "Graft PCI (text)",
        "cto": "CTO (text)",
        "restenosis": "Restenosis (text)",
        "des": "DES (text)",
        "bms": "BMS (text)",
        "one_year_death": "1-year death",
    }
    table1_full = table1_full.rename(columns=table1_labels)
    table1 = table1_full[
        ["Index period", "Patients, n", "Age, y", "Women", "Primary PCI (text)", "DES (text)", "1-year death"]
    ]
    table1_md = markdown_table(table1)
    missing_sex = int(audit.get("missing_sex_patients", 0))
    missing_age = int(audit.get("missing_age_patients", 0))

    # Table 2: fixed-horizon estimates with the descriptive peak contrasts.
    rows = []
    for label, cohort, n_key, e_key, srow, last_year in [
        ("30-day", "2015–2025", "early_n", "early_events", s30, 2025),
        ("1-year", "2015–2024", "primary_n", "primary_events", s1, 2024),
    ]:
        rows.append(
            {
                "Outcome (cohort)": f"{label} ({cohort})",
                "N": int(float(k[n_key])),
                "Deaths": int(float(k[e_key])),
                "2015, %": f"{100 * srow['start_risk']:.2f}",
                "Peak, % (year)": f"{100 * srow['peak_risk']:.2f} ({int(srow['peak_year'])})",
                f"Final year, %": f"{100 * srow['end_risk']:.2f} ({last_year})",
                "Peak − 2015": ci(srow, "rd_peak_start"),
                "Final − peak": ci(srow, "rd_end_peak"),
                "Final − 2015": ci(srow, "rd_end_start"),
                "*P*": fp(srow["p_year"]),
            }
        )
    rows.append(
        {
            "Outcome (cohort)": "5-year (2015–2020)",
            "N": int(float(k["long_n"])),
            "Deaths": int(float(k["long_events"])),
            "2015, %": f"{100 * long['start_risk']:.2f}",
            "Peak, % (year)": "—",
            "Final year, %": f"{100 * long['end_risk']:.2f} (2020)",
            "Peak − 2015": "—",
            "Final − peak": "—",
            "Final − 2015": f"{num(100 * long['risk_difference'])} ({num(100 * long['rd_lower'])} to {num(100 * long['rd_upper'])})",
            "*P*": fp(long["p_year"]),
        }
    )
    table2_md = markdown_table(pd.DataFrame(rows))

    # Table 3: sensitivity analyses for the primary outcome.
    def short_label(label: str) -> str:
        if label.startswith("Primary"):
            return "Primary model"
        if label.startswith("Age/sex + text-derived"):
            return "Plus text-derived case mix"
        if label.startswith("Index years"):
            return "Excluding 2015"
        if label.startswith("Adjusted for report source"):
            return "Report-source adjusted"
        if label.startswith("Fixed 2-year lookback"):
            return "2-year look-back episodes"
        if label == "Year as categorical":
            return "Categorical year"
        return label

    t3 = pd.DataFrame(
        {
            "Analysis": sens_one["analysis"].map(short_label),
            "N": sens_one["model_n"].astype(int),
            "Deaths": sens_one["events"].astype(int),
            "First year, %": [f"{100 * x:.2f} ({int(y)})" for x, y in zip(sens_one["start_risk"], sens_one["start_year"])],
            "Peak, % (year)": [f"{100 * x:.2f} ({int(y)})" for x, y in zip(sens_one["peak_risk"], sens_one["peak_year"])],
            "2024, %": [f"{100 * x:.2f}" for x in sens_one["end_risk"]],
            "Peak − first": [ci(r, "rd_peak_start") for _, r in sens_one.iterrows()],
            "2024 − peak": [ci(r, "rd_end_peak") for _, r in sens_one.iterrows()],
            "2024 − first": [ci(r, "rd_end_start") for _, r in sens_one.iterrows()],
            "*P*": [fp(x) for x in sens_one["p_year"]],
        }
    )
    table3_md = markdown_table(t3)

    adjustment_methods = (
        "The selected model adjusted for age, sex, clinical presentation, cardiogenic shock, cardiac arrest, diabetes, hypertension, smoking, heart failure, previous cardiovascular disease, dialysis, eGFR, hemoglobin, and LVEF when the prespecified completeness threshold was met."
        if enriched
        else "Because the EHR enrichment extract was not available for this working analysis, the selected preliminary model adjusted for age (3-df natural spline) and sex. A separate exploratory model additionally used unvalidated report-text classifications for clinical presentation, left-main involvement, and graft PCI; it was not treated as the primary estimate."
    )

    paper = f"""---
title: "Temporal Trends in One-Year All-Cause Mortality After Percutaneous Coronary Intervention: A Decade-Long Cohort Linked to Civil Registry Data"
author:
  - "Fatih Köksal, MD"
  - "Fatih Levent, MD"
  - "Fatih Koca, MD"
  - "Kübra Severgün, MD"
  - "Tolga Doğan, MD"
  - "Mahmut Kapsız, MD"
  - "Mehmet Melek, MD"
  - "Fahriye Vatansever Ağca, MD"
  - "Erhan Tenekecioğlu, MD"
  - "Hasan Arı, MD"
date: "2026-10-01"
bibliography: references.bib
link-citations: true
format:
  docx: default
  html: default
---

> **{internal_status}**

# Abstract

## Background

Temporal trends in mortality after percutaneous coronary intervention (PCI) differ across registries and may have been disrupted by the COVID-19 pandemic. We examined annual all-cause mortality after PCI in a Turkish cohort linked to civil-registry data.

## Methods

This retrospective cohort included each adult patient's first observed PCI during 2015–2024; the primary outcome was 1-year all-cause death. Secondary cohorts assessed 30-day (2015–2025) and 5-year (2015–2020) mortality. Annual risks were standardized for {adjustment} from spline models of calendar year. Sensitivity analyses addressed the year specification, a 2022 change in reporting system, and unequal look-back.

## Results

Among {fint(k['primary_n'])} patients (median age, {fnum(k['median_age'], 0)} years; {fnum(k['female_percent'])}% women), {fint(k['primary_events'])} ({fnum(k['primary_crude_percent'])}%) died within 1 year. Adjusted 1-year mortality rose from {100 * s1['start_risk']:.2f}% in 2015 to {100 * s1['peak_risk']:.2f}% in {int(s1['peak_year'])} (difference, {cit_in(s1, 'rd_peak_start')}) and fell to {100 * s1['end_risk']:.2f}% in 2024 (decrease, {drop_in(s1, 'rd_end_peak')}); both changes persisted {robust_phrase}. The net 2024-versus-2015 difference ({cit_in(s1, 'rd_end_start')}) was not robust, with confidence intervals including zero in {one_net_null} of {n_alt} sensitivity analyses. Thirty-day mortality rose similarly, but its decline was attenuated after adjustment for report source. Five-year mortality was not associated with calendar year ({pstr(long['p_year'])}). Among 1-year survivors of 2015–2018 procedures, mortality was higher in calendar year 2020 than in 2019 (rate ratio, {rr(2020)}).

## Conclusions

Mortality after PCI did not improve steadily: it peaked around 2020–2022 and then declined, and the net change across the decade is uncertain. Pandemic-era background mortality probably contributed to the peak. Full clinical case-mix adjustment is required.

**Keywords:** percutaneous coronary intervention; mortality; temporal trends; COVID-19; registry; routinely collected health data

# Introduction

Percutaneous coronary intervention has evolved through changes in stent technology, vascular access, intracoronary imaging, antithrombotic treatment, and the treatment of increasingly complex and high-risk patients [@lawton2022]. Risk-adjusted outcome surveillance is therefore an important component of quality assessment in coronary revascularization programs [@dehmer2023; @castrodominguez2021]. Crude mortality trends may mislead when the age, comorbidity, acuity, and anatomy of treated patients change over time.

Prior evidence has not shown a uniform temporal pattern. In Western Denmark, 1-year mortality after primary PCI for ST-segment elevation myocardial infarction declined steadily from 2003 through 2018 [@thrane2023]. By contrast, all-comer registries from Australia, Canada, Korea, the United States, and England and Wales reported stable or rising crude and predicted mortality as treated patients became older, more acutely ill, and more comorbid [@dawson2021; @tran2019; @choi2022; @alkhouli2020; @kataruka2020; @ayayo2024]. Calendar period was frequently not independently associated with mortality after adjustment [@dawson2021], and demographic change explained part, but not all, of temporal survival differences [@hulme2019]. A systematic review of all-comers trials found declining cardiac death but not a consistent decline in all-cause mortality [@asano2022], consistent with a shift from cardiac to noncardiac causes of death after PCI [@spoon2014].

The COVID-19 pandemic added a further disruption. International registries documented fewer primary PCI procedures, longer ischemic times, and higher short-term mortality in 2020 [@deluca2022; @chew2021; @rodriguezleor2020], and Türkiye experienced substantial excess all-cause mortality during 2020–2022 [@keskin2025; @wang2022]. National Turkish data have characterized outcomes after acute myocardial infarction [@kilickap2021], but long, contemporary series spanning all PCI presentations, linked to civil-registry mortality, and covering the pre-pandemic, pandemic, and post-pandemic periods remain limited. We therefore examined temporal changes in 1-year all-cause mortality after first observed PCI from 2015 through 2024, together with 30-day mortality through 2025 and 5-year mortality in mature cohorts.

# Methods

## Study design and data sources

This retrospective observational cohort used routinely collected PCI reports from the catheterization laboratory of the University of Health Sciences, Bursa Yüksek İhtisas Training and Research Hospital, Bursa, Türkiye. Reports were linked at the institution to MERNİS civil-registry vital status. The source workbook contained structured fields and free-text procedure documentation generated from PDF reports until August 2022 and from the PARS system thereafter. The raw source was retained unchanged; patient identity was used only for in-memory linkage, after which a study identifier replaced all direct identifiers.

The study protocol was approved by the Clinical Research Ethics Committee of the University of Health Sciences, Bursa Yüksek İhtisas Training and Research Hospital (decision no. 2024-TBEK 2026/05-22; 20 May 2026). Because of the retrospective design and use of de-identified data, the requirement for informed consent was waived. The study was conducted in accordance with the Declaration of Helsinki. Reporting follows STROBE and its RECORD extension for routinely collected health data [@vonelm2007; @benchimol2015].

## Study population

Reports were consolidated by patient and calendar date, so that same-day reports formed one procedure day. The primary analysis used each patient's earliest observed PCI during 1 January 2015 through 31 December 2024. Because the source contained almost no reports before 2015, "first observed" does not imply first lifetime PCI. Patients whose recorded death date preceded the index PCI were considered to have irreconcilable linkage errors and were excluded from all analyses. Implausible ages (<18 or >110 years) were set to missing and therefore omitted from adjusted complete-case models. The 2025 cohort was used only for the 30-day analysis; 2026 was excluded as a partial calendar year.

## Variables and outcomes

Age, sex, procedure date, report source, procedure description, clinical diagnosis, procedural note, and result text were obtained from the PCI reports. Preliminary rule-based classifications identified primary PCI, clinical presentation, left-main or graft involvement, chronic total occlusion, restenosis, and stent type. These classifications were used descriptively and in one exploratory model; they require locked clinician validation before journal submission.

The outcome was all-cause death recorded in MERNİS, which is not subject to cause-of-death misclassification [@gaudino2020]. For decedents, follow-up ended on the death date. Patients without a recorded death were censored at their latest individually verified vital-status query date. Dates were normalized to calendar days; same-day deaths were assigned 0.5 days in time-to-event models. The primary outcome was death within 365 days, and secondary outcomes were death within 30 and 1,825 days. A patient contributed to a fixed-horizon analysis only if vital status was known through the horizon or death occurred within it.

## Statistical analysis

Continuous variables are presented as median and interquartile range and categorical variables as counts and percentages. The primary fixed-horizon model was logistic regression with index year represented by a 3-degree-of-freedom natural spline. Annual mortality was marginally standardized to the covariate distribution of the pooled analysis cohort [@muller2014], with delta-method confidence intervals. The global contribution of calendar year was tested by likelihood-ratio comparison with an otherwise identical model omitting year. {adjustment_methods}

The prespecified estimand was the standardized risk in each calendar year and the absolute risk difference between the last and first years. Because the modeled curve was non-monotonic, we also report, as descriptive post hoc contrasts, the difference between the peak year and 2015 and between the final year and the peak year. Analogous models evaluated 30-day mortality in 2015–2025 and 5-year mortality in 2015–2020. Cox regression censored at each horizon, with the same covariates, was a secondary analysis; proportional hazards were assessed with Schoenfeld residuals. Five-year Kaplan–Meier survival was plotted for 2015–2017 versus 2018–2020.

## Sensitivity analyses

Sensitivity analyses, all added after review of the preliminary results, examined threats specific to this data source. First, calendar year was modeled with 2-, 4-, and 5-df splines and as a categorical variable. Second, the first index year (2015), in which capture and text structure differed from later years, was excluded. Third, report source (PDF vs PARS) was added to the model and risks were standardized as if all procedures had been documented in PDF reports; the source effect is identified by the abrupt switch in September 2022 relative to the smooth calendar trend. In addition, outcomes were compared between the two systems within the 2022 bridge year. Fourth, selecting each patient's first observed procedure produces a look-back that lengthens with calendar year: a 2015 patient with an earlier, unrecorded PCI is included, whereas a 2024 patient with a recorded 2018 PCI is not. This is analogous to prevalent-user bias and left truncation [@danaei2012; @applebaum2011]. We therefore constructed a fixed-look-back cohort of all procedure episodes during 2017–2025 (2017–2024 for 1-year mortality) that had no recorded PCI in the preceding 730 days, so every episode had the same observable 2-year look-back; patients could contribute more than one episode, and variance was estimated with a patient-clustered sandwich estimator and a robust Wald test.

Finally, to distinguish procedure-cohort effects from calendar-period (background) mortality, we followed patients treated during 2015–2018 who survived 1 year and estimated their death rate in each calendar year from 2016 through 2025 using Poisson regression with person-time offsets, adjusted for attained age and sex. Time since PCI was not included because it is nearly collinear with calendar year in this closed cohort.

Analyses were performed in R version 4.3.3. All tests were two-sided, and 95% confidence intervals are reported. No causal interpretation was assigned to calendar-year associations.

# Results

## Cohort formation and temporal case mix

The source comprised {fint(k['source_records'])} PCI reports, consolidated into {fint(audit['procedure_days'])} procedure days from {fint(k['unique_patients'])} patients. After selection of the first observed PCI, {fint(k['outside_primary_period'])} patients had an index procedure outside 2015–2024 (mostly in 2025–2026), and {fint(k['invalid_death_before_index_primary'])} patients in the primary period had a recorded death date preceding the index procedure and were excluded, leaving {fint(k['primary_n'])} patients (Figure 1). Across all index years, {fint(k['invalid_death_before_index'])} patients with this inconsistency were excluded from every analysis. The adjusted model included {fint(k['primary_model_n'])} patients with complete age and sex.

Median age was {fnum(k['median_age'], 0)} years and {fnum(k['female_percent'])}% were women. The annual cohort size increased from {fint(case_mix_annual.loc[2015, 'n'])} in 2015 to {fint(case_mix_annual.loc[2024, 'n'])} in 2024, with lower volumes in 2018 ({fint(case_mix_annual.loc[2018, 'n'])}) and 2020 ({fint(case_mix_annual.loc[2020, 'n'])}). Median age increased from 60 to 62 years. Text classification suggested increasing DES and decreasing BMS use, but these trends span the change in reporting system and are not interpreted as validated clinical estimates (Table 1; Figure 2).

## One-year mortality

Overall, {fint(k['primary_events'])} patients died within 1 year ({fnum(k['primary_crude_percent'])}%). Crude annual risk ranged from {100 * crude_min['crude_risk']:.1f}% in {int(crude_min['index_year'])} to {100 * crude_max['crude_risk']:.1f}% in {int(crude_max['index_year'])}. After adjustment for {adjustment}, standardized risk was nonlinearly associated with index year ({pstr(s1['p_year'])}). Modeled risk was {100 * s1['start_risk']:.2f}% in 2015, peaked at {100 * s1['peak_risk']:.2f}% in {int(s1['peak_year'])}, and declined to {100 * s1['end_risk']:.2f}% in 2024 (Figure 3; Table 2). The peak exceeded the 2015 risk by {cit(s1, 'rd_peak_start')}, and the 2024 risk was lower than the peak by {drop(s1, 'rd_end_peak')}. The prespecified 2024-versus-2015 difference was {cit(s1, 'rd_end_start')}. With additional adjustment for text-derived presentation, left-main, and graft PCI, the corresponding difference was {case_mix_text(case_mix_one)}. In the secondary Cox analysis, calendar year was associated with 1-year mortality ({pstr(cox1['p_year_lrt'])}) without evidence that its effect varied over follow-up (proportional-hazards {pstr(cox1['ph_p_year'])}).

## Thirty-day and five-year mortality

The 30-day cohort included {fint(k['early_n'])} patients and {fint(k['early_events'])} deaths ({fnum(k['early_crude_percent'])}%). Adjusted risk increased from {100 * s30['start_risk']:.2f}% in 2015 to a peak of {100 * s30['peak_risk']:.2f}% in {int(s30['peak_year'])} and declined to {100 * s30['end_risk']:.2f}% in 2025 ({pstr(s30['p_year'])}; Table 2). The 2025-versus-2015 difference was {cit(s30, 'rd_end_start')}; with text-derived case-mix adjustment it was {case_mix_text(case_mix_early)}.

Among {fint(k['long_n'])} patients with a mature 5-year horizon, {fint(k['long_events'])} deaths occurred ({fnum(k['long_crude_percent'])}%). Five-year adjusted mortality was {100 * long['start_risk']:.2f}% for 2015 and {100 * long['end_risk']:.2f}% for 2020; calendar year was not significantly associated with 5-year mortality ({pstr(long['p_year'])}; Figure 4). In the corresponding Cox model, the calendar-year association was likewise not significant ({pstr(cox5['p_year_lrt'])}), but its hazard ratio changed over follow-up (proportional-hazards {pstr(cox5['ph_p_year'])}), consistent with a calendar-period rather than an index-year effect.

## Sensitivity analyses

The rise from 2015 to the peak was significant in {one_rise_robust} of {n_alt} sensitivity analyses of 1-year mortality, and the decline from the peak to 2024 in {one_fall_robust} of {n_alt}. The peak year ranged from {peak_years_one[0]} to {peak_years_one[-1]} (Table 3; Supplementary Figure S1). In contrast, the 2024-versus-2015 difference varied from {num(100 * alt_one['rd_end_start'].min())} to {num(100 * alt_one['rd_end_start'].max())} percentage points, and its confidence interval included zero in {one_net_null} of {n_alt} analyses. With year as a categorical variable, the 2024 risk was {100 * categorical_one['end_risk']:.2f}% versus {100 * categorical_one['start_risk']:.2f}% in 2015.

In the fixed 2-year look-back cohort ({fint(lookback_one['model_n'])} episodes in {fint(lookback_one['patients'])} patients), 1-year risk rose from {100 * lookback_one['start_risk']:.2f}% in 2017 to {100 * lookback_one['peak_risk']:.2f}% in {int(lookback_one['peak_year'])} and fell to {100 * lookback_one['end_risk']:.2f}% in 2024 (robust Wald {pstr(lookback_one['p_year'])}). Thirty-day results were similar (robust Wald {pstr(lookback_early['p_year'])}).

The reporting system changed in a single step: PDF reports covered January–August 2022 and PARS reports covered August–December 2022. Within 2022, 1-year mortality was {b1.loc['PDF', 'crude_percent']:.1f}% with PDF and {b1.loc['PARS', 'crude_percent']:.1f}% with PARS documentation (adjusted odds ratio for PARS, {b1.loc['PARS', 'adjusted_or_pars']:.2f}; 95% CI, {b1.loc['PARS', 'or_lower']:.2f} to {b1.loc['PARS', 'or_upper']:.2f}). For 30-day mortality, the values were {b30.loc['PDF', 'crude_percent']:.1f}% and {b30.loc['PARS', 'crude_percent']:.1f}% (odds ratio, {b30.loc['PARS', 'adjusted_or_pars']:.2f}; 95% CI, {b30.loc['PARS', 'or_lower']:.2f} to {b30.loc['PARS', 'or_upper']:.2f}). After adjustment for report source, the 1-year curve was essentially unchanged (2024 versus peak, {cit_in(source_one, 'rd_end_peak')}). The post-peak decline in 30-day mortality, however, was attenuated and no longer significant ({cit_in(source_early, 'rd_end_peak')}).

## Calendar-period mortality among earlier cohorts

Among 1-year survivors of PCI performed during 2015–2018, the attained-age- and sex-adjusted death rate was higher in calendar year 2020 than in 2019 (rate ratio, {rr(2020)}). The rate ratio was {rr(2021)} in 2021 and returned to values near 1 from 2022 onward (Supplementary Table S7).

# Discussion

## Principal findings

In this decade-long cohort linked to civil-registry mortality, four findings emerged. First, 1-year mortality after PCI did not improve steadily: adjusted risk rose from 2015 to a peak around 2020–2022 and declined thereafter. Second, this rise-and-fall shape was robust to the calendar-year specification, exclusion of the first year, adjustment for the change in reporting system, and a fixed-look-back cohort definition. By contrast, the prespecified net comparison of 2024 with 2015 was sensitive to modeling choices and should not be interpreted as evidence of a decade-long deterioration. Third, 30-day mortality showed a similar rise, but its post-peak decline was not robust to adjustment for the reporting system. Fourth, long-term survivors of earlier procedures experienced excess mortality in calendar year 2020, suggesting that part of the peak reflects period mortality during the pandemic rather than changes in procedural care.

## Comparison with previous studies

The absence of a monotonic decline is consistent with all-comer registries in which increasing age, acuity, and comorbidity offset procedural progress [@dawson2021; @tran2019; @choi2022; @alkhouli2020; @kataruka2020; @ayayo2024]. It contrasts with the steady improvement after primary PCI in Western Denmark [@thrane2023], where the study period ended before the pandemic and the population was restricted to a single presentation. Relative-survival analyses from British registry data showed that demographic change explains part of temporal survival differences [@hulme2019]. In addition, the proportion of noncardiac deaths after PCI has increased over time [@spoon2014]; all-cause mortality is therefore sensitive to background mortality, which is exactly what changed during 2020–2022.

The timing of the peak coincides with the COVID-19 pandemic. Multinational registries reported fewer primary PCI procedures, longer total ischemic times, and higher short-term mortality in 2020 [@deluca2022; @chew2021; @rodriguezleor2020]. In our cohort, the annual volume fell in 2020, and the share of text-classified primary PCI increased, consistent with selective deferral of elective procedures. Türkiye experienced an estimated 247,640 excess deaths during 2020–2022, with a peak in 2021 [@keskin2025; see also @wang2022]. The excess calendar-period mortality among stable, long-term survivors of earlier PCI supports a contribution of background pandemic mortality, separate from any change in the quality of PCI. This remains an interpretation: without cause-of-death and infection data, we cannot separate COVID-19 deaths from indirect effects.

## Methodological implications

Two features of routinely collected PCI data warrant attention in temporal analyses. First, selecting each patient's first observed procedure from a data source that begins at a fixed date creates a look-back that lengthens over calendar time, so early years contain more patients with previous, unrecorded revascularization. This is analogous to prevalent-user bias [@danaei2012; @applebaum2011]. A fixed-look-back episode definition is a simple remedy, and in our data it did not change the conclusions. Second, information-system transitions can create artificial discontinuities. Here the transition was abrupt and outcome rates were similar on either side of it, but the 30-day decline after the transition was sensitive to source adjustment and should be confirmed with validated clinical variables.

## Clinical and quality implications

Annual PCI outcome review should combine absolute mortality with transparent case-mix measurement and should be interpreted alongside background population mortality. A single early-versus-late comparison would have been misleading in this cohort: depending on modeling choices, it suggested either a modest increase or no change, and it hid the larger rise and fall in between. Linking catheterization-laboratory data to complete vital status and maintaining stable definitions across information-system transitions provide a practical foundation for local quality surveillance [@dehmer2023; @castrodominguez2021].

## Limitations

This study has several limitations. It was observational and single-center, leaving residual confounding and limited generalizability. The present working analysis lacks the planned EHR comorbidities, renal function, hemoglobin, LVEF, shock, cardiac arrest, and complete angiographic complexity; its adjusted findings are therefore preliminary. Text-derived classifications have not yet undergone the prespecified blinded clinical validation. Only all-cause death was available, so cardiac and noncardiac deaths and COVID-19 infection could not be distinguished. Vital status was assessed on individual query dates rather than one common administrative date, although all fixed-horizon cohorts were required to have mature follow-up. The sensitivity analyses, the peak contrasts, and the calendar-period analysis were not prespecified and should be regarded as exploratory. The calendar-period analysis was confined to 2015–2018 cohorts and cannot fully separate age, period, and cohort effects. Finally, calendar-year associations cannot identify the contribution of individual devices, medications, the pandemic, or changes in referral patterns.

# Conclusions

All-cause mortality after first observed PCI varied nonlinearly during 2015–2024. One-year risk rose to a peak around 2020–2022 and declined thereafter, a pattern robust to model specification, reporting-system change, and cohort definition; 30-day risk rose similarly, but its subsequent decline was less certain. The net change across the decade was uncertain, and five-year mortality showed no significant temporal association in mature cohorts. Excess calendar-period mortality during the pandemic probably contributed to the peak. Complete EHR case-mix enrichment and source validation are required before these findings can support a submission-ready assessment of temporal PCI outcomes.

```{{=openxml}}
<w:p><w:r><w:br w:type="page"/></w:r></w:p>
```

# Tables

## Table 1. Characteristics of the primary cohort by index period

{table1_md}

Values are median (interquartile range) or n (%). Percentages use patients with non-missing data as the denominator (sex missing in {missing_sex} and age missing or implausible in {missing_age} patients across all years). Variables marked "text" were derived from report text by preliminary rules, are unvalidated, and span the September 2022 change from PDF to PARS reports. The complete descriptive table is provided in Supplementary Table S2. DES indicates drug-eluting stent; PCI, percutaneous coronary intervention.

```{{=openxml}}
<w:p><w:pPr><w:sectPr><w:pgSz w:w="12240" w:h="15840"/><w:pgMar w:top="1440" w:right="1440" w:bottom="1440" w:left="1440" w:header="720" w:footer="720" w:gutter="0"/></w:sectPr></w:pPr></w:p>
```

## Table 2. Fixed-horizon mortality: standardized risks and calendar-year contrasts

{table2_md}

Risks are standardized to the {adjustment} distribution of each pooled cohort. Differences are in percentage points (95% CI). The 5-year model did not have an interior maximum, so peak contrasts are not shown. *P* is the likelihood-ratio test for all calendar-year spline terms. Peak contrasts are descriptive and post hoc.

```{{=openxml}}
<w:p><w:r><w:br w:type="page"/></w:r></w:p>
```

## Table 3. Sensitivity analyses of 1-year mortality

{table3_md}

Risks are standardized; differences are in percentage points (95% CI). Unless stated, models adjust for {adjustment} with a 3-df spline for calendar year. "Plus text-derived case mix" adds unvalidated presentation, left-main, and graft PCI classifications. The source-adjusted model standardizes all procedures to PDF documentation. The look-back analysis includes every procedure episode without a recorded PCI in the preceding 730 days (N counts episodes, from {fint(lookback_one['patients'])} patients); patients may contribute more than one episode, and confidence intervals and the Wald test use patient-clustered robust variance. *P* is the global test for calendar year. Corresponding 30-day results are in Supplementary Table S5.

```{{=openxml}}
<w:p><w:pPr><w:sectPr><w:pgSz w:w="15840" w:h="12240" w:orient="landscape"/><w:pgMar w:top="1080" w:right="1080" w:bottom="1080" w:left="1080" w:header="720" w:footer="720" w:gutter="0"/></w:sectPr></w:pPr></w:p>
```

# Figure legends

**Figure 1. Study cohort flow.** Consolidation of reports into procedure days, selection of each patient's first observed PCI, exclusions, and formation of the primary and secondary fixed-horizon cohorts.

![Study cohort flow](../outputs/figures/figure1_cohort_flow.png){{width=85%}}

**Figure 2. Temporal changes in cohort size and selected characteristics.** Text-classified primary PCI is preliminary and spans the change in reporting system.

![Case-mix trends](../outputs/figures/figure2_case_mix_trends.png){{width=95%}}

**Figure 3. Crude and adjusted 30-day and 1-year all-cause mortality.** Lines show annual standardized estimates; shaded areas are 95% confidence intervals.

![Adjusted mortality trends](../outputs/figures/figure3_adjusted_mortality_trends.png){{width=95%}}

**Figure 4. Five-year all-cause survival in mature cohorts.** Kaplan–Meier estimates compare index periods 2015–2017 and 2018–2020.

![Five-year survival](../outputs/figures/figure4_five_year_survival.png){{width=85%}}

# End matter

**Ethics approval:** Approved by the Clinical Research Ethics Committee of the University of Health Sciences, Bursa Yüksek İhtisas Training and Research Hospital (decision no. 2024-TBEK 2026/05-22; 20 May 2026). The requirement for informed consent was waived because of the retrospective design.

**Funding:** This research received no specific grant from any funding agency in the public, commercial, or not-for-profit sectors.

**Conflicts of interest:** The authors declare no conflicts of interest.

**Author contributions (CRediT):** Fatih Köksal: Conceptualization, Methodology, Formal analysis, Investigation, Writing – original draft, Writing – review & editing. Fatih Levent, Fatih Koca, Kübra Severgün, Tolga Doğan, and Mahmut Kapsız: Data curation, Writing – review & editing. Mehmet Melek, Fahriye Vatansever Ağca, and Erhan Tenekecioğlu: Supervision, Writing – review & editing. Hasan Arı: Supervision, Project administration, Writing – review & editing.

**Data availability:** Patient-level data contain protected health information and cannot be shared publicly. De-identified aggregate outputs, the complete analysis code, and the variable dictionary will be made publicly available at **[repository DOI/URL required]**.

**Declaration of generative AI and AI-assisted technologies:** During manuscript preparation, the authors used AI coding assistants (OpenAI Codex and Anthropic Claude Code) to assist with reproducible code generation, sensitivity-analysis implementation, literature retrieval, language drafting, and consistency checks. The authors reviewed and verified all analyses, citations, and text and take full responsibility for the content.

# References
"""

    # Supplement.
    def pct(x: float) -> str:
        return f"{100 * x:.1f}"

    maturity = pd.read_csv(AUDIT / "followup_maturity_by_year.csv")
    maturity.columns = [
        "Index year", "Patients, n", "Known at 30 d", "Deaths by 30 d",
        "Known at 1 y", "Deaths by 1 y", "Known at 5 y", "Deaths by 5 y",
    ]

    crude_s = crude.copy()
    crude_s["outcome"] = crude_s["outcome"].str.replace(" all-cause mortality", "", regex=False)
    crude_s["crude_risk"] = crude_s["crude_risk"].map(pct)
    crude_s.columns = ["Outcome", "Index year", "Patients, n", "Deaths, n", "Crude mortality, %"]

    contrasts = pd.read_csv(TABLES / "model_contrasts.csv")
    model_labels = {
        "unadjusted": "Unadjusted",
        "age_sex_adjusted": "Age/sex",
        "available_case_mix": "Age/sex + text-derived case mix",
    }
    contrasts_s = pd.DataFrame(
        {
            "Outcome": contrasts["outcome"].str.replace(" all-cause mortality", "", regex=False),
            "Model": contrasts["model"].map(model_labels),
            "Patients, n": contrasts["model_n"].astype(int),
            "First year, %": [f"{pct(r)} ({int(y)})" for r, y in zip(contrasts["start_risk"], contrasts["start_year"])],
            "Final year, %": [f"{pct(r)} ({int(y)})" for r, y in zip(contrasts["end_risk"], contrasts["end_year"])],
            "Difference, points (95% CI)": [
                f"{100 * a:.2f} ({100 * b:.2f} to {100 * c:.2f})"
                for a, b, c in zip(contrasts["risk_difference"], contrasts["rd_lower"], contrasts["rd_upper"])
            ],
            "Global *P*": contrasts["p_year"].map(fp),
        }
    )

    t5 = t3.copy()
    t5_early = pd.DataFrame(
        {
            "Analysis": sens_early["analysis"].map(short_label),
            "N": sens_early["model_n"].astype(int),
            "Deaths": sens_early["events"].astype(int),
            "First year, %": [f"{100 * x:.2f} ({int(y)})" for x, y in zip(sens_early["start_risk"], sens_early["start_year"])],
            "Peak, % (year)": [f"{100 * x:.2f} ({int(y)})" for x, y in zip(sens_early["peak_risk"], sens_early["peak_year"])],
            "2025, %": [f"{100 * x:.2f}" for x in sens_early["end_risk"]],
            "Peak − first": [ci(r, "rd_peak_start") for _, r in sens_early.iterrows()],
            "2025 − peak": [ci(r, "rd_end_peak") for _, r in sens_early.iterrows()],
            "2025 − first": [ci(r, "rd_end_start") for _, r in sens_early.iterrows()],
            "*P*": [fp(x) for x in sens_early["p_year"]],
        }
    )

    bridge_s = bridge.copy()
    bridge_s = pd.DataFrame(
        {
            "Outcome": bridge_s["outcome"].str.replace(" all-cause mortality", "", regex=False),
            "Report source": bridge_s["source"],
            "Months": bridge_s["first_month"] + " to " + bridge_s["last_month"],
            "Patients, n": bridge_s["n"].astype(int),
            "Median age, y": bridge_s["median_age"].map(lambda x: f"{x:.0f}"),
            "Women, %": bridge_s["female_percent"].map(lambda x: f"{x:.1f}"),
            "STEMI or primary PCI (text), %": bridge_s["stemi_primary_percent"].map(lambda x: f"{x:.1f}"),
            "Deaths, n": bridge_s["deaths"].astype(int),
            "Crude mortality, %": bridge_s["crude_percent"].map(lambda x: f"{x:.1f}"),
            "Adjusted OR, PARS vs PDF (95% CI)": [
                f"{o:.2f} ({lo:.2f} to {hi:.2f})" if src == "PARS" else "Reference"
                for o, lo, hi, src in zip(bridge_s["adjusted_or_pars"], bridge_s["or_lower"], bridge_s["or_upper"], bridge_s["source"])
            ],
        }
    )

    period_s = period.reset_index()
    period_s = pd.DataFrame(
        {
            "Calendar year": period_s["calendar_year"].astype(int),
            "Person-years": period_s["person_years"].map(lambda x: f"{x:,.0f}"),
            "Deaths, n": period_s["deaths"].astype(int),
            "Crude rate per 100 person-years": period_s["crude_rate_per_100py"].map(lambda x: f"{x:.2f}"),
            "Adjusted rate ratio (95% CI)": [
                "1.00 (reference)" if y == 2019 else f"{r:.2f} ({lo:.2f} to {hi:.2f})"
                for y, r, lo, hi in zip(period_s["calendar_year"], period_s["adjusted_rate_ratio_vs_2019"], period_s["rr_lower"], period_s["rr_upper"])
            ],
        }
    )
    period_p = fp(period["p_calendar_year"].iloc[0])

    cox_s = cox.reset_index()
    cox_s = pd.DataFrame(
        {
            "Outcome": cox_s["outcome"].str.replace(" all-cause mortality", "", regex=False),
            "Patients, n": cox_s["n"].astype(int),
            "Deaths, n": cox_s["events"].astype(int),
            "Calendar year, LRT *P*": cox_s["p_year_lrt"].map(fp),
            "PH test for calendar year, *P*": cox_s["ph_p_year"].map(fp),
            "Global PH test, *P*": cox_s["ph_p_global"].map(fp),
        }
    )

    annual_cm = pd.read_csv(TABLES / "annual_case_mix.csv")
    annual_cm = pd.DataFrame(
        {
            "Index year": annual_cm["index_year"].astype(int),
            "Patients, n": annual_cm["n"].astype(int),
            "Median age, y": annual_cm["median_age"].map(lambda x: f"{x:.0f}"),
            "Women, %": annual_cm["female_percent"].map(lambda x: f"{x:.1f}"),
            "Primary PCI (text), %": annual_cm["primary_pci_percent"].map(lambda x: f"{x:.1f}"),
            "Left main (text), %": annual_cm["left_main_percent"].map(lambda x: f"{x:.1f}"),
            "Graft PCI (text), %": annual_cm["graft_pci_percent"].map(lambda x: f"{x:.1f}"),
            "DES (text), %": annual_cm["des_percent"].map(lambda x: f"{x:.1f}"),
        }
    )

    supplement = f"""# Supplementary appendix

> **Status:** {internal_status}

## Supplementary Methods

The raw workbook was hashed before analysis. Direct identifiers were used only in memory to consolidate repeated reports. Same-day reports were concatenated before selection of the first observed PCI. Death and index datetimes were normalized to calendar dates so that deaths on the procedure date were not misclassified as negative follow-up. The protected patient-level analysis files contain no name, national identifier, catheter number, or operator name.

The primary fixed-horizon regression and standardization are described in the main manuscript and `docs/SAP.md`. This run was classified as **{analysis_mode}**, and the model selected for manuscript estimates was `{model}`. The sensitivity analyses in Supplementary Tables S5–S8 and Supplementary Figure S1 were added after review of the preliminary results; they are documented as Amendment 1 to the statistical analysis plan.

## Supplementary Table S1. Follow-up maturity by index year

Patients are counted as "known" at a horizon if vital status was verified beyond the horizon or death occurred within it.

{markdown_table(maturity, digits=0)}

## Supplementary Table S2. Detailed primary-cohort characteristics by period

{markdown_table(table1_full)}

Values are median (interquartile range) or n (%). "STEMI or primary PCI (text)" includes every report classified as primary PCI plus reports with STEMI terms in other fields, so it is necessarily at least as large as "Primary PCI (text)".

## Supplementary Table S3. Annual crude mortality

{markdown_table(crude_s)}

## Supplementary Table S4. Standardized first- and final-year risks by adjustment tier

{markdown_table(contrasts_s)}

## Supplementary Table S5. Sensitivity analyses of 30-day mortality

{markdown_table(t5_early)}

## Supplementary Table S6. Outcomes by report source in the 2022 bridge year

{markdown_table(bridge_s)}

Odds ratios are adjusted for {adjustment}.

## Supplementary Table S7. Calendar-period mortality among 1-year survivors of PCI performed in 2015–2018

{markdown_table(period_s)}

Poisson regression with a person-time offset, adjusted for attained age (3-df natural spline) and sex; reference year 2019. Global test for calendar year {pstr(period["p_calendar_year"].iloc[0])}.

## Supplementary Table S8. Secondary Cox models

{markdown_table(cox_s)}

Cox models censored at each horizon with the same covariates as the selected logistic model. LRT indicates likelihood-ratio test; PH, proportional hazards (Schoenfeld residuals).

## Supplementary Table S9. Annual case-mix indicators

{markdown_table(annual_cm)}

## Supplementary Figure S1. Calendar-year pattern across sensitivity analyses

![Sensitivity analyses](../outputs/figures/figureS1_sensitivity_trends.png){{width=95%}}

Standardized 30-day and 1-year mortality from the primary model, the model with categorical year, the model standardized to PDF documentation, and the fixed 2-year look-back episode cohort.

## Outstanding pre-submission analyses

1. Merge the approved EHR enrichment extract and run 40-fold multiple imputation; repeat all sensitivity analyses with the enriched model.
2. Complete blinded validation of text-derived presentation and anatomy variables, stratified by report source.
3. Run the admission-level repeated-PCI sensitivity analysis after encounter IDs are supplied.
4. Obtain cause-of-death and COVID-19 infection data, if available, to separate direct and indirect pandemic mortality.
5. Add the missingness table and imputation diagnostics.
"""

    cover_letter = f"""# Cover letter — working draft

**To:** Editor, Catheterization & Cardiovascular Interventions  
**Article type:** Original Article — Clinical Science

Dear Editor,

Please consider our manuscript, “Temporal Trends in One-Year All-Cause Mortality
After Percutaneous Coronary Intervention: A Decade-Long Cohort Linked to Civil
Registry Data.”

We examined {fint(k['primary_n'])} patients undergoing a first observed PCI from
2015 through 2024, with all-cause mortality ascertained through civil-registry
records. Adjusted 1-year mortality did not improve steadily: it rose to a peak
around 2020–2022 and declined thereafter. This pattern was robust to the
calendar-year specification, a change in the hospital reporting system, and a
fixed-look-back cohort definition. Long-term survivors of earlier procedures had
excess mortality in calendar year 2020, which points to a contribution of
pandemic-era background mortality. A mature 30-day cohort included
{fint(k['early_n'])} patients, and {fint(k['long_n'])} patients contributed to
the 5-year analysis.

The manuscript is relevant to interventional clinicians and registry leads
because it shows how a simple early-versus-late comparison can mislead. It also
shows how cohort-selection and information-system artifacts can be tested in
routinely collected PCI data. All analyses are reproducible, and the code and
data dictionary will be shared while protected patient-level data remain secure.

This manuscript is original, is not under consideration elsewhere, and has been
approved by all authors. **[Replace this sentence only after all authors confirm.]**
Ethics approval was provided by the Clinical Research Ethics Committee of the
University of Health Sciences, Bursa Yüksek İhtisas Training and Research
Hospital (decision no. 2024-TBEK 2026/05-22; 20 May 2026). Funding and
conflict-of-interest statements are included in the manuscript.

Sincerely,

Fatih Köksal, MD  
Department of Cardiology, University of Health Sciences, Bursa Yüksek İhtisas
Training and Research Hospital, Mimarsinan Mah., Emniyet Cad., 16310 Yıldırım,
Bursa, Türkiye  
Email: dr.fatihkoksal@hotmail.com · ORCID: 0000-0002-4197-4683
"""

    abstract_words = len(
        re.findall(
            r"\b[\w%–.-]+\b",
            re.sub(r"[#*]", "", paper.split("# Abstract", 1)[1].split("**Keywords:**", 1)[0]),
        )
    )
    main_words = len(
        re.findall(
            r"\b[\w%–.-]+\b",
            re.sub(r"\[@[^\]]+\]", "", re.sub(r"[#*]", "", paper.split("# Introduction", 1)[1].split("# Tables", 1)[0])),
        )
    )

    title_page = f"""# Title page

## Full title

Temporal Trends in One-Year All-Cause Mortality After Percutaneous Coronary
Intervention: A Decade-Long Cohort Linked to Civil Registry Data

## Running title

Temporal Mortality Trends After PCI

## Authors and affiliations

Fatih Köksal, MD^a^ · Fatih Levent, MD^a^ · Fatih Koca, MD^a^ · Kübra Severgün, MD^a^ · Tolga Doğan, MD^a^ · Mahmut Kapsız, MD^a^ · Mehmet Melek, MD^a^ · Fahriye Vatansever Ağca, MD^a^ · Erhan Tenekecioğlu, MD^a^ · Hasan Arı, MD^a^

^a^ Department of Cardiology, University of Health Sciences, Bursa Yüksek İhtisas Training and Research Hospital, Mimarsinan Mah., Emniyet Cad., 16310 Yıldırım, Bursa, Türkiye

| Author | ORCID |
| ---------------- | ------------------- |
| Fatih Köksal | 0000-0002-4197-4683 |
| Fatih Levent | 0000-0002-7160-4050 |
| Fatih Koca | 0000-0002-6824-1017 |
| Kübra Severgün | 0000-0001-5789-5691 |
| Tolga Doğan | 0000-0001-7074-1419 |
| Mahmut Kapsız | 0000-0001-5412-594X |
| Mehmet Melek | 0000-0002-9510-9206 |
| Fahriye Vatansever Ağca | 0000-0002-1401-2980 |
| Erhan Tenekecioğlu | 0000-0003-4376-2833 |
| Hasan Arı | 0000-0002-9681-2374 |

## Corresponding author

Fatih Köksal, MD  
Department of Cardiology, University of Health Sciences, Bursa Yüksek İhtisas Training and Research Hospital, Mimarsinan Mah., Emniyet Cad., 16310 Yıldırım, Bursa, Türkiye  
Email: dr.fatihkoksal@hotmail.com  
ORCID: 0000-0002-4197-4683  
Telephone: [omitted from public repository]

## Word count and contents

Abstract: {abstract_words} words  
Main text (Introduction through Conclusions, excluding citations): approximately {main_words} words  
Tables: 3  
Figures: 4  
Supplementary appendix: 1 (9 tables, 1 figure)

## Funding, disclosures, and acknowledgments

Funding: This research received no specific grant from any funding agency in the public, commercial, or not-for-profit sectors.  
Conflicts of interest: The authors declare no conflicts of interest.  
Ethics: Clinical Research Ethics Committee of the University of Health Sciences, Bursa Yüksek İhtisas Training and Research Hospital, decision no. 2024-TBEK 2026/05-22 (20 May 2026); informed consent waived.
"""

    (MANUSCRIPT / "paper.qmd").write_text(paper, encoding="utf-8")
    (MANUSCRIPT / "supplement.md").write_text(supplement, encoding="utf-8")
    (SUBMISSION / "cover_letter.md").write_text(cover_letter, encoding="utf-8")
    (SUBMISSION / "title_page.md").write_text(title_page, encoding="utf-8")
    print(f"Wrote manuscript in {analysis_mode} mode")
    print(f"Abstract words: {abstract_words}; main text words: {main_words}")


if __name__ == "__main__":
    build_documents()
