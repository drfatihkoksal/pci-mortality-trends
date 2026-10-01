#!/usr/bin/env python3
"""Build a de-identified, first-observed PCI cohort from the source workbook.

Direct identifiers are held only in memory. The output uses a sequential study
identifier and deliberately omits names, national identifiers, catheter numbers,
and physician names.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import re
from pathlib import Path

import numpy as np
import pandas as pd


REQUIRED_COLUMNS = {
    "patientID",
    "RaporTarihi",
    "yas_anjiyo",
    "cinsiyet",
    "Kaynak",
    "Islem",
    "KlinikTani",
    "IslemNotu",
    "Sonuc",
    "MernisDurum",
    "VefatTarihi",
    "TakipSuresi_Gun_Sansurlu",
}

DIRECT_IDENTIFIERS = {
    "TcKimlik",
    "patientID",
    "HastaAdSoyad",
    "KataterNo",
    "DoktorAdi1",
    "DoktorAdi2",
}

OPTIONAL_EHR_COLUMNS = [
    "diabetes",
    "hypertension",
    "current_smoker",
    "heart_failure",
    "prior_mi",
    "prior_pci",
    "prior_cabg",
    "prior_stroke",
    "peripheral_artery_disease",
    "dialysis",
    "egfr",
    "hemoglobin",
    "lvef",
    "cardiogenic_shock",
    "cardiac_arrest",
    "ehr_presentation",
    "radial_access",
    "multivessel_disease",
    "procedural_success",
    "intravascular_imaging",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--input",
        default="KoronerArterMudahale_birlesik_mortalite.xlsx",
        type=Path,
    )
    parser.add_argument(
        "--output", default="data/private/analysis_cohort.csv", type=Path
    )
    parser.add_argument(
        "--audit",
        default="outputs/audit/data_quality_summary.json",
        type=Path,
    )
    parser.add_argument(
        "--ehr", default="data/private/ehr_enrichment.csv", type=Path
    )
    return parser.parse_args()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def normalize_patient_id(series: pd.Series) -> pd.Series:
    out = series.astype("string").str.strip()
    return out.replace({"": pd.NA, "nan": pd.NA, "None": pd.NA})


def text_join(series: pd.Series) -> str:
    values = [str(v).strip() for v in series.dropna() if str(v).strip()]
    return " | ".join(dict.fromkeys(values))


def any_true(series: pd.Series) -> bool:
    return bool(series.fillna(False).astype(bool).any())


def classify_text(frame: pd.DataFrame) -> pd.DataFrame:
    text = (
        frame["procedure_text"].fillna("")
        + " "
        + frame["diagnosis_text"].fillna("")
        + " "
        + frame["note_text"].fillna("")
        + " "
        + frame["result_text"].fillna("")
    ).str.casefold()

    def contains(pattern: str) -> pd.Series:
        return text.str.contains(pattern, regex=True, na=False)

    primary = frame["procedure_text"].fillna("").str.casefold().str.contains(
        r"primer|primary", regex=True, na=False
    )
    stemi = primary | contains(
        r"\bstemi\b|akut\s+(?:anterior|inferior|posterior|lateral).*\bmi\b|"
        r"\ba\.\s*(?:anterior|inferior|posterior|lateral)\s*mi\b"
    )
    nstemi = contains(r"\bnstemi\b|non[- ]?st")
    other_acs = contains(r"\baks\b|akut koroner|unstable|kararsız")
    presentation = np.select(
        [stemi, nstemi, other_acs],
        ["STEMI/primary PCI", "NSTEMI", "Other ACS"],
        default="Chronic/uncertain",
    )

    frame = frame.copy()
    frame["primary_pci"] = primary.astype(int)
    frame["clinical_presentation"] = presentation
    frame["left_main"] = contains(
        r"sol ana koroner|left main|\blmca\b|(?<![a-z])lm(?![a-z])"
    ).astype(int)
    frame["graft_pci"] = contains(r"safen|\blima\b|greft|graft|cabg").astype(int)
    frame["cto"] = contains(r"kronik total|\bcto\b|total okl").astype(int)
    frame["restenosis"] = contains(r"restenoz|in[- ]?stent|instent").astype(int)
    frame["des"] = contains(
        r"ilaç salınımlı|drug[- ]?eluting|\bdes\b|sirolimus|everolimus|"
        r"zotarolimus|biolimus"
    ).astype(int)
    frame["bms"] = contains(r"çıplak stent|bare[- ]?metal|\bbms\b").astype(int)
    frame["target_lad"] = contains(r"(?<![a-z])lad(?![a-z])").astype(int)
    frame["target_rca"] = contains(r"(?<![a-z])rca(?![a-z])").astype(int)
    frame["target_cx"] = contains(r"(?<![a-z])cx(?![a-z])|sirkumfleks").astype(int)
    return frame


def build_cohort(df: pd.DataFrame) -> tuple[pd.DataFrame, dict]:
    missing = sorted(REQUIRED_COLUMNS.difference(df.columns))
    if missing:
        raise ValueError(f"Missing required columns: {', '.join(missing)}")

    audit: dict[str, object] = {
        "source_rows": int(len(df)),
        "source_columns": int(len(df.columns)),
        "exact_duplicate_rows": int(df.duplicated().sum()),
        "direct_identifier_columns_present": sorted(
            DIRECT_IDENTIFIERS.intersection(df.columns)
        ),
    }

    df = df.copy()
    df["_pid"] = normalize_patient_id(df["patientID"])
    df["_report_datetime"] = pd.to_datetime(
        df["RaporTarihi"], errors="coerce", dayfirst=True
    )
    df["_index_date"] = df["_report_datetime"].dt.normalize()
    df["_death_date"] = pd.to_datetime(
        df["VefatTarihi"], errors="coerce", dayfirst=True
    ).dt.normalize()
    df["_dead"] = (
        df["MernisDurum"]
        .astype("string")
        .str.strip()
        .str.casefold()
        .eq("ölüm")
    )
    df["_followup_days"] = pd.to_numeric(
        df["TakipSuresi_Gun_Sansurlu"], errors="coerce"
    )
    df["_inferred_status_date"] = df["_index_date"] + pd.to_timedelta(
        df["_followup_days"], unit="D"
    )
    df["_age"] = pd.to_numeric(df["yas_anjiyo"], errors="coerce")
    invalid_age = df["_age"].notna() & ~df["_age"].between(18, 110)
    audit["invalid_age_rows"] = int(invalid_age.sum())
    df.loc[invalid_age, "_age"] = np.nan

    sex = df["cinsiyet"].astype("string").str.strip().str.upper()
    df["_sex"] = sex.map({"E": "Male", "M": "Male", "K": "Female", "F": "Female"})

    audit["missing_patient_id_rows"] = int(df["_pid"].isna().sum())
    audit["missing_index_date_rows"] = int(df["_index_date"].isna().sum())
    audit["unique_source_patients"] = int(df["_pid"].nunique(dropna=True))

    valid = df.dropna(subset=["_pid", "_index_date"]).copy()
    patient_order = pd.Series(valid["_pid"].drop_duplicates().tolist())
    study_map = {pid: f"P{i:06d}" for i, pid in enumerate(patient_order, start=1)}

    # Consolidate multiple reports on the same calendar day before choosing the
    # first observed PCI. Text fields are combined, not selected opportunistically.
    day_level = (
        valid.sort_values(["_pid", "_report_datetime"])
        .groupby(["_pid", "_index_date"], as_index=False)
        .agg(
            age=("_age", "first"),
            sex=("_sex", "first"),
            source=("Kaynak", text_join),
            report_type=("RaporTuru", text_join),
            procedure_text=("Islem", text_join),
            diagnosis_text=("KlinikTani", text_join),
            note_text=("IslemNotu", text_join),
            result_text=("Sonuc", text_join),
            dead=("_dead", any_true),
            death_date=("_death_date", "first"),
            inferred_status_date=("_inferred_status_date", "max"),
            reports_same_day=("_pid", "size"),
        )
    )

    patient_status = (
        valid.groupby("_pid", as_index=False)
        .agg(
            ever_dead=("_dead", any_true),
            patient_death_date=("_death_date", "first"),
            last_alive_status_date=("_inferred_status_date", "max"),
            source_record_count=("_pid", "size"),
        )
    )
    index = (
        day_level.sort_values(["_pid", "_index_date"])
        .drop_duplicates("_pid", keep="first")
        .merge(patient_status, on="_pid", how="left", validate="one_to_one")
    )
    index["study_id"] = index["_pid"].map(study_map)
    index["dead"] = index["ever_dead"].astype(int)
    index["death_date"] = index["patient_death_date"]
    index["last_vital_status_date"] = np.where(
        index["dead"].eq(1),
        index["death_date"],
        index["last_alive_status_date"],
    )
    index["last_vital_status_date"] = pd.to_datetime(
        index["last_vital_status_date"], errors="coerce"
    )
    index["followup_days"] = (
        index["last_vital_status_date"] - index["_index_date"]
    ).dt.days
    index["death_before_index"] = (
        index["dead"].eq(1) & index["followup_days"].lt(0)
    ).astype(int)
    index["same_day_death"] = (
        index["dead"].eq(1) & index["followup_days"].eq(0)
    ).astype(int)
    index["index_year"] = index["_index_date"].dt.year

    index = classify_text(index)

    # Every consolidated procedure day, used for the fixed-lookback sensitivity
    # analysis that does not depend on the start of the source data.
    days = day_level.merge(patient_status, on="_pid", how="left", validate="many_to_one")
    days = days.sort_values(["_pid", "_index_date"]).reset_index(drop=True)
    days["study_id"] = days["_pid"].map(study_map)
    days["dead"] = days["ever_dead"].astype(int)
    days["death_date"] = days["patient_death_date"]
    days["last_vital_status_date"] = pd.to_datetime(
        np.where(days["dead"].eq(1), days["death_date"], days["last_alive_status_date"]),
        errors="coerce",
    )
    days["followup_days"] = (days["last_vital_status_date"] - days["_index_date"]).dt.days
    days["death_before_index"] = (days["dead"].eq(1) & days["followup_days"].lt(0)).astype(int)
    days["index_year"] = days["_index_date"].dt.year
    days["days_since_prior_procedure"] = (
        days.groupby("_pid")["_index_date"].diff().dt.days
    )
    days = classify_text(days).rename(columns={"_index_date": "index_date"})
    days = days[
        [
            "study_id",
            "index_date",
            "index_year",
            "age",
            "sex",
            "source",
            "dead",
            "death_date",
            "last_vital_status_date",
            "followup_days",
            "death_before_index",
            "days_since_prior_procedure",
            "primary_pci",
            "clinical_presentation",
            "left_main",
            "graft_pci",
        ]
    ]

    output_columns = [
        "study_id",
        "_index_date",
        "index_year",
        "age",
        "sex",
        "source",
        "report_type",
        "dead",
        "death_date",
        "last_vital_status_date",
        "followup_days",
        "death_before_index",
        "same_day_death",
        "source_record_count",
        "reports_same_day",
        "primary_pci",
        "clinical_presentation",
        "left_main",
        "graft_pci",
        "cto",
        "restenosis",
        "des",
        "bms",
        "target_lad",
        "target_rca",
        "target_cx",
    ]
    index = index.rename(columns={"_index_date": "index_date"})
    output_columns[1] = "index_date"

    audit.update(
        {
            "same_day_report_rows_collapsed": int(
                (day_level["reports_same_day"] - 1).clip(lower=0).sum()
            ),
            "procedure_days": int(len(day_level)),
            "first_observed_patients": int(len(index)),
            "patients_with_multiple_source_records": int(
                index["source_record_count"].gt(1).sum()
            ),
            "death_before_index_patients": int(index["death_before_index"].sum()),
            "same_day_death_patients": int(index["same_day_death"].sum()),
            "index_date_min": str(index["index_date"].min().date()),
            "index_date_max": str(index["index_date"].max().date()),
            "first_observed_by_year": {
                str(year): int(count)
                for year, count in index["index_year"].value_counts().sort_index().items()
            },
            "missing_sex_patients": int(index["sex"].isna().sum()),
            "missing_age_patients": int(index["age"].isna().sum()),
            "source_system_by_year": {
                str(year): {
                    str(source): int(count)
                    for source, count in group["source"].value_counts().items()
                }
                for year, group in index.groupby("index_year")
            },
        }
    )

    return index[["_pid"] + output_columns], days, audit


def merge_ehr(index: pd.DataFrame, ehr_path: Path, audit: dict) -> pd.DataFrame:
    if not ehr_path.exists():
        audit["ehr_enrichment_status"] = "not supplied; limited-adjustment run"
        audit["ehr_columns_added"] = []
        return index

    ehr = pd.read_csv(ehr_path, dtype={"patientID": "string"})
    if "patientID" not in ehr.columns:
        raise ValueError("ehr_enrichment.csv must contain patientID")
    ehr["_pid"] = normalize_patient_id(ehr["patientID"])
    if ehr["_pid"].duplicated().any():
        raise ValueError("ehr_enrichment.csv must contain one row per patientID")
    available = [column for column in OPTIONAL_EHR_COLUMNS if column in ehr.columns]
    merged = index.merge(
        ehr[["_pid"] + available], on="_pid", how="left", validate="one_to_one"
    )
    audit["ehr_enrichment_status"] = "supplied"
    audit["ehr_columns_added"] = available
    audit["ehr_match_rate"] = round(
        float(merged[available].notna().any(axis=1).mean()) if available else 0.0, 4
    )
    return merged


def main() -> None:
    args = parse_args()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.audit.parent.mkdir(parents=True, exist_ok=True)

    source = pd.read_excel(args.input, sheet_name="Sheet1", engine="openpyxl")
    index, days, audit = build_cohort(source)
    audit["source_file"] = args.input.name
    audit["source_sha256"] = sha256(args.input)
    index = merge_ehr(index, args.ehr, audit)

    # Direct patient identifiers are removed before any persistent write.
    index = index.drop(columns=["_pid"])
    leaked = DIRECT_IDENTIFIERS.intersection(index.columns)
    if leaked:
        raise RuntimeError(f"Direct identifiers reached output: {sorted(leaked)}")

    index.to_csv(args.output, index=False, date_format="%Y-%m-%d")
    days_path = args.output.with_name("procedure_days.csv")
    if DIRECT_IDENTIFIERS.intersection(days.columns):
        raise RuntimeError("Direct identifiers reached procedure-day output")
    days.to_csv(days_path, index=False, date_format="%Y-%m-%d")
    with args.audit.open("w", encoding="utf-8") as handle:
        json.dump(audit, handle, ensure_ascii=False, indent=2)

    print(f"Wrote protected cohort: {args.output} ({len(index):,} patients)")
    print(f"Wrote protected procedure days: {days_path} ({len(days):,} days)")
    print(f"Wrote non-identifying audit: {args.audit}")
    print(f"EHR status: {audit['ehr_enrichment_status']}")


if __name__ == "__main__":
    main()
