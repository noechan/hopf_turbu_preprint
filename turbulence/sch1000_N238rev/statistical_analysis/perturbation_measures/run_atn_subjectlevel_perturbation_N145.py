#!/usr/bin/env python3
"""Relate subject-level Hopf perturbation measures to AT(N) biomarkers.

Information capacity and susceptibility are fitted in separate HC3 robust
ordinary least-squares models because they are almost perfectly correlated in
the current N145 export. Every model adjusts for age, sex, and years of
education from the canonical local ComBat metadata. BH-FDR is applied across
the five prespecified AT(N) outcomes separately within each predictor.

Only aggregate coefficients, diagnostics, prediction lines, and provenance
are written. The participant-level merged analysis table remains in memory.
"""

from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
from typing import Final

import numpy as np
import pandas as pd
import statsmodels.api as sm
from statsmodels.stats.multitest import multipletests


SCRIPT_DIR: Final = Path(__file__).resolve().parent
STATISTICAL_DIR: Final = SCRIPT_DIR.parent
SCH1000_ROOT: Final = STATISTICAL_DIR.parent
DEFAULT_ADNI3_ROOT: Final = Path(
    os.environ.get(
        "ADNI3_ROOT",
        "/path/to/ADNI3",
    )
)
MOUNTED_EXTERNAL_DIR: Final = (
    DEFAULT_ADNI3_ROOT
    / "code"
    / "HPC_Hopf_SUB_DTI_1000_Staging"
    / "visualization"
    / "data"
)
LOCAL_EXTERNAL_DIR: Final = (
    SCH1000_ROOT.parents[1]
    / "HPC_Hopf_SUB_DTI_1000_Staging"
    / "visualization"
    / "data"
)
DEFAULT_EXTERNAL_DIR: Final = (
    LOCAL_EXTERNAL_DIR if LOCAL_EXTERNAL_DIR.is_dir() else MOUNTED_EXTERNAL_DIR
)
LEGACY_HOPF_CSV: Final = (
    SCH1000_ROOT
    / "data_export"
    / "ML_Input_ADNI3_4STAGINGBYABETA_ComBat_N145_with_infocap_suscep_sch1000.csv"
)
LOCAL_HOPF_XLSX: Final = (
    SCH1000_ROOT.parents[1]
    / "machine_learning"
    / "Data"
    / "turbu_hopf"
    / "ML_Input_ADNI3_4STAGINGBYABETA_ComBat_N145_with_infocap_suscep_sch1000.xlsx"
)
DEFAULT_HOPF_FILE: Final = (
    LOCAL_HOPF_XLSX if LOCAL_HOPF_XLSX.is_file() else LEGACY_HOPF_CSV
)
DEFAULT_TAU_FILE: Final = (
    DEFAULT_EXTERNAL_DIR / "ADNI3_N238rev_with_ABETA_Status_CL24_tau_regional.xlsx"
)
DEFAULT_VBM_FILE: Final = DEFAULT_EXTERNAL_DIR / "ADNI3_VBM_postCOMBAT.csv"
LOCAL_METADATA_FILE: Final = (
    SCH1000_ROOT
    / "data"
    / "covariates"
    / "covariates_ADNI3_ABeta_N152.csv"
)
MOUNTED_METADATA_FILE: Final = (
    DEFAULT_ADNI3_ROOT
    / "code"
    / "ADNI3_neuroHarmonize_site"
    / "data"
    / "raw"
    / "turbu"
    / "covariates_ADNI3_ABeta_N152.csv"
)
DEFAULT_METADATA_FILE: Final = (
    LOCAL_METADATA_FILE if LOCAL_METADATA_FILE.is_file() else MOUNTED_METADATA_FILE
)
DEFAULT_OUTPUT_DIR: Final = (
    SCRIPT_DIR / "results" / "N145_ATN_subjectlevel_perturbation"
)

GROUP_ORDER: Final = ["HC_ABneg", "HC_ABpos", "MCI_ABpos", "AD_ABpos"]
EXPECTED_GROUP_COUNTS: Final = {
    "HC_ABneg": 51,
    "HC_ABpos": 37,
    "MCI_ABpos": 31,
    "AD_ABpos": 26,
}
PREDICTORS: Final = {
    "InformationCapacity": "Info_Cap",
    "Susceptibility": "Susceptibility",
}
OUTCOMES: Final = [
    {"panel": "Fig4c", "column": "CL_pvc", "label": "Amyloid burden (CL-PVC)"},
    {
        "panel": "Fig4d",
        "column": "Mean_GMV_HIP_BI",
        "label": "Bilateral hippocampal volume",
    },
    {"panel": "Fig4j", "column": "tau_mesial_pvc", "label": "Mesial tau (PVC)"},
    {
        "panel": "Fig4k",
        "column": "tau_metatemporal_pvc",
        "label": "Metatemporal tau (PVC)",
    },
    {
        "panel": "Fig4l",
        "column": "tau_temporoparietal_pvc",
        "label": "Temporoparietal tau (PVC)",
    },
]


def parse_args(cli_args: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--hopf-file", type=Path, default=DEFAULT_HOPF_FILE)
    parser.add_argument("--tau-file", type=Path, default=DEFAULT_TAU_FILE)
    parser.add_argument("--vbm-file", type=Path, default=DEFAULT_VBM_FILE)
    parser.add_argument("--metadata-file", type=Path, default=DEFAULT_METADATA_FILE)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument(
        "--validate-only",
        action="store_true",
        help="Validate inputs and samples without fitting or writing",
    )
    return parser.parse_args(cli_args)


def require_file(path: Path, label: str) -> None:
    if not path.is_file():
        raise FileNotFoundError(
            f"Missing {label}: {path}\nMount ADNI or pass the corresponding path option."
        )


def require_columns(frame: pd.DataFrame, columns: list[str], label: str) -> None:
    missing = [column for column in columns if column not in frame.columns]
    if missing:
        raise ValueError(f"{label} is missing columns: {missing}")


def normalize_ptid(frame: pd.DataFrame, label: str) -> pd.DataFrame:
    require_columns(frame, ["PTID"], label)
    result = frame.copy()
    result["PTID"] = result["PTID"].astype("string").str.strip()
    if result["PTID"].isna().any() or result["PTID"].eq("").any():
        raise ValueError(f"{label} contains an empty PTID")
    duplicates = result.loc[result["PTID"].duplicated(keep=False), "PTID"].unique()
    if len(duplicates):
        raise ValueError(f"{label} contains duplicate PTIDs: {duplicates.tolist()}")
    return result


def numeric(series: pd.Series, label: str) -> pd.Series:
    if pd.api.types.is_numeric_dtype(series):
        result = pd.to_numeric(series, errors="coerce")
    else:
        result = pd.to_numeric(
            series.astype("string").str.replace(",", ".", regex=False),
            errors="coerce",
        )
    unexpected = series.notna() & result.isna()
    if unexpected.any():
        values = sorted(series.loc[unexpected].astype(str).unique().tolist())
        raise ValueError(f"Could not convert {label} to numeric: {values[:10]}")
    return result.astype(float)


def read_hopf_table(path: Path) -> pd.DataFrame:
    """Read the canonical Hopf table from its distributed CSV or XLSX form."""

    suffix = path.suffix.lower()
    if suffix == ".csv":
        return pd.read_csv(path)
    if suffix in {".xlsx", ".xls"}:
        return pd.read_excel(path, engine="openpyxl")
    raise ValueError(f"Unsupported participant-level Hopf table format: {path}")


def portable_path(path: Path) -> str:
    resolved = path.resolve()
    for root, label in (
        (SCH1000_ROOT.resolve(), "$SCH1000_ROOT"),
        (DEFAULT_ADNI3_ROOT.resolve(), "$ADNI3_ROOT"),
    ):
        try:
            return str(Path(label) / resolved.relative_to(root))
        except ValueError:
            continue
    return str(resolved)


def load_analysis_data(args: argparse.Namespace) -> tuple[pd.DataFrame, float]:
    for path, label in (
        (args.hopf_file, "participant-level Hopf input"),
        (args.tau_file, "amyloid/regional-tau input"),
        (args.vbm_file, "post-ComBat VBM input"),
        (args.metadata_file, "canonical harmonisation metadata"),
    ):
        require_file(path, label)

    hopf = normalize_ptid(read_hopf_table(args.hopf_file), "Hopf input")
    tau = normalize_ptid(pd.read_excel(args.tau_file), "Amyloid/tau input")
    vbm = normalize_ptid(pd.read_csv(args.vbm_file), "VBM input")
    metadata = normalize_ptid(pd.read_csv(args.metadata_file), "Metadata input")

    tau_columns = [
        "PTID",
        "CL_pvc",
        "tau_mesial_pvc",
        "tau_metatemporal_pvc",
        "tau_temporoparietal_pvc",
    ]
    vbm_columns = ["PTID", "Mean_GMV_HIP_LH", "Mean_GMV_HIP_RH"]
    require_columns(hopf, ["PTID", "Group", *PREDICTORS.values()], "Hopf input")
    require_columns(tau, tau_columns, "Amyloid/tau input")
    require_columns(vbm, vbm_columns, "VBM input")
    require_columns(
        metadata,
        ["PTID", "Group", "age", "gender", "edu"],
        "Metadata input",
    )

    counts = hopf["Group"].value_counts().reindex(GROUP_ORDER, fill_value=0).to_dict()
    if counts != EXPECTED_GROUP_COUNTS:
        raise ValueError(
            "Hopf input is not the canonical N145 cohort. "
            f"Expected {EXPECTED_GROUP_COUNTS}; found {counts}."
        )

    merged = hopf[["PTID", "Group", *PREDICTORS.values()]].merge(
        tau[tau_columns], on="PTID", how="left", validate="one_to_one",
        indicator="_tau_match",
    )
    missing = merged.loc[merged["_tau_match"] != "both", "PTID"].tolist()
    if missing:
        raise ValueError(f"PTIDs missing from amyloid/tau input: {missing}")
    merged = merged.drop(columns="_tau_match")

    merged = merged.merge(
        vbm[vbm_columns], on="PTID", how="left", validate="one_to_one",
        indicator="_vbm_match",
    )
    missing = merged.loc[merged["_vbm_match"] != "both", "PTID"].tolist()
    if missing:
        raise ValueError(f"PTIDs missing from VBM input: {missing}")
    merged = merged.drop(columns="_vbm_match")

    merged = merged.merge(
        metadata[["PTID", "Group", "age", "gender", "edu"]],
        on="PTID", how="left", validate="one_to_one", suffixes=("", "_metadata"),
        indicator="_metadata_match",
    )
    missing = merged.loc[merged["_metadata_match"] != "both", "PTID"].tolist()
    if missing:
        raise ValueError(f"PTIDs missing from metadata input: {missing}")
    merged = merged.drop(columns="_metadata_match")
    mismatch = merged.loc[merged["Group"] != merged["Group_metadata"], "PTID"].tolist()
    if mismatch:
        raise ValueError(f"Group labels disagree between Hopf and metadata: {mismatch}")

    numeric_columns = [
        *PREDICTORS.values(), "age", "gender", "edu",
        "CL_pvc", "tau_mesial_pvc", "tau_metatemporal_pvc",
        "tau_temporoparietal_pvc", "Mean_GMV_HIP_LH", "Mean_GMV_HIP_RH",
    ]
    for column in numeric_columns:
        merged[column] = numeric(merged[column], column)
    if not set(merged["gender"].dropna().unique()).issubset({0.0, 1.0}):
        raise ValueError("Metadata gender must contain only 0 and 1")

    merged["AGE"] = merged["age"]
    merged["SEX_NUM"] = merged["gender"]
    merged["EDUCATION"] = merged["edu"]
    merged["Mean_GMV_HIP_BI"] = merged[
        ["Mean_GMV_HIP_LH", "Mean_GMV_HIP_RH"]
    ].mean(axis=1, skipna=False)

    required_complete = [
        *PREDICTORS.values(), "AGE", "SEX_NUM", "EDUCATION", "CL_pvc",
        "Mean_GMV_HIP_BI",
    ]
    missing_complete = merged[required_complete].isna().sum()
    if missing_complete.any():
        raise ValueError(
            "Unexpected missing values in complete N145 variables: "
            f"{missing_complete[missing_complete > 0].to_dict()}"
        )

    tau_outcomes = [
        "tau_mesial_pvc", "tau_metatemporal_pvc", "tau_temporoparietal_pvc",
    ]
    tau_complete = merged[tau_outcomes].notna()
    if not tau_complete.nunique(axis=1).eq(1).all():
        raise ValueError("The three regional tau outcomes have inconsistent missingness")
    tau_n = int(tau_complete.all(axis=1).sum())
    if tau_n != 134:
        raise ValueError(f"Expected 134 tau-complete participants; found {tau_n}")

    correlation = float(merged[list(PREDICTORS.values())].corr().iloc[0, 1])
    return merged, correlation


def fit_models(frame: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    covariates = ["AGE", "SEX_NUM", "EDUCATION"]
    results: list[dict[str, object]] = []
    predictions: list[pd.DataFrame] = []

    for predictor_label, predictor in PREDICTORS.items():
        for outcome in OUTCOMES:
            outcome_column = outcome["column"]
            model_frame = frame[[outcome_column, predictor, *covariates]].dropna().copy()
            predictor_mean = float(model_frame[predictor].mean())
            predictor_sd = float(model_frame[predictor].std(ddof=1))
            if not np.isfinite(predictor_sd) or predictor_sd <= 0:
                raise ValueError(f"{predictor} has zero or invalid SD for {outcome_column}")
            model_frame["PREDICTOR_Z"] = (
                model_frame[predictor] - predictor_mean
            ) / predictor_sd

            design = sm.add_constant(
                model_frame[["PREDICTOR_Z", *covariates]], has_constant="add"
            )
            model = sm.OLS(model_frame[outcome_column], design).fit(cov_type="HC3")
            reduced_design = sm.add_constant(model_frame[covariates], has_constant="add")
            reduced_model = sm.OLS(model_frame[outcome_column], reduced_design).fit()
            confidence = model.conf_int(alpha=0.05).loc["PREDICTOR_Z"]
            delta_r2 = float(model.rsquared - reduced_model.rsquared)
            partial_r2 = float(delta_r2 / (1.0 - reduced_model.rsquared))

            result = {
                "Predictor_Label": predictor_label,
                "Predictor": predictor,
                "Formula": f"{outcome_column} ~ z({predictor}) + Age + Sex + Education",
                "Panel": outcome["panel"],
                "Outcome": outcome_column,
                "Outcome_Label": outcome["label"],
                "N": int(model.nobs),
                "Residual_DF": float(model.df_resid),
                "Beta_Per_Predictor_SD": float(model.params["PREDICTOR_Z"]),
                "Beta_Per_Raw_Predictor_Unit": float(model.params["PREDICTOR_Z"] / predictor_sd),
                "HC3_SE_Per_Predictor_SD": float(model.bse["PREDICTOR_Z"]),
                "CI95_Lower_Per_Predictor_SD": float(confidence.iloc[0]),
                "CI95_Upper_Per_Predictor_SD": float(confidence.iloc[1]),
                "T": float(model.tvalues["PREDICTOR_Z"]),
                "P_Raw": float(model.pvalues["PREDICTOR_Z"]),
                "R2_Full_Model": float(model.rsquared),
                "Adjusted_R2_Full_Model": float(model.rsquared_adj),
                "R2_Covariates_Only": float(reduced_model.rsquared),
                "Delta_R2_Predictor": delta_r2,
                "Partial_R2_Predictor": partial_r2,
                "AIC": float(model.aic),
                "BIC": float(model.bic),
                "Predictor_Mean": predictor_mean,
                "Predictor_SD": predictor_sd,
                "Reference_AGE": float(model_frame["AGE"].mean()),
                "Reference_SEX_NUM": float(model_frame["SEX_NUM"].mode().iloc[0]),
                "Reference_EDUCATION": float(model_frame["EDUCATION"].mean()),
            }
            results.append(result)

            raw_values = np.linspace(
                float(model_frame[predictor].min()),
                float(model_frame[predictor].max()), 200,
            )
            prediction_design = pd.DataFrame(
                {
                    "const": np.ones(len(raw_values)),
                    "PREDICTOR_Z": (raw_values - predictor_mean) / predictor_sd,
                    "AGE": np.repeat(result["Reference_AGE"], len(raw_values)),
                    "SEX_NUM": np.repeat(result["Reference_SEX_NUM"], len(raw_values)),
                    "EDUCATION": np.repeat(result["Reference_EDUCATION"], len(raw_values)),
                }
            )
            predicted = model.predict(prediction_design[design.columns])
            predictions.append(
                pd.DataFrame(
                    {
                        "Predictor_Label": predictor_label,
                        "Predictor": predictor,
                        "Panel": outcome["panel"],
                        "Outcome": outcome_column,
                        "Predictor_Raw_Value": raw_values,
                        "Predicted_Outcome": predicted.to_numpy(),
                        "Reference_AGE": result["Reference_AGE"],
                        "Reference_SEX_NUM": result["Reference_SEX_NUM"],
                        "Reference_EDUCATION": result["Reference_EDUCATION"],
                    }
                )
            )

    results_frame = pd.DataFrame(results)
    results_frame["P_FDR_BH_Across_5_Outcomes"] = np.nan
    results_frame["Significant_FDR_0_05"] = False
    for predictor in PREDICTORS.values():
        mask = results_frame["Predictor"] == predictor
        rejected, adjusted, _, _ = multipletests(
            results_frame.loc[mask, "P_Raw"].to_numpy(), alpha=0.05, method="fdr_bh"
        )
        results_frame.loc[mask, "P_FDR_BH_Across_5_Outcomes"] = adjusted
        results_frame.loc[mask, "Significant_FDR_0_05"] = rejected
    return results_frame, pd.concat(predictions, ignore_index=True)


def main(cli_args: list[str] | None = None) -> None:
    args = parse_args(cli_args)
    frame, predictor_correlation = load_analysis_data(args)
    sample_sizes = {
        outcome["panel"]: int(
            frame[[outcome["column"], *PREDICTORS.values(), "AGE", "SEX_NUM", "EDUCATION"]]
            .dropna().shape[0]
        )
        for outcome in OUTCOMES
    }
    print("Validated canonical N145 cohort; tau-complete N=134.")
    print("Covariates: age, sex, and education from local ComBat metadata.")
    print(f"Information capacity--susceptibility correlation: r={predictor_correlation:.8f}")
    print(f"Outcome sample sizes: {sample_sizes}")
    if args.validate_only:
        print("Validation completed; no models were fitted and no files were written.")
        return

    results, predictions = fit_models(frame)
    output_dir = args.output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    result_path = output_dir / "N145_ATN_subjectlevel_perturbation_HC3_age_sex_education.csv"
    prediction_path = output_dir / "N145_ATN_subjectlevel_perturbation_prediction_lines_age_sex_education.csv"
    results.to_csv(result_path, index=False)
    predictions.to_csv(prediction_path, index=False)

    provenance = {
        "analysis": "N145 subject-level Hopf perturbation-measure--AT(N) associations",
        "hopf_input": portable_path(args.hopf_file),
        "amyloid_regional_tau_input": portable_path(args.tau_file),
        "postcombat_vbm_input": portable_path(args.vbm_file),
        "covariate_metadata_input": portable_path(args.metadata_file),
        "covariates": {"age": "age", "sex": "gender", "education": "edu"},
        "sex_coding": "female=0, male=1",
        "cohort_group_counts": EXPECTED_GROUP_COUNTS,
        "predictors": PREDICTORS,
        "predictor_correlation": predictor_correlation,
        "joint_predictor_model_fitted": False,
        "joint_model_reason": (
            "Information capacity and susceptibility are nearly collinear; "
            "separate models preserve interpretable total associations."
        ),
        "outcomes": OUTCOMES,
        "covariance": "HC3 heteroscedasticity-consistent",
        "multiple_testing": (
            "BH-FDR across five prespecified AT(N) outcomes separately within "
            "information capacity and susceptibility"
        ),
        "beta_definition": (
            "Outcome-unit change per one-SD increase in the tested Hopf measure; "
            "raw-unit coefficients are also saved"
        ),
        "group_role": "Cohort validation only; group is not a regression covariate",
        "participant_level_merged_data_saved": False,
        "aggregate_output_files": [str(result_path), str(prediction_path)],
    }
    provenance_path = output_dir / "N145_ATN_subjectlevel_perturbation_provenance.json"
    provenance_path.write_text(json.dumps(provenance, indent=2) + "\n", encoding="utf-8")

    print("\nEducation-adjusted AT(N) results:")
    for row in results.itertuples(index=False):
        print(
            f"  {row.Predictor_Label}, {row.Panel}: N={row.N}, "
            f"beta/SD={row.Beta_Per_Predictor_SD:.8g}, "
            f"p={row.P_Raw:.8g}, pFDR={row.P_FDR_BH_Across_5_Outcomes:.8g}, "
            f"deltaR2={row.Delta_R2_Predictor:.6g}"
        )
    print(f"\nSaved aggregate outputs and provenance under: {output_dir}")


if __name__ == "__main__":
    main()
