#!/usr/bin/env python3
"""Run the N145 turbulence--AT(N) regressions with two covariate sets.

The default ``both`` mode fits the prespecified age/sex model and the
age/sex/education model. Participant-level merged data remain in memory; only
aggregate coefficients, model diagnostics, prediction lines, and provenance
are written to the repository.
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
SCH1000_ROOT: Final = SCRIPT_DIR.parent.parent
DEFAULT_HARMONIZED_FILE: Final = (
    SCH1000_ROOT
    / "harmonization_allfeat"
    / "recomputed"
    / "Turbu_ComBat_ADNI3_allfeatures_N145.xlsx"
)
DEFAULT_ADNI3_ROOT: Final = Path(
    os.environ.get(
        "ADNI3_ROOT",
        "/path/to/ADNI3",
    )
)
DEFAULT_EXTERNAL_DIR: Final = (
    DEFAULT_ADNI3_ROOT
    / "code"
    / "HPC_Hopf_SUB_DTI_1000_Staging"
    / "visualization"
    / "data"
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
DEFAULT_OUTPUT_DIR: Final = SCRIPT_DIR / "results" / "N145_atn_regressions"

PREDICTOR: Final = "Turbu_lam_0_01"
MODEL_COVARIATES: Final = {
    "age_sex": ["AGE", "SEX_NUM"],
    "age_sex_education": ["AGE", "SEX_NUM", "EDUCATION"],
}
GROUP_ORDER: Final = ["HC_ABneg", "HC_ABpos", "MCI_ABpos", "AD_ABpos"]
EXPECTED_GROUP_COUNTS: Final = {
    "HC_ABneg": 51,
    "HC_ABpos": 37,
    "MCI_ABpos": 31,
    "AD_ABpos": 26,
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
    parser.add_argument(
        "--model",
        choices=("age_sex", "age_sex_education", "both"),
        default="both",
        help="Covariate specification to fit (default: both)",
    )
    parser.add_argument(
        "--harmonized-file",
        type=Path,
        default=DEFAULT_HARMONIZED_FILE,
        help="Current N145 all-feature ComBat workbook",
    )
    parser.add_argument(
        "--tau-file",
        type=Path,
        default=DEFAULT_TAU_FILE,
        help="Amyloid/demographic/regional-tau workbook",
    )
    parser.add_argument(
        "--vbm-file",
        type=Path,
        default=DEFAULT_VBM_FILE,
        help="Post-ComBat VBM table",
    )
    parser.add_argument(
        "--metadata-file",
        type=Path,
        default=DEFAULT_METADATA_FILE,
        help="ComBat metadata CSV containing PTID and edu",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=DEFAULT_OUTPUT_DIR,
        help="Destination for aggregate model outputs",
    )
    parser.add_argument(
        "--validate-only",
        action="store_true",
        help="Validate all inputs and model samples without fitting or writing files",
    )
    return parser.parse_args(cli_args)


def require_file(path: Path, label: str) -> None:
    if not path.is_file():
        raise FileNotFoundError(f"Missing {label}: {path}")


def portable_input_path(path: Path) -> str:
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


def to_numeric(series: pd.Series, label: str) -> pd.Series:
    if pd.api.types.is_numeric_dtype(series):
        return pd.to_numeric(series, errors="coerce")
    converted = pd.to_numeric(
        series.astype("string").str.replace(",", ".", regex=False),
        errors="coerce",
    )
    unexpected = series.notna() & converted.isna()
    if unexpected.any():
        values = sorted(series.loc[unexpected].astype(str).unique().tolist())
        raise ValueError(f"Could not convert {label} to numeric: {values[:10]}")
    return converted


def encode_sex(series: pd.Series) -> pd.Series:
    numeric = pd.to_numeric(series, errors="coerce")
    if numeric.notna().all() and set(numeric.unique()).issubset({0, 1}):
        return numeric.astype(float)
    normalized = series.astype("string").str.strip().str.lower()
    mapped = normalized.map({"f": 0.0, "female": 0.0, "m": 1.0, "male": 1.0})
    if mapped.isna().any():
        values = sorted(series.loc[mapped.isna()].astype(str).unique().tolist())
        raise ValueError(f"Unrecognized Sex values: {values}")
    return mapped


def load_and_merge_inputs(
    harmonized_file: Path,
    tau_file: Path,
    vbm_file: Path,
    metadata_file: Path,
) -> pd.DataFrame:
    for path, label in (
        (harmonized_file, "harmonized N145 input"),
        (tau_file, "amyloid/demographic/tau input"),
        (vbm_file, "VBM input"),
        (metadata_file, "ComBat metadata input"),
    ):
        require_file(path, label)

    harmonized = normalize_ptid(pd.read_excel(harmonized_file), "harmonized input")
    tau = normalize_ptid(pd.read_excel(tau_file), "amyloid/demographic/tau input")
    vbm = normalize_ptid(pd.read_csv(vbm_file), "VBM input")
    metadata = normalize_ptid(pd.read_csv(metadata_file), "ComBat metadata input")

    tau_columns = [
        "PTID",
        "Age",
        "Sex",
        "CL_pvc",
        "tau_mesial_pvc",
        "tau_metatemporal_pvc",
        "tau_temporoparietal_pvc",
    ]
    vbm_columns = ["PTID", "Mean_GMV_HIP_LH", "Mean_GMV_HIP_RH"]
    require_columns(harmonized, ["PTID", "Group", PREDICTOR], "harmonized input")
    require_columns(tau, tau_columns, "amyloid/demographic/tau input")
    require_columns(vbm, vbm_columns, "VBM input")
    require_columns(metadata, ["PTID", "edu"], "ComBat metadata input")

    group_counts = harmonized["Group"].value_counts().to_dict()
    if group_counts != EXPECTED_GROUP_COUNTS:
        raise ValueError(
            "The harmonized input is not the canonical N145 cohort. "
            f"Expected {EXPECTED_GROUP_COUNTS}; found {group_counts}."
        )

    merged = harmonized[["PTID", "Group", PREDICTOR]].merge(
        tau[tau_columns], on="PTID", how="left", validate="one_to_one", indicator="_tau"
    )
    missing = merged.loc[merged["_tau"] != "both", "PTID"].tolist()
    if missing:
        raise ValueError(f"PTIDs missing from amyloid/tau source: {missing}")
    merged = merged.drop(columns="_tau")

    merged = merged.merge(
        vbm[vbm_columns], on="PTID", how="left", validate="one_to_one", indicator="_vbm"
    )
    missing = merged.loc[merged["_vbm"] != "both", "PTID"].tolist()
    if missing:
        raise ValueError(f"PTIDs missing from VBM source: {missing}")
    merged = merged.drop(columns="_vbm")

    merged = merged.merge(
        metadata[["PTID", "edu"]],
        on="PTID",
        how="left",
        validate="one_to_one",
        indicator="_metadata",
    )
    missing = merged.loc[merged["_metadata"] != "both", "PTID"].tolist()
    if missing:
        raise ValueError(f"PTIDs missing from ComBat metadata: {missing}")
    merged = merged.drop(columns="_metadata")

    numeric_columns = [
        PREDICTOR,
        "Age",
        "edu",
        "CL_pvc",
        "tau_mesial_pvc",
        "tau_metatemporal_pvc",
        "tau_temporoparietal_pvc",
        "Mean_GMV_HIP_LH",
        "Mean_GMV_HIP_RH",
    ]
    for column in numeric_columns:
        merged[column] = to_numeric(merged[column], column)

    merged["AGE"] = merged["Age"]
    merged["SEX_NUM"] = encode_sex(merged["Sex"])
    merged["EDUCATION"] = merged["edu"]
    merged["Mean_GMV_HIP_BI"] = merged[
        ["Mean_GMV_HIP_LH", "Mean_GMV_HIP_RH"]
    ].mean(axis=1, skipna=False)

    complete_columns = [
        PREDICTOR,
        "AGE",
        "SEX_NUM",
        "EDUCATION",
        "CL_pvc",
        "Mean_GMV_HIP_BI",
    ]
    missing_complete = merged[complete_columns].isna().sum()
    if missing_complete.any():
        raise ValueError(
            "Unexpected missing values in complete N145 variables: "
            f"{missing_complete[missing_complete > 0].to_dict()}"
        )

    tau_columns_outcome = [
        "tau_mesial_pvc",
        "tau_metatemporal_pvc",
        "tau_temporoparietal_pvc",
    ]
    tau_complete = merged[tau_columns_outcome].notna()
    if not tau_complete.nunique(axis=1).eq(1).all():
        raise ValueError("The three regional tau outcomes have inconsistent missingness")
    if int(tau_complete.all(axis=1).sum()) != 134:
        raise ValueError(
            "Expected 134 participants with complete regional tau; found "
            f"{int(tau_complete.all(axis=1).sum())}"
        )
    return merged


def model_formula(model_name: str) -> str:
    return "outcome ~ " + " + ".join([PREDICTOR, *MODEL_COVARIATES[model_name]])


def fit_model_set(
    frame: pd.DataFrame,
    model_name: str,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    covariates = MODEL_COVARIATES[model_name]
    results: list[dict[str, object]] = []
    predictions: list[pd.DataFrame] = []

    for outcome in OUTCOMES:
        column = outcome["column"]
        model_frame = frame[[column, PREDICTOR, *covariates]].dropna().copy()
        design = sm.add_constant(model_frame[[PREDICTOR, *covariates]], has_constant="add")
        model = sm.OLS(model_frame[column], design).fit(cov_type="HC3")
        confidence = model.conf_int(alpha=0.05).loc[PREDICTOR]

        result: dict[str, object] = {
            "Model": model_name,
            "Formula": model_formula(model_name),
            "Panel": outcome["panel"],
            "Outcome": column,
            "Outcome_Label": outcome["label"],
            "N": int(model.nobs),
            "Residual_DF": float(model.df_resid),
            "Predictor": PREDICTOR,
            "Beta_Unstandardized": float(model.params[PREDICTOR]),
            "HC3_SE": float(model.bse[PREDICTOR]),
            "CI95_Lower": float(confidence.iloc[0]),
            "CI95_Upper": float(confidence.iloc[1]),
            "T": float(model.tvalues[PREDICTOR]),
            "P_Raw": float(model.pvalues[PREDICTOR]),
            "R2_Full_Model": float(model.rsquared),
            "Adjusted_R2_Full_Model": float(model.rsquared_adj),
            "AIC": float(model.aic),
            "BIC": float(model.bic),
            "Intercept": float(model.params["const"]),
            "Beta_AGE": float(model.params["AGE"]),
            "Beta_SEX_NUM": float(model.params["SEX_NUM"]),
            "Beta_EDUCATION": (
                float(model.params["EDUCATION"])
                if "EDUCATION" in model.params.index
                else np.nan
            ),
            "Reference_AGE": float(model_frame["AGE"].mean()),
            "Reference_SEX_NUM": float(model_frame["SEX_NUM"].mode().iloc[0]),
            "Reference_EDUCATION": (
                float(model_frame["EDUCATION"].mean())
                if "EDUCATION" in model_frame.columns
                else np.nan
            ),
        }
        results.append(result)

        x_values = np.linspace(model_frame[PREDICTOR].min(), model_frame[PREDICTOR].max(), 200)
        prediction_design: dict[str, np.ndarray] = {
            "const": np.ones(len(x_values)),
            PREDICTOR: x_values,
            "AGE": np.repeat(result["Reference_AGE"], len(x_values)),
            "SEX_NUM": np.repeat(result["Reference_SEX_NUM"], len(x_values)),
        }
        if "EDUCATION" in covariates:
            prediction_design["EDUCATION"] = np.repeat(
                result["Reference_EDUCATION"], len(x_values)
            )
        predicted = model.predict(pd.DataFrame(prediction_design)[design.columns])
        predictions.append(
            pd.DataFrame(
                {
                    "Model": model_name,
                    "Panel": outcome["panel"],
                    "Outcome": column,
                    PREDICTOR: x_values,
                    "Predicted_Outcome": predicted.to_numpy(),
                    "Reference_AGE": result["Reference_AGE"],
                    "Reference_SEX_NUM": result["Reference_SEX_NUM"],
                    "Reference_EDUCATION": result["Reference_EDUCATION"],
                }
            )
        )

    results_frame = pd.DataFrame(results)
    rejected, adjusted, _, _ = multipletests(
        results_frame["P_Raw"].to_numpy(), alpha=0.05, method="fdr_bh"
    )
    results_frame["P_FDR_BH_Across_5_Outcomes"] = adjusted
    results_frame["Significant_FDR_0_05"] = rejected
    return results_frame, pd.concat(predictions, ignore_index=True)


def make_comparison(model_outputs: dict[str, pd.DataFrame]) -> pd.DataFrame:
    left = model_outputs["age_sex"].copy()
    right = model_outputs["age_sex_education"].copy()
    columns = [
        "Panel",
        "Outcome",
        "N",
        "Beta_Unstandardized",
        "HC3_SE",
        "P_Raw",
        "P_FDR_BH_Across_5_Outcomes",
        "Significant_FDR_0_05",
        "R2_Full_Model",
        "Adjusted_R2_Full_Model",
    ]
    comparison = left[columns].merge(
        right[columns],
        on=["Panel", "Outcome"],
        suffixes=("_Age_Sex", "_Age_Sex_Education"),
        validate="one_to_one",
    )
    comparison["Delta_Beta_After_Education"] = (
        comparison["Beta_Unstandardized_Age_Sex_Education"]
        - comparison["Beta_Unstandardized_Age_Sex"]
    )
    comparison["Percent_Signed_Beta_Change_After_Education"] = (
        comparison["Delta_Beta_After_Education"]
        / comparison["Beta_Unstandardized_Age_Sex"]
        * 100
    )
    comparison["Percent_Absolute_Beta_Change_After_Education"] = (
        (
            comparison["Beta_Unstandardized_Age_Sex_Education"].abs()
            - comparison["Beta_Unstandardized_Age_Sex"].abs()
        )
        / comparison["Beta_Unstandardized_Age_Sex"].abs()
        * 100
    )
    comparison["Delta_R2_After_Education"] = (
        comparison["R2_Full_Model_Age_Sex_Education"]
        - comparison["R2_Full_Model_Age_Sex"]
    )
    return comparison


def main(cli_args: list[str] | None = None) -> None:
    args = parse_args(cli_args)
    frame = load_and_merge_inputs(
        harmonized_file=args.harmonized_file,
        tau_file=args.tau_file,
        vbm_file=args.vbm_file,
        metadata_file=args.metadata_file,
    )
    requested_models = (
        list(MODEL_COVARIATES) if args.model == "both" else [args.model]
    )
    print(f"Validated canonical N145 cohort; tau-complete N=134.")
    print("Education source: ComBat metadata column 'edu' joined by PTID.")
    for model_name in requested_models:
        sample_sizes = {
            outcome["panel"]: int(
                frame[[outcome["column"], PREDICTOR, *MODEL_COVARIATES[model_name]]]
                .dropna()
                .shape[0]
            )
            for outcome in OUTCOMES
        }
        print(f"{model_name}: formula={model_formula(model_name)}; N={sample_sizes}")
    if args.validate_only:
        print("Validation completed; no models were fitted and no files were written.")
        return

    output_dir = args.output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    result_outputs: dict[str, pd.DataFrame] = {}
    saved: list[str] = []
    for model_name in requested_models:
        results, predictions = fit_model_set(frame, model_name)
        result_outputs[model_name] = results
        result_path = output_dir / f"N145_ATN_regression_{model_name}.csv"
        prediction_path = output_dir / f"N145_ATN_prediction_lines_{model_name}.csv"
        results.to_csv(result_path, index=False)
        predictions.to_csv(prediction_path, index=False)
        saved.extend([str(result_path), str(prediction_path)])
        print(f"\n{model_name} results:")
        for row in results.itertuples(index=False):
            print(
                f"  {row.Panel}: N={row.N}, beta={row.Beta_Unstandardized:.8g}, "
                f"p={row.P_Raw:.8g}, "
                f"pFDR={row.P_FDR_BH_Across_5_Outcomes:.8g}"
            )

    comparison_path: Path | None = None
    if set(result_outputs) == set(MODEL_COVARIATES):
        comparison = make_comparison(result_outputs)
        comparison_path = output_dir / "N145_ATN_regression_model_comparison.csv"
        comparison.to_csv(comparison_path, index=False)
        saved.append(str(comparison_path))

    provenance = {
        "analysis": "N145 turbulence--AT(N) ordinary least-squares regressions",
        "harmonized_turbulence_input": portable_input_path(args.harmonized_file),
        "external_amyloid_demographic_tau_input": portable_input_path(args.tau_file),
        "external_postcombat_vbm_input": portable_input_path(args.vbm_file),
        "education_metadata_input": portable_input_path(args.metadata_file),
        "education_column": "edu",
        "cohort_group_counts": EXPECTED_GROUP_COUNTS,
        "models": {
            name: {
                "formula": model_formula(name),
                "covariates": MODEL_COVARIATES[name],
            }
            for name in requested_models
        },
        "covariance": "HC3 heteroscedasticity-consistent",
        "sex_coding": "female=0, male=1",
        "group_role": "point colour only; group is not a regression covariate",
        "multiple_testing": "BH-FDR across five outcomes separately within each model",
        "beta_definition": f"unstandardized coefficient for {PREDICTOR}",
        "prediction_lines": (
            "age fixed at model-sample mean, sex at model-sample mode, and "
            "education at model-sample mean when included"
        ),
        "participant_level_merged_data_saved": False,
        "aggregate_output_files": saved,
    }
    provenance_path = output_dir / "N145_ATN_regression_provenance.json"
    provenance_path.write_text(json.dumps(provenance, indent=2) + "\n", encoding="utf-8")
    print(f"\nSaved aggregate outputs and provenance under: {output_dir}")


if __name__ == "__main__":
    main()
