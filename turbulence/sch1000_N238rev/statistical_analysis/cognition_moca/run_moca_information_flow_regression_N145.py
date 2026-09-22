#!/usr/bin/env python3
"""Fit AT(N)-style N145 brain-dynamics--MOCA regressions.

The model is cross-sectional and estimates the association between MOCA and a
prespecified brain-dynamics measure at lambda=0.01 while adjusting for age,
sex, and education. The default information-flow model is:

    MOCA ~ z(InfoFlow_lam_0_01) + AGE + SEX_NUM + EDUCATION

HC3 heteroscedasticity-consistent inference matches the active AT(N) analysis.
Group is retained for visualization only and is not a model covariate. The
script saves aggregate statistics, diagnostics, prediction lines, and
provenance; participant-level merged data are not written.
"""

from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
from typing import Final

import numpy as np
import pandas as pd
from scipy import stats
import statsmodels.api as sm
from statsmodels.stats.diagnostic import het_breuschpagan, linear_reset
from statsmodels.stats.outliers_influence import variance_inflation_factor


SCRIPT_DIR: Final = Path(__file__).resolve().parent
SCH1000_ROOT: Final = SCRIPT_DIR.parent.parent
DEFAULT_ADNI3_ROOT: Final = Path(
    os.environ.get(
        "ADNI3_ROOT",
        "/path/to/ADNI3",
    )
)
DEFAULT_HARMONIZED_FILE: Final = (
    SCH1000_ROOT
    / "harmonization_allfeat"
    / "recomputed"
    / "Turbu_ComBat_ADNI3_allfeatures_N145.xlsx"
)
LOCAL_CLINICAL_FILE: Final = (
    SCH1000_ROOT
    / "visualization"
    / "python"
    / "data"
    / "Turbu_ComBat_clin_ADNI3_allfeatures_N145.csv"
)
MOUNTED_CLINICAL_FILE: Final = (
    DEFAULT_ADNI3_ROOT
    / "code"
    / "HPC_Hopf_SUB_DTI_1000_Staging"
    / "visualization"
    / "data"
    / "Turbu_ComBat_clin_ADNI3_allfeatures_N145.csv"
)
DEFAULT_CLINICAL_FILE: Final = (
    LOCAL_CLINICAL_FILE if LOCAL_CLINICAL_FILE.is_file() else MOUNTED_CLINICAL_FILE
)
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
PREDICTOR_SPECS: Final = {
    "information-flow": {
        "source": "InfoFlow_lam_0_01",
        "value": "InfoFlow_lam_0_01",
        "z": "InformationFlow_z",
        "slug": "information_flow",
        "label": "Information flow",
        "transform": "identity",
    },
    "turbulence": {
        "source": "Turbu_lam_0_01",
        "value": "Turbu_lam_0_01",
        "z": "Turbulence_z",
        "slug": "turbulence",
        "label": "Turbulence",
        "transform": "identity",
    },
    "one-minus-information-transfer": {
        "source": "InfoTransfer_lam_0_01",
        "value": "OneMinus_InfoTransfer_lam_0_01",
        "z": "OneMinusInformationTransfer_z",
        "slug": "one_minus_information_transfer",
        "label": "1 - Information transfer",
        "transform": "one_minus",
    },
}
GROUP_ORDER: Final = ["HC_ABneg", "HC_ABpos", "MCI_ABpos", "AD_ABpos"]
EXPECTED_GROUP_COUNTS: Final = {
    "HC_ABneg": 51,
    "HC_ABpos": 37,
    "MCI_ABpos": 31,
    "AD_ABpos": 26,
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--predictor",
        choices=tuple(PREDICTOR_SPECS),
        default="information-flow",
        help="Dynamical predictor to fit (default: information-flow)",
    )
    parser.add_argument(
        "--harmonized-file", type=Path, default=DEFAULT_HARMONIZED_FILE
    )
    parser.add_argument("--clinical-file", type=Path, default=DEFAULT_CLINICAL_FILE)
    parser.add_argument("--metadata-file", type=Path, default=DEFAULT_METADATA_FILE)
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=None,
        help="Destination (default: predictor-specific results directory)",
    )
    parser.add_argument(
        "--validate-only",
        action="store_true",
        help="Validate inputs and the N144 model sample without fitting or writing",
    )
    return parser.parse_args()


def require_file(path: Path, label: str) -> None:
    if not path.is_file():
        raise FileNotFoundError(f"Missing {label}: {path}")


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
    converted = pd.to_numeric(series, errors="coerce")
    unexpected = series.notna() & converted.isna()
    if unexpected.any():
        values = sorted(series.loc[unexpected].astype(str).unique().tolist())
        raise ValueError(f"Could not convert {label} to numeric: {values[:10]}")
    return converted


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


def load_analysis_data(
    args: argparse.Namespace,
    predictor_source: str,
    predictor_value: str,
    predictor_z: str,
    predictor_label: str,
    predictor_transform: str,
) -> tuple[pd.DataFrame, list[str]]:
    for path, label in (
        (args.harmonized_file, "harmonized N145 input"),
        (args.clinical_file, "clinical input"),
        (args.metadata_file, "harmonization metadata"),
    ):
        require_file(path, label)

    harmonized = normalize_ptid(pd.read_excel(args.harmonized_file), "harmonized input")
    clinical = normalize_ptid(pd.read_csv(args.clinical_file), "clinical input")
    metadata = normalize_ptid(pd.read_csv(args.metadata_file), "metadata input")
    require_columns(harmonized, ["PTID", "Group", predictor_source], "harmonized input")
    require_columns(clinical, ["PTID", "Group", "MOCA"], "clinical input")
    require_columns(
        metadata,
        ["PTID", "Group", "age", "gender", "edu"],
        "metadata input",
    )

    group_counts = harmonized["Group"].value_counts().to_dict()
    if group_counts != EXPECTED_GROUP_COUNTS:
        raise ValueError(
            "The harmonized input is not canonical N145. "
            f"Expected {EXPECTED_GROUP_COUNTS}; found {group_counts}."
        )

    merged = harmonized[["PTID", "Group", predictor_source]].merge(
        clinical[["PTID", "Group", "MOCA"]].rename(columns={"Group": "Clinical_Group"}),
        on="PTID",
        how="left",
        validate="one_to_one",
        indicator="_clinical",
    )
    missing = merged.loc[merged["_clinical"] != "both", "PTID"].tolist()
    if missing:
        raise ValueError(f"PTIDs missing from clinical input: {missing}")
    merged = merged.drop(columns="_clinical")
    merged = merged.merge(
        metadata[["PTID", "Group", "age", "gender", "edu"]].rename(
            columns={"Group": "Metadata_Group"}
        ),
        on="PTID",
        how="left",
        validate="one_to_one",
        indicator="_metadata",
    )
    missing = merged.loc[merged["_metadata"] != "both", "PTID"].tolist()
    if missing:
        raise ValueError(f"PTIDs missing from metadata input: {missing}")
    merged = merged.drop(columns="_metadata")

    mismatched = merged.loc[
        (merged["Group"] != merged["Clinical_Group"])
        | (merged["Group"] != merged["Metadata_Group"]),
        "PTID",
    ].tolist()
    if mismatched:
        raise ValueError(f"Group labels disagree across inputs for: {mismatched}")

    for column in [predictor_source, "MOCA", "age", "gender", "edu"]:
        merged[column] = to_numeric(merged[column], column)
    if predictor_transform == "one_minus":
        merged[predictor_value] = 1 - merged[predictor_source]
    elif predictor_transform == "identity":
        merged[predictor_value] = merged[predictor_source]
    else:
        raise ValueError(f"Unsupported predictor transform: {predictor_transform}")
    if not set(merged["gender"].dropna().unique()).issubset({0, 1}):
        raise ValueError("Metadata gender must contain only 0 and 1")
    complete_covariates = merged[[predictor_value, "age", "gender", "edu"]].isna().sum()
    if complete_covariates.any():
        raise ValueError(
            "Unexpected missing model inputs: "
            f"{complete_covariates[complete_covariates > 0].to_dict()}"
        )

    missing_moca_ptids = merged.loc[merged["MOCA"].isna(), "PTID"].tolist()
    analysis = merged.loc[merged["MOCA"].notna()].copy()
    observed_counts = analysis["Group"].value_counts().reindex(GROUP_ORDER).to_dict()
    expected_complete = {**EXPECTED_GROUP_COUNTS, "HC_ABpos": 36}
    if observed_counts != expected_complete:
        raise ValueError(
            f"Expected N144 MOCA groups {expected_complete}; found {observed_counts}."
        )
    if len(missing_moca_ptids) != 1:
        raise ValueError(
            "Expected exactly one participant with missing MOCA; found "
            f"{len(missing_moca_ptids)}."
        )

    analysis["AGE"] = analysis["age"]
    analysis["SEX_NUM"] = analysis["gender"].astype(float)
    analysis["EDUCATION"] = analysis["edu"]
    predictor_mean = float(analysis[predictor_value].mean())
    predictor_sd = float(analysis[predictor_value].std(ddof=1))
    if not np.isfinite(predictor_sd) or predictor_sd == 0:
        raise ValueError(f"{predictor_label} has zero or invalid variance")
    analysis[predictor_z] = (analysis[predictor_value] - predictor_mean) / predictor_sd
    return analysis, missing_moca_ptids


def calculate_vif(design: pd.DataFrame) -> pd.DataFrame:
    return pd.DataFrame(
        {
            "Term": design.columns,
            "VIF": [
                variance_inflation_factor(design.to_numpy(), index)
                for index in range(design.shape[1])
            ],
        }
    )


def main() -> None:
    args = parse_args()
    spec = PREDICTOR_SPECS[args.predictor]
    predictor_source = spec["source"]
    predictor_value = spec["value"]
    predictor_z = spec["z"]
    predictor_slug = spec["slug"]
    predictor_label = spec["label"]
    if args.output_dir is None:
        args.output_dir = (
            SCRIPT_DIR / "results" / f"N145_MOCA_{predictor_slug}_regression"
        )
    analysis, missing_moca_ptids = load_analysis_data(
        args,
        predictor_source,
        predictor_value,
        predictor_z,
        predictor_label,
        spec["transform"],
    )
    group_counts = analysis["Group"].value_counts().reindex(GROUP_ORDER).to_dict()
    print(f"Validated MOCA {predictor_slug} sample: N={len(analysis)}; {group_counts}")
    print(f"Participants with missing MOCA: {len(missing_moca_ptids)}")
    if args.validate_only:
        print("Validation completed; no model was fitted and no files were written.")
        return

    model_columns = [predictor_z, "AGE", "SEX_NUM", "EDUCATION"]
    design = sm.add_constant(analysis[model_columns], has_constant="add")
    reduced_design = sm.add_constant(
        analysis[["AGE", "SEX_NUM", "EDUCATION"]], has_constant="add"
    )
    model = sm.OLS(analysis["MOCA"], design).fit(cov_type="HC3")
    reduced = sm.OLS(analysis["MOCA"], reduced_design).fit(cov_type="HC3")
    confidence = model.conf_int(alpha=0.05).loc[predictor_z]
    predictor_beta = float(model.params[predictor_z])
    residual_df = float(model.df_resid)
    predictor_t = float(model.tvalues[predictor_z])
    partial_r2 = (model.rsquared - reduced.rsquared) / (1 - reduced.rsquared)

    result = pd.DataFrame(
        [
            {
                "Model": "ATN_style_age_sex_education",
                "Formula": f"MOCA ~ z({predictor_value}) + AGE + SEX_NUM + EDUCATION",
                "N": int(model.nobs),
                "Residual_DF": residual_df,
                "Outcome": "MOCA",
                "Predictor": predictor_value,
                "Source_variable": predictor_source,
                "Transform": spec["transform"],
                "Predictor_scale": "standardized within complete N144 sample",
                "Beta_MOCA_per_predictor_SD": predictor_beta,
                "HC3_SE": float(model.bse[predictor_z]),
                "CI95_Lower": float(confidence.iloc[0]),
                "CI95_Upper": float(confidence.iloc[1]),
                "T": predictor_t,
                "P_HC3": float(model.pvalues[predictor_z]),
                "R2_Full_Model": float(model.rsquared),
                "Adjusted_R2_Full_Model": float(model.rsquared_adj),
                "R2_Reduced_Covariate_Model": float(reduced.rsquared),
                "Delta_R2_Predictor": float(model.rsquared - reduced.rsquared),
                "Partial_R2_Predictor": partial_r2,
                "AIC": float(model.aic),
                "BIC": float(model.bic),
                "Predictor_mean": float(analysis[predictor_value].mean()),
                "Predictor_SD": float(analysis[predictor_value].std(ddof=1)),
                "Reference_AGE": float(analysis["AGE"].mean()),
                "Reference_SEX_NUM": float(analysis["SEX_NUM"].mode().iloc[0]),
                "Reference_EDUCATION": float(analysis["EDUCATION"].mean()),
            }
        ]
    )

    ordinary_model = sm.OLS(analysis["MOCA"], design).fit()
    bp_lm, bp_lm_p, bp_f, bp_f_p = het_breuschpagan(
        ordinary_model.resid, ordinary_model.model.exog
    )
    reset = linear_reset(ordinary_model, power=2, use_f=True)
    shapiro_w, shapiro_p = stats.shapiro(ordinary_model.resid)
    vif = calculate_vif(design)
    diagnostics = pd.DataFrame(
        {
            "Diagnostic": [
                "N",
                "Unique_PTIDs",
                "Design_rank",
                "Design_columns",
                "Predictor_covariate_max_VIF",
                "Breusch_Pagan_LM",
                "Breusch_Pagan_LM_p",
                "Breusch_Pagan_F",
                "Breusch_Pagan_F_p",
                "RESET_F",
                "RESET_p",
                "Shapiro_Wilk_W",
                "Shapiro_Wilk_p",
            ],
            "Value": [
                len(analysis),
                analysis["PTID"].nunique(),
                int(np.linalg.matrix_rank(design.to_numpy())),
                design.shape[1],
                float(vif.loc[vif["Term"] != "const", "VIF"].max()),
                bp_lm,
                bp_lm_p,
                bp_f,
                bp_f_p,
                float(reset.fvalue),
                float(reset.pvalue),
                shapiro_w,
                shapiro_p,
            ],
        }
    )

    raw_values = np.linspace(
        analysis[predictor_value].min(), analysis[predictor_value].max(), 200
    )
    predictor_mean = result.loc[0, "Predictor_mean"]
    predictor_sd = result.loc[0, "Predictor_SD"]
    prediction_design = pd.DataFrame(
        {
            "const": np.ones(len(raw_values)),
            predictor_z: (raw_values - predictor_mean) / predictor_sd,
            "AGE": np.repeat(result.loc[0, "Reference_AGE"], len(raw_values)),
            "SEX_NUM": np.repeat(
                result.loc[0, "Reference_SEX_NUM"], len(raw_values)
            ),
            "EDUCATION": np.repeat(
                result.loc[0, "Reference_EDUCATION"], len(raw_values)
            ),
        }
    )
    prediction = pd.DataFrame(
        {
            predictor_value: raw_values,
            predictor_z: prediction_design[predictor_z],
            "Predicted_MOCA": model.predict(prediction_design[design.columns]),
            "Reference_AGE": result.loc[0, "Reference_AGE"],
            "Reference_SEX_NUM": result.loc[0, "Reference_SEX_NUM"],
            "Reference_EDUCATION": result.loc[0, "Reference_EDUCATION"],
        }
    )

    args.output_dir.mkdir(parents=True, exist_ok=True)
    file_stem = f"N145_MOCA_{predictor_slug}_HC3_age_sex_education"
    result_path = args.output_dir / f"{file_stem}.csv"
    coefficient_path = (
        args.output_dir
        / f"{file_stem}_coefficients.csv"
    )
    prediction_path = (
        args.output_dir
        / f"{file_stem}_prediction_line.csv"
    )
    diagnostics_path = (
        args.output_dir
        / f"{file_stem}_diagnostics.csv"
    )
    vif_path = (
        args.output_dir / f"{file_stem}_VIF.csv"
    )
    provenance_path = args.output_dir / f"N145_MOCA_{predictor_slug}_provenance.json"
    result.to_csv(result_path, index=False)
    pd.DataFrame(
        {
            "Term": model.params.index,
            "Estimate": model.params.to_numpy(),
            "HC3_SE": model.bse.to_numpy(),
            "T": model.tvalues.to_numpy(),
            "P_HC3": model.pvalues.to_numpy(),
            "CI95_Lower": model.conf_int().iloc[:, 0].to_numpy(),
            "CI95_Upper": model.conf_int().iloc[:, 1].to_numpy(),
        }
    ).to_csv(coefficient_path, index=False)
    prediction.to_csv(prediction_path, index=False)
    diagnostics.to_csv(diagnostics_path, index=False)
    vif.to_csv(vif_path, index=False)
    provenance_path.write_text(
        json.dumps(
            {
                "analysis": f"N145 AT(N)-style {predictor_slug}--MOCA regression",
                "estimand": "Total cross-sectional association across disease stages",
                "formula": result.loc[0, "Formula"],
                "covariance": "HC3 heteroscedasticity-consistent",
                "group_role": "visualization only; not a regression covariate",
                "harmonized_input": portable_path(args.harmonized_file),
                "clinical_input": portable_path(args.clinical_file),
                "clinical_columns_used": ["PTID", "Group", "MOCA"],
                "metadata_input": portable_path(args.metadata_file),
                "metadata_columns_used": ["PTID", "Group", "age", "gender", "edu"],
                "participants_missing_moca": len(missing_moca_ptids),
                "group_counts_complete_MOCA": group_counts,
                "multiplicity": "None: one prespecified outcome and predictor",
            },
            indent=2,
        )
        + "\n",
        encoding="utf-8",
    )

    print(f"\nAT(N)-style {predictor_slug}--MOCA result:")
    print(result.to_string(index=False))
    print("\nVIF:")
    print(vif.to_string(index=False))
    print(f"\nSaved results under: {args.output_dir.resolve()}")


if __name__ == "__main__":
    main()
