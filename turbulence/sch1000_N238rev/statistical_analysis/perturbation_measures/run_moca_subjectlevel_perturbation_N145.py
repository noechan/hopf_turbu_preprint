#!/usr/bin/env python3
"""Fit education-adjusted MOCA regressions for subject-level Hopf measures.

The two participant-level predictors are tested in separate multiple linear
regression models, matching the legacy estimand while adding education and
using the canonical N145 inputs. Two scientifically distinct specifications
are reported:

1. Total association: MOCA ~ z(measure) + age + sex + education.
2. Beyond-stage association: the same model additionally adjusted for the
   four-level amyloid-status group.

HC3 heteroscedasticity-consistent inference is used for coefficients. BH-FDR
is applied across information capacity and susceptibility separately within
each specification. Participant-level merged data remain in memory.
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
from statsmodels.stats.outliers_influence import variance_inflation_factor


SCRIPT_DIR: Final = Path(__file__).resolve().parent
STATISTICAL_DIR: Final = SCRIPT_DIR.parent
SCH1000_ROOT: Final = STATISTICAL_DIR.parent
DEFAULT_ADNI3_ROOT: Final = Path(
    os.environ.get(
        "ADNI3_ROOT",
        "/path/to/ADNI3",
    )
)
DEFAULT_HOPF_FILE: Final = (
    SCH1000_ROOT
    / "data_export"
    / "ML_Input_ADNI3_4STAGINGBYABETA_ComBat_N145_with_infocap_suscep_sch1000.csv"
)
DEFAULT_CLINICAL_FILE: Final = (
    SCH1000_ROOT
    / "visualization"
    / "python"
    / "data"
    / "Turbu_ComBat_clin_ADNI3_allfeatures_N145.csv"
)
DEFAULT_METADATA_FILE: Final = (
    SCH1000_ROOT
    / "data"
    / "covariates"
    / "covariates_ADNI3_ABeta_N152.csv"
)
DEFAULT_OUTPUT_DIR: Final = (
    SCRIPT_DIR / "results" / "N145_MOCA_subjectlevel_perturbation"
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
MODEL_DEFINITIONS: Final = {
    "Total_association": {
        "include_group": False,
        "estimand": "Association adjusted for age, sex, and education",
    },
    "Group_adjusted": {
        "include_group": True,
        "estimand": "Association beyond four-level amyloid-status stage",
    },
}


def parse_args(cli_args: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--hopf-file",
        type=Path,
        default=DEFAULT_HOPF_FILE,
        help="N145 PTID/Group/Info_Cap/Susceptibility table",
    )
    parser.add_argument(
        "--clinical-file",
        type=Path,
        default=DEFAULT_CLINICAL_FILE,
        help="N145 clinical table containing PTID, Group, and MOCA",
    )
    parser.add_argument(
        "--metadata-file",
        type=Path,
        default=DEFAULT_METADATA_FILE,
        help="Canonical ComBat metadata containing age, gender, and edu",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=DEFAULT_OUTPUT_DIR,
        help="Destination for aggregate model outputs",
    )
    parser.add_argument(
        "--model-set",
        choices=("both", "total", "group-adjusted"),
        default="both",
        help="Regression specification(s) to fit (default: both)",
    )
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
    result = pd.to_numeric(series, errors="coerce")
    unexpected = series.notna() & result.isna()
    if unexpected.any():
        values = sorted(series.loc[unexpected].astype(str).unique().tolist())
        raise ValueError(f"Could not convert {label} to numeric: {values[:10]}")
    return result.astype(float)


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


def load_analysis_data(args: argparse.Namespace) -> tuple[pd.DataFrame, dict[str, int]]:
    for path, label in (
        (args.hopf_file, "participant-level Hopf input"),
        (args.clinical_file, "N145 clinical input"),
        (args.metadata_file, "canonical harmonisation metadata"),
    ):
        require_file(path, label)

    hopf = normalize_ptid(pd.read_csv(args.hopf_file), "Hopf input")
    clinical = normalize_ptid(pd.read_csv(args.clinical_file), "Clinical input")
    metadata = normalize_ptid(pd.read_csv(args.metadata_file), "Metadata input")
    require_columns(
        hopf,
        ["PTID", "Group", *PREDICTORS.values()],
        "Hopf input",
    )
    require_columns(clinical, ["PTID", "Group", "MOCA"], "Clinical input")
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
        clinical[["PTID", "Group", "MOCA"]],
        on="PTID",
        how="left",
        validate="one_to_one",
        suffixes=("", "_clinical"),
        indicator="_clinical_match",
    )
    missing = merged.loc[merged["_clinical_match"] != "both", "PTID"].tolist()
    if missing:
        raise ValueError(f"PTIDs missing from clinical input: {missing}")
    merged = merged.drop(columns="_clinical_match")

    merged = merged.merge(
        metadata[["PTID", "Group", "age", "gender", "edu"]],
        on="PTID",
        how="left",
        validate="one_to_one",
        suffixes=("", "_metadata"),
        indicator="_metadata_match",
    )
    missing = merged.loc[merged["_metadata_match"] != "both", "PTID"].tolist()
    if missing:
        raise ValueError(f"PTIDs missing from metadata input: {missing}")
    merged = merged.drop(columns="_metadata_match")

    mismatch = merged.loc[
        (merged["Group"] != merged["Group_clinical"])
        | (merged["Group"] != merged["Group_metadata"]),
        "PTID",
    ].tolist()
    if mismatch:
        raise ValueError(f"Group labels disagree across inputs for: {mismatch}")

    for column in [*PREDICTORS.values(), "MOCA", "age", "gender", "edu"]:
        merged[column] = numeric(merged[column], column)
    if merged[[*PREDICTORS.values(), "age", "gender", "edu"]].isna().any().any():
        missing_counts = merged[
            [*PREDICTORS.values(), "age", "gender", "edu"]
        ].isna().sum()
        raise ValueError(f"Missing required analysis values: {missing_counts.to_dict()}")
    if not set(merged["gender"].unique()).issubset({0.0, 1.0}):
        raise ValueError("Metadata gender must contain only 0 and 1")

    merged["Group"] = pd.Categorical(
        merged["Group"], categories=GROUP_ORDER, ordered=True
    )
    return merged, counts


def selected_models(model_set: str) -> dict[str, dict[str, object]]:
    if model_set == "total":
        return {"Total_association": MODEL_DEFINITIONS["Total_association"]}
    if model_set == "group-adjusted":
        return {"Group_adjusted": MODEL_DEFINITIONS["Group_adjusted"]}
    return dict(MODEL_DEFINITIONS)


def design_matrix(data: pd.DataFrame, predictor: str, include_group: bool) -> pd.DataFrame:
    predictor_sd = data[predictor].std(ddof=1)
    if not np.isfinite(predictor_sd) or predictor_sd == 0:
        raise ValueError(f"Predictor has zero or invalid variance: {predictor}")
    design = pd.DataFrame(
        {
            "Predictor_z": (data[predictor] - data[predictor].mean()) / predictor_sd,
            "Age": data["age"],
            "Sex_male": data["gender"],
            "Education": data["edu"],
        },
        index=data.index,
    )
    if include_group:
        group_dummies = pd.get_dummies(
            data["Group"],
            prefix="Group",
            drop_first=True,
            dtype=float,
        )
        design = pd.concat([design, group_dummies], axis=1)
    return sm.add_constant(design, has_constant="add").astype(float)


def predictor_vif(design: pd.DataFrame) -> float:
    predictor_index = design.columns.get_loc("Predictor_z")
    return float(variance_inflation_factor(design.to_numpy(), predictor_index))


def fit_models(
    data: pd.DataFrame,
    model_set: str,
) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for model_name, definition in selected_models(model_set).items():
        include_group = bool(definition["include_group"])
        for predictor_label, predictor_column in PREDICTORS.items():
            analysis = data[
                [
                    "PTID",
                    "Group",
                    predictor_column,
                    "MOCA",
                    "age",
                    "gender",
                    "edu",
                ]
            ].dropna()
            design = design_matrix(analysis, predictor_column, include_group)
            outcome = analysis["MOCA"].astype(float)
            full_fit = sm.OLS(outcome, design).fit(cov_type="HC3", use_t=True)
            reduced_design = design.drop(columns="Predictor_z")
            reduced_fit = sm.OLS(outcome, reduced_design).fit()
            coefficient = float(full_fit.params["Predictor_z"])
            robust_se = float(full_fit.bse["Predictor_z"])
            confidence_interval = full_fit.conf_int(alpha=0.05).loc["Predictor_z"]
            t_value = float(full_fit.tvalues["Predictor_z"])
            residual_df = float(full_fit.df_resid)
            partial_r2 = t_value**2 / (t_value**2 + residual_df)
            sample_counts = (
                analysis["Group"].value_counts().reindex(GROUP_ORDER, fill_value=0)
            )

            rows.append(
                {
                    "Model": model_name,
                    "Estimand": definition["estimand"],
                    "Predictor": predictor_label,
                    "Source_variable": predictor_column,
                    "N": len(analysis),
                    "HC_ABneg_N": int(sample_counts["HC_ABneg"]),
                    "HC_ABpos_N": int(sample_counts["HC_ABpos"]),
                    "MCI_ABpos_N": int(sample_counts["MCI_ABpos"]),
                    "AD_ABpos_N": int(sample_counts["AD_ABpos"]),
                    "MOCA_mean": float(outcome.mean()),
                    "MOCA_SD": float(outcome.std(ddof=1)),
                    "Predictor_mean": float(analysis[predictor_column].mean()),
                    "Predictor_SD": float(analysis[predictor_column].std(ddof=1)),
                    "Beta_MOCA_per_predictor_SD": coefficient,
                    "HC3_SE": robust_se,
                    "HC3_CI95_low": float(confidence_interval.iloc[0]),
                    "HC3_CI95_high": float(confidence_interval.iloc[1]),
                    "HC3_t": t_value,
                    "DF_residual": int(residual_df),
                    "HC3_P_two_sided": float(full_fit.pvalues["Predictor_z"]),
                    "Partial_R2_from_HC3_t": float(partial_r2),
                    "Full_model_R2": float(full_fit.rsquared),
                    "Adjusted_R2": float(full_fit.rsquared_adj),
                    "Delta_R2_vs_covariates": float(
                        full_fit.rsquared - reduced_fit.rsquared
                    ),
                    "Predictor_VIF": predictor_vif(design),
                }
            )

    results = pd.DataFrame(rows)
    results["P_FDR_BH_across_2_predictors_within_model"] = np.nan
    for model_name, indexes in results.groupby("Model").groups.items():
        p_values = results.loc[indexes, "HC3_P_two_sided"].to_numpy(float)
        results.loc[indexes, "P_FDR_BH_across_2_predictors_within_model"] = (
            multipletests(p_values, alpha=0.05, method="fdr_bh")[1]
        )
    results["Significant_FDR_0_05"] = (
        results["P_FDR_BH_across_2_predictors_within_model"] < 0.05
    )
    return results


def write_results(
    results: pd.DataFrame,
    data: pd.DataFrame,
    args: argparse.Namespace,
) -> tuple[Path, Path]:
    args.output_dir.mkdir(parents=True, exist_ok=True)
    results_file = (
        args.output_dir
        / "N145_MOCA_subjectlevel_perturbation_HC3_age_sex_education.csv"
    )
    results.to_csv(results_file, index=False)

    missing_moca = data.loc[data["MOCA"].isna(), "PTID"].astype(str).tolist()
    predictor_correlation = float(
        data[["Info_Cap", "Susceptibility"]].corr().iloc[0, 1]
    )
    provenance = {
        "participant_level_hopf_input": portable_path(args.hopf_file),
        "clinical_input": portable_path(args.clinical_file),
        "clinical_columns_used": ["PTID", "Group", "MOCA"],
        "metadata_input": portable_path(args.metadata_file),
        "metadata_columns_used": ["PTID", "Group", "age", "gender", "edu"],
        "outcome": "MOCA",
        "predictors_tested_separately": PREDICTORS,
        "predictor_pearson_correlation": predictor_correlation,
        "joint_model_rationale": (
            "The two predictors are not entered jointly because their near-perfect "
            "correlation would make predictor-specific coefficients unstable."
        ),
        "total_association_model": "MOCA ~ z(measure) + age + sex + education",
        "group_adjusted_model": (
            "MOCA ~ z(measure) + age + sex + education + amyloid-status group"
        ),
        "coefficient_inference": "OLS with HC3 heteroscedasticity-consistent SE",
        "multiplicity": (
            "BH-FDR across information capacity and susceptibility separately "
            "within each model specification"
        ),
        "canonical_n": int(len(data)),
        "moca_complete_n": int(data["MOCA"].notna().sum()),
        "moca_missing_n": int(data["MOCA"].isna().sum()),
        "moca_missing_ptids": missing_moca,
        "group_reference": "HC_ABneg",
        "sex_coding": "0=Female, 1=Male",
    }
    provenance_file = args.output_dir / "N145_MOCA_subjectlevel_perturbation_provenance.json"
    provenance_file.write_text(json.dumps(provenance, indent=2) + "\n")
    return results_file, provenance_file


def main(cli_args: list[str] | None = None) -> None:
    args = parse_args(cli_args)
    data, counts = load_analysis_data(args)
    complete_counts = (
        data.loc[data["MOCA"].notna(), "Group"]
        .value_counts()
        .reindex(GROUP_ORDER, fill_value=0)
        .to_dict()
    )
    print(
        "Validated canonical N145 participant-level perturbation cohort: "
        + " | ".join(f"{group}={counts[group]}" for group in GROUP_ORDER)
    )
    print(
        f"MOCA-complete sample: N={data['MOCA'].notna().sum()} "
        f"(missing={data['MOCA'].isna().sum()}) | "
        + " | ".join(f"{group}={complete_counts[group]}" for group in GROUP_ORDER)
    )
    predictor_correlation = data[["Info_Cap", "Susceptibility"]].corr().iloc[0, 1]
    print(
        "Information capacity--susceptibility Pearson r: "
        f"{predictor_correlation:.6f}"
    )
    if args.validate_only:
        print("Validation completed; no models were fitted and no files were written.")
        return

    results = fit_models(data, args.model_set)
    results_file, provenance_file = write_results(results, data, args)
    columns = [
        "Model",
        "Predictor",
        "N",
        "Beta_MOCA_per_predictor_SD",
        "HC3_CI95_low",
        "HC3_CI95_high",
        "HC3_P_two_sided",
        "P_FDR_BH_across_2_predictors_within_model",
        "Delta_R2_vs_covariates",
    ]
    print("\nMOCA regression results:")
    print(results[columns].to_string(index=False))
    print(f"\nSaved aggregate results: {results_file}")
    print(f"Saved provenance: {provenance_file}")


if __name__ == "__main__":
    main()
