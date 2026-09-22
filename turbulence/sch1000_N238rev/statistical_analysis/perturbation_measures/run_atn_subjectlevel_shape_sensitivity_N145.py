#!/usr/bin/env python3
"""Test group-specific linear and pooled quadratic Figure 4 associations.

This exploratory sensitivity analysis leaves the prespecified pooled linear
models unchanged.  For the four AT(N) outcomes displayed against information
capability in manuscript Figure 4, it fits:

1. a covariate-adjusted HC3 linear model separately in each diagnostic group;
2. a pooled group-by-information-capability interaction model; and
3. pooled covariate-adjusted linear and quadratic models.

Only aggregate coefficients, model comparisons, prediction lines, and
provenance are saved.  The protected participant-level merged table remains
in memory.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Final

import numpy as np
import pandas as pd
import statsmodels.api as sm
from statsmodels.stats.multitest import multipletests

import run_atn_subjectlevel_perturbation_N145 as atn


SCRIPT_DIR: Final = Path(__file__).resolve().parent
DEFAULT_OUTPUT_DIR: Final = (
    SCRIPT_DIR / "results" / "N145_ATN_subjectlevel_shape_sensitivity"
)
PREDICTOR: Final = "Info_Cap"
PREDICTOR_LABEL: Final = "Information capability"
COVARIATES: Final = ["AGE", "SEX_NUM", "EDUCATION"]
FIGURE4_OUTCOMES: Final = [
    {"panel": "Fig4c", "column": "CL_pvc", "label": "Amyloid burden (CL-PVC)"},
    {"panel": "Fig4d", "column": "tau_mesial_pvc", "label": "Mesial tau burden (SUVR)"},
    {
        "panel": "Fig4e",
        "column": "tau_metatemporal_pvc",
        "label": "Metatemporal tau burden (SUVR)",
    },
    {
        "panel": "Fig4f",
        "column": "tau_temporoparietal_pvc",
        "label": "Temporoparietal tau burden (SUVR)",
    },
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--hopf-file", type=Path, default=atn.DEFAULT_HOPF_FILE)
    parser.add_argument("--tau-file", type=Path, default=atn.DEFAULT_TAU_FILE)
    parser.add_argument("--vbm-file", type=Path, default=atn.DEFAULT_VBM_FILE)
    parser.add_argument("--metadata-file", type=Path, default=atn.DEFAULT_METADATA_FILE)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument(
        "--validate-only",
        action="store_true",
        help="Validate inputs and samples without fitting or writing outputs.",
    )
    return parser.parse_args()


def reference_covariates(frame: pd.DataFrame) -> dict[str, float]:
    return {
        "AGE": float(frame["AGE"].mean()),
        "SEX_NUM": float(frame["SEX_NUM"].mode().iloc[0]),
        "EDUCATION": float(frame["EDUCATION"].mean()),
    }


def standardize_predictor(frame: pd.DataFrame) -> tuple[pd.DataFrame, float, float]:
    result = frame.copy()
    mean = float(result[PREDICTOR].mean())
    sd = float(result[PREDICTOR].std(ddof=1))
    if not np.isfinite(sd) or sd <= 0:
        raise ValueError("Information capability has zero or invalid variance.")
    result["PREDICTOR_Z"] = (result[PREDICTOR] - mean) / sd
    result["PREDICTOR_Z2"] = result["PREDICTOR_Z"] ** 2
    return result, mean, sd


def fit_group_models(
    frame: pd.DataFrame,
    outcome: dict[str, str],
) -> tuple[list[dict[str, object]], list[pd.DataFrame], dict[str, object]]:
    columns = ["Group", PREDICTOR, outcome["column"], *COVARIATES]
    complete = frame[columns].dropna().copy()
    complete, predictor_mean, predictor_sd = standardize_predictor(complete)
    reference = reference_covariates(complete)
    rows: list[dict[str, object]] = []
    predictions: list[pd.DataFrame] = []

    for group in atn.GROUP_ORDER:
        group_frame = complete.loc[complete["Group"] == group].copy()
        design = sm.add_constant(
            group_frame[["PREDICTOR_Z", *COVARIATES]], has_constant="add"
        )
        if np.linalg.matrix_rank(design.to_numpy(dtype=float)) != design.shape[1]:
            raise ValueError(f"Rank-deficient within-group design for {group}/{outcome['column']}.")
        model = sm.OLS(group_frame[outcome["column"]], design).fit(cov_type="HC3")
        confidence = model.conf_int(alpha=0.05).loc["PREDICTOR_Z"]
        reduced_design = sm.add_constant(group_frame[COVARIATES], has_constant="add")
        reduced = sm.OLS(group_frame[outcome["column"]], reduced_design).fit()
        delta_r2 = float(model.rsquared - reduced.rsquared)
        partial_r2 = float(delta_r2 / (1.0 - reduced.rsquared))
        rows.append(
            {
                "Panel": outcome["panel"],
                "Outcome": outcome["column"],
                "Outcome_Label": outcome["label"],
                "Group": group,
                "N": int(model.nobs),
                "Predictor_Mean_Pooled": predictor_mean,
                "Predictor_SD_Pooled": predictor_sd,
                "Beta_Per_Pooled_SD": float(model.params["PREDICTOR_Z"]),
                "Beta_Per_Raw_Unit": float(model.params["PREDICTOR_Z"] / predictor_sd),
                "HC3_SE_Per_Pooled_SD": float(model.bse["PREDICTOR_Z"]),
                "CI95_Lower_Per_Pooled_SD": float(confidence.iloc[0]),
                "CI95_Upper_Per_Pooled_SD": float(confidence.iloc[1]),
                "T_HC3": float(model.tvalues["PREDICTOR_Z"]),
                "P_Raw": float(model.pvalues["PREDICTOR_Z"]),
                "R2_Full_Model": float(model.rsquared),
                "Adjusted_R2_Full_Model": float(model.rsquared_adj),
                "Delta_R2_Predictor": delta_r2,
                "Partial_R2_Predictor": partial_r2,
                "Predictor_Min": float(group_frame[PREDICTOR].min()),
                "Predictor_Max": float(group_frame[PREDICTOR].max()),
            }
        )

        raw_values = np.linspace(
            float(group_frame[PREDICTOR].min()),
            float(group_frame[PREDICTOR].max()),
            200,
        )
        prediction_design = pd.DataFrame(
            {
                "const": np.ones(raw_values.size),
                "PREDICTOR_Z": (raw_values - predictor_mean) / predictor_sd,
                **{name: np.repeat(value, raw_values.size) for name, value in reference.items()},
            }
        )
        predicted = model.predict(prediction_design[design.columns])
        predictions.append(
            pd.DataFrame(
                {
                    "Panel": outcome["panel"],
                    "Outcome": outcome["column"],
                    "Group": group,
                    "Predictor_Raw_Value": raw_values,
                    "Predicted_Outcome": predicted.to_numpy(),
                    **{f"Reference_{name}": value for name, value in reference.items()},
                }
            )
        )

    group_dummies = pd.get_dummies(
        pd.Categorical(complete["Group"], categories=atn.GROUP_ORDER),
        prefix="Group",
        drop_first=True,
        dtype=float,
    )
    interaction = group_dummies.mul(complete["PREDICTOR_Z"].to_numpy(), axis=0)
    interaction.columns = [f"PREDICTOR_Z_x_{column}" for column in group_dummies.columns]
    parallel_design = pd.concat(
        [
            pd.Series(1.0, index=complete.index, name="const"),
            complete[["PREDICTOR_Z"]],
            group_dummies.set_axis(complete.index),
            complete[COVARIATES],
        ],
        axis=1,
    )
    interaction_design = pd.concat(
        [
            parallel_design.drop(columns=COVARIATES),
            interaction.set_axis(complete.index),
            complete[COVARIATES],
        ],
        axis=1,
    )
    parallel_model = sm.OLS(
        complete[outcome["column"]], parallel_design
    ).fit(cov_type="HC3")
    interaction_model = sm.OLS(
        complete[outcome["column"]], interaction_design
    ).fit(cov_type="HC3")
    restriction = np.zeros((interaction.shape[1], interaction_design.shape[1]))
    for row_index, column in enumerate(interaction.columns):
        restriction[row_index, interaction_design.columns.get_loc(column)] = 1.0
    joint_test = interaction_model.wald_test(restriction, scalar=True)
    parallel_confidence = parallel_model.conf_int(alpha=0.05).loc["PREDICTOR_Z"]
    interaction_row = {
        "Panel": outcome["panel"],
        "Outcome": outcome["column"],
        "Outcome_Label": outcome["label"],
        "N": int(interaction_model.nobs),
        "Test": "joint_group_by_information_capability_interaction",
        "Parallel_Slope_Beta_Per_Pooled_SD": float(
            parallel_model.params["PREDICTOR_Z"]
        ),
        "Parallel_Slope_HC3_SE": float(parallel_model.bse["PREDICTOR_Z"]),
        "Parallel_Slope_CI95_Lower": float(parallel_confidence.iloc[0]),
        "Parallel_Slope_CI95_Upper": float(parallel_confidence.iloc[1]),
        "Parallel_Slope_HC3_P_Raw": float(
            parallel_model.pvalues["PREDICTOR_Z"]
        ),
        "Parallel_Slope_R2": float(parallel_model.rsquared),
        "Parallel_Slope_Adjusted_R2": float(parallel_model.rsquared_adj),
        "Wald_Chi2": float(np.asarray(joint_test.statistic).squeeze()),
        "DF": int(interaction.shape[1]),
        "P_Raw": float(np.asarray(joint_test.pvalue).squeeze()),
        "R2_Full_Model": float(interaction_model.rsquared),
        "Adjusted_R2_Full_Model": float(interaction_model.rsquared_adj),
        "Delta_R2_Interactions": float(
            interaction_model.rsquared - parallel_model.rsquared
        ),
        "Delta_AIC_Interactions_minus_Parallel": float(
            interaction_model.aic - parallel_model.aic
        ),
        "Delta_BIC_Interactions_minus_Parallel": float(
            interaction_model.bic - parallel_model.bic
        ),
    }
    return rows, predictions, interaction_row


def fit_pooled_shape_models(
    frame: pd.DataFrame,
    outcome: dict[str, str],
) -> tuple[dict[str, object], pd.DataFrame]:
    columns = [PREDICTOR, outcome["column"], *COVARIATES]
    complete = frame[columns].dropna().copy()
    complete, predictor_mean, predictor_sd = standardize_predictor(complete)
    reference = reference_covariates(complete)

    linear_design = sm.add_constant(
        complete[["PREDICTOR_Z", *COVARIATES]], has_constant="add"
    )
    quadratic_design = sm.add_constant(
        complete[["PREDICTOR_Z", "PREDICTOR_Z2", *COVARIATES]],
        has_constant="add",
    )
    linear = sm.OLS(complete[outcome["column"]], linear_design).fit(cov_type="HC3")
    quadratic = sm.OLS(
        complete[outcome["column"]], quadratic_design
    ).fit(cov_type="HC3")
    quadratic_confidence = quadratic.conf_int(alpha=0.05).loc["PREDICTOR_Z2"]

    # Polynomial fits can be dominated by a single extreme predictor value.
    # Report a prespecified |z| <= 3 diagnostic without replacing the full fit.
    central = complete.loc[complete["PREDICTOR_Z"].abs() <= 3].copy()
    central_design = sm.add_constant(
        central[["PREDICTOR_Z", "PREDICTOR_Z2", *COVARIATES]],
        has_constant="add",
    )
    central_quadratic = sm.OLS(
        central[outcome["column"]], central_design
    ).fit(cov_type="HC3")
    influence = sm.OLS(complete[outcome["column"]], quadratic_design).fit().get_influence()

    result = {
        "Panel": outcome["panel"],
        "Outcome": outcome["column"],
        "Outcome_Label": outcome["label"],
        "N": int(quadratic.nobs),
        "Predictor_Mean": predictor_mean,
        "Predictor_SD": predictor_sd,
        "Linear_Beta_Per_SD": float(linear.params["PREDICTOR_Z"]),
        "Linear_HC3_P": float(linear.pvalues["PREDICTOR_Z"]),
        "Linear_R2": float(linear.rsquared),
        "Linear_Adjusted_R2": float(linear.rsquared_adj),
        "Linear_AIC": float(linear.aic),
        "Linear_BIC": float(linear.bic),
        "Quadratic_Linear_Beta_Per_SD": float(quadratic.params["PREDICTOR_Z"]),
        "Quadratic_Beta_Per_SD2": float(quadratic.params["PREDICTOR_Z2"]),
        "Quadratic_HC3_SE": float(quadratic.bse["PREDICTOR_Z2"]),
        "Quadratic_CI95_Lower": float(quadratic_confidence.iloc[0]),
        "Quadratic_CI95_Upper": float(quadratic_confidence.iloc[1]),
        "Quadratic_HC3_P_Raw": float(quadratic.pvalues["PREDICTOR_Z2"]),
        "Quadratic_R2": float(quadratic.rsquared),
        "Quadratic_Adjusted_R2": float(quadratic.rsquared_adj),
        "Quadratic_AIC": float(quadratic.aic),
        "Quadratic_BIC": float(quadratic.bic),
        "Delta_R2_Quadratic_vs_Linear": float(quadratic.rsquared - linear.rsquared),
        "Partial_R2_Quadratic": float(
            (quadratic.rsquared - linear.rsquared) / (1.0 - linear.rsquared)
        ),
        "Delta_AIC_Quadratic_minus_Linear": float(quadratic.aic - linear.aic),
        "Delta_BIC_Quadratic_minus_Linear": float(quadratic.bic - linear.bic),
        "N_Abs_Predictor_Z_LE_3": int(central_quadratic.nobs),
        "Quadratic_Beta_Abs_Predictor_Z_LE_3": float(
            central_quadratic.params["PREDICTOR_Z2"]
        ),
        "Quadratic_HC3_P_Abs_Predictor_Z_LE_3": float(
            central_quadratic.pvalues["PREDICTOR_Z2"]
        ),
        "Maximum_Cooks_D_Quadratic": float(np.max(influence.cooks_distance[0])),
        "Maximum_Leverage_Quadratic": float(np.max(influence.hat_matrix_diag)),
    }

    raw_values = np.linspace(
        float(complete[PREDICTOR].min()), float(complete[PREDICTOR].max()), 250
    )
    prediction_base = {
        "const": np.ones(raw_values.size),
        "PREDICTOR_Z": (raw_values - predictor_mean) / predictor_sd,
        **{name: np.repeat(value, raw_values.size) for name, value in reference.items()},
    }
    linear_prediction_design = pd.DataFrame(prediction_base)
    quadratic_prediction_design = linear_prediction_design.copy()
    quadratic_prediction_design["PREDICTOR_Z2"] = (
        quadratic_prediction_design["PREDICTOR_Z"] ** 2
    )
    predictions = pd.concat(
        [
            pd.DataFrame(
                {
                    "Panel": outcome["panel"],
                    "Outcome": outcome["column"],
                    "Model": "linear",
                    "Predictor_Raw_Value": raw_values,
                    "Predicted_Outcome": linear.predict(
                        linear_prediction_design[linear_design.columns]
                    ).to_numpy(),
                }
            ),
            pd.DataFrame(
                {
                    "Panel": outcome["panel"],
                    "Outcome": outcome["column"],
                    "Model": "quadratic",
                    "Predictor_Raw_Value": raw_values,
                    "Predicted_Outcome": quadratic.predict(
                        quadratic_prediction_design[quadratic_design.columns]
                    ).to_numpy(),
                }
            ),
        ],
        ignore_index=True,
    )
    for name, value in reference.items():
        predictions[f"Reference_{name}"] = value
    return result, predictions


def apply_multiplicity(
    group_results: pd.DataFrame,
    interaction_results: pd.DataFrame,
    shape_results: pd.DataFrame,
) -> None:
    rejected, adjusted, _, _ = multipletests(
        group_results["P_Raw"].to_numpy(), alpha=0.05, method="fdr_bh"
    )
    group_results["P_FDR_BH_Across_16_Group_Slopes"] = adjusted
    group_results["Significant_FDR_Across_16"] = rejected
    group_results["P_FDR_BH_Within_Outcome_4_Groups"] = np.nan
    for outcome in group_results["Outcome"].unique():
        mask = group_results["Outcome"] == outcome
        _, adjusted_outcome, _, _ = multipletests(
            group_results.loc[mask, "P_Raw"].to_numpy(),
            alpha=0.05,
            method="fdr_bh",
        )
        group_results.loc[mask, "P_FDR_BH_Within_Outcome_4_Groups"] = adjusted_outcome

    for table, p_column, output_column in (
        (interaction_results, "P_Raw", "P_FDR_BH_Across_4_Outcomes"),
        (
            interaction_results,
            "Parallel_Slope_HC3_P_Raw",
            "Parallel_Slope_HC3_P_FDR_BH_Across_4_Outcomes",
        ),
        (
            shape_results,
            "Quadratic_HC3_P_Raw",
            "Quadratic_HC3_P_FDR_BH_Across_4_Outcomes",
        ),
    ):
        rejected, adjusted, _, _ = multipletests(
            table[p_column].to_numpy(), alpha=0.05, method="fdr_bh"
        )
        table[output_column] = adjusted
        table[f"Significant_{output_column}"] = rejected


def main() -> None:
    args = parse_args()
    frame, predictor_correlation = atn.load_analysis_data(args)
    sample_sizes = {
        outcome["panel"]: int(
            frame[[outcome["column"], PREDICTOR, *COVARIATES]].dropna().shape[0]
        )
        for outcome in FIGURE4_OUTCOMES
    }
    print("Validated canonical N145 cohort; tau-complete N=134.")
    print(f"Figure 4 sensitivity sample sizes: {sample_sizes}")
    if args.validate_only:
        print("Validation completed; no models were fitted and no files were written.")
        return

    group_rows: list[dict[str, object]] = []
    group_predictions: list[pd.DataFrame] = []
    interaction_rows: list[dict[str, object]] = []
    shape_rows: list[dict[str, object]] = []
    shape_predictions: list[pd.DataFrame] = []
    for outcome in FIGURE4_OUTCOMES:
        rows, predictions, interaction = fit_group_models(frame, outcome)
        shape, pooled_predictions = fit_pooled_shape_models(frame, outcome)
        group_rows.extend(rows)
        group_predictions.extend(predictions)
        interaction_rows.append(interaction)
        shape_rows.append(shape)
        shape_predictions.append(pooled_predictions)

    group_results = pd.DataFrame(group_rows)
    interaction_results = pd.DataFrame(interaction_rows)
    shape_results = pd.DataFrame(shape_rows)
    apply_multiplicity(group_results, interaction_results, shape_results)

    output_dir = args.output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    paths = {
        "group_slopes": output_dir / "N145_Figure4_information_capability_group_slopes_HC3.csv",
        "group_predictions": output_dir / "N145_Figure4_information_capability_group_prediction_lines.csv",
        "interaction_tests": output_dir / "N145_Figure4_information_capability_group_interactions_HC3.csv",
        "pooled_shape": output_dir / "N145_Figure4_information_capability_pooled_quadratic_HC3.csv",
        "pooled_predictions": output_dir / "N145_Figure4_information_capability_pooled_shape_prediction_lines.csv",
    }
    group_results.to_csv(paths["group_slopes"], index=False)
    pd.concat(group_predictions, ignore_index=True).to_csv(
        paths["group_predictions"], index=False
    )
    interaction_results.to_csv(paths["interaction_tests"], index=False)
    shape_results.to_csv(paths["pooled_shape"], index=False)
    pd.concat(shape_predictions, ignore_index=True).to_csv(
        paths["pooled_predictions"], index=False
    )

    provenance = {
        "analysis": "Exploratory Figure 4 regression-shape sensitivity analysis",
        "primary_analysis_unchanged": True,
        "predictor": {"column": PREDICTOR, "label": PREDICTOR_LABEL},
        "outcomes": FIGURE4_OUTCOMES,
        "models": {
            "within_group": (
                "Outcome ~ pooled-z(Info_Cap) + age + sex + education, fitted "
                "separately in each group with HC3 robust covariance"
            ),
            "slope_heterogeneity": (
                "Outcome ~ pooled-z(Info_Cap) * group + age + sex + education; "
                "joint HC3 Wald test of three interaction terms"
            ),
            "group_adjusted_parallel_slope": (
                "Outcome ~ pooled-z(Info_Cap) + group + age + sex + education; "
                "HC3 inference for the common within-group slope"
            ),
            "pooled_quadratic": (
                "Outcome ~ z(Info_Cap) + z(Info_Cap)^2 + age + sex + education, "
                "with HC3 robust covariance"
            ),
        },
        "multiplicity": {
            "group_slopes": "BH-FDR across all 16 slopes; within-outcome four-group FDR also reported",
            "interaction_tests": "BH-FDR across four Figure 4 outcomes",
            "group_adjusted_parallel_slopes": "BH-FDR across four Figure 4 outcomes",
            "quadratic_terms": "BH-FDR across four Figure 4 outcomes",
        },
        "quadratic_diagnostic": (
            "The quadratic term is additionally refitted after excluding observations "
            "with |standardized information capability| > 3; this is an influence "
            "diagnostic and does not replace the full-sample model."
        ),
        "information_capacity_susceptibility_correlation": predictor_correlation,
        "participant_level_merged_data_saved": False,
        "inputs": {
            "hopf": atn.portable_path(args.hopf_file),
            "amyloid_tau": atn.portable_path(args.tau_file),
            "vbm": atn.portable_path(args.vbm_file),
            "metadata": atn.portable_path(args.metadata_file),
        },
        "outputs": {name: str(path) for name, path in paths.items()},
    }
    provenance_path = output_dir / "N145_Figure4_shape_sensitivity_provenance.json"
    provenance_path.write_text(json.dumps(provenance, indent=2) + "\n", encoding="utf-8")

    print("\nWithin-group slopes:")
    print(
        group_results[
            ["Panel", "Group", "N", "Beta_Per_Pooled_SD", "P_Raw", "P_FDR_BH_Across_16_Group_Slopes"]
        ].to_string(index=False)
    )
    print("\nGroup-by-predictor interaction tests:")
    print(
        interaction_results[
            [
                "Panel",
                "Parallel_Slope_Beta_Per_Pooled_SD",
                "Parallel_Slope_HC3_P_Raw",
                "Parallel_Slope_HC3_P_FDR_BH_Across_4_Outcomes",
                "Wald_Chi2",
                "DF",
                "P_Raw",
                "P_FDR_BH_Across_4_Outcomes",
            ]
        ].to_string(index=False)
    )
    print("\nPooled quadratic tests:")
    print(
        shape_results[
            [
                "Panel",
                "N",
                "Quadratic_Beta_Per_SD2",
                "Quadratic_HC3_P_Raw",
                "Quadratic_HC3_P_FDR_BH_Across_4_Outcomes",
                "Delta_AIC_Quadratic_minus_Linear",
                "Quadratic_HC3_P_Abs_Predictor_Z_LE_3",
            ]
        ].to_string(index=False)
    )
    print(f"\nSaved aggregate sensitivity outputs under: {output_dir}")


if __name__ == "__main__":
    main()
