#!/usr/bin/env python3
"""Test turbulence--AT(N) associations within and across diagnostic stages.

This sensitivity analysis complements the original pooled regressions. For
each AT(N) outcome it fits: (i) the original pooled model, (ii) four separate
stage-specific models, (iii) a group-adjusted common-slope model, and (iv) a
turbulence-by-group interaction model. All models adjust for age, sex, and
education and use HC3 heteroscedasticity-consistent inference.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import statsmodels.api as sm
from statsmodels.stats.multitest import multipletests

from run_atn_regressions_N145 import (
    DEFAULT_HARMONIZED_FILE,
    DEFAULT_METADATA_FILE,
    DEFAULT_TAU_FILE,
    DEFAULT_VBM_FILE,
    GROUP_ORDER,
    OUTCOMES,
    PREDICTOR,
    SCRIPT_DIR,
    load_and_merge_inputs,
    portable_input_path,
)


GROUP_LABELS = {
    "HC_ABneg": "HC−",
    "HC_ABpos": "HC+",
    "MCI_ABpos": "MCI+",
    "AD_ABpos": "AD+",
}
GROUP_COLORS = {
    "HC_ABneg": "#3B82B8",
    "HC_ABpos": "#60A5FA",
    "MCI_ABpos": "#F59E0B",
    "AD_ABpos": "#C94C4C",
}
COVARIATES = ["AGE", "SEX_NUM", "EDUCATION"]
DEFAULT_OUTPUT_DIR = SCRIPT_DIR / "results" / "N145_ATN_turbulence_stage_sensitivity"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--harmonized-file", type=Path, default=DEFAULT_HARMONIZED_FILE)
    parser.add_argument("--tau-file", type=Path, default=DEFAULT_TAU_FILE)
    parser.add_argument("--vbm-file", type=Path, default=DEFAULT_VBM_FILE)
    parser.add_argument("--metadata-file", type=Path, default=DEFAULT_METADATA_FILE)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    return parser.parse_args()


def add_fdr(frame: pd.DataFrame, p_column: str, output_column: str) -> None:
    valid = frame[p_column].notna()
    frame[output_column] = np.nan
    if valid.any():
        frame.loc[valid, output_column] = multipletests(
            frame.loc[valid, p_column], method="fdr_bh"
        )[1]


def fit_hc3(data: pd.DataFrame, outcome: str, columns: list[str]):
    design = sm.add_constant(data[columns], has_constant="add")
    return sm.OLS(data[outcome], design).fit(cov_type="HC3")


def coefficient_row(model, term: str, **metadata) -> dict[str, object]:
    ci = model.conf_int().loc[term]
    return {
        **metadata,
        "N": int(model.nobs),
        "Beta_per_turbulence_SD": float(model.params[term]),
        "HC3_SE": float(model.bse[term]),
        "CI95_Lower": float(ci.iloc[0]),
        "CI95_Upper": float(ci.iloc[1]),
        "T": float(model.tvalues[term]),
        "P_HC3": float(model.pvalues[term]),
        "R2_Full_Model": float(model.rsquared),
        "Adjusted_R2_Full_Model": float(model.rsquared_adj),
    }


def analyse_outcome(frame: pd.DataFrame, spec: dict[str, str]):
    outcome = spec["column"]
    data = frame[["Group", outcome, PREDICTOR, *COVARIATES]].dropna().copy()
    mean = float(data[PREDICTOR].mean())
    sd = float(data[PREDICTOR].std(ddof=1))
    data["Turbulence_z"] = (data[PREDICTOR] - mean) / sd

    pooled = fit_hc3(data, outcome, ["Turbulence_z", *COVARIATES])
    pooled_row = coefficient_row(
        pooled,
        "Turbulence_z",
        Panel=spec["panel"],
        Outcome=outcome,
        Outcome_Label=spec["label"],
        Model="Pooled_original",
        Group="All",
    )

    within_rows = []
    for group in GROUP_ORDER:
        subset = data.loc[data["Group"] == group].copy()
        model = fit_hc3(subset, outcome, ["Turbulence_z", *COVARIATES])
        within_rows.append(
            coefficient_row(
                model,
                "Turbulence_z",
                Panel=spec["panel"],
                Outcome=outcome,
                Outcome_Label=spec["label"],
                Model="Within_stage",
                Group=group,
            )
        )

    group_dummies = pd.get_dummies(
        data["Group"].astype(pd.CategoricalDtype(GROUP_ORDER)),
        prefix="Group",
        drop_first=True,
        dtype=float,
    )
    common_data = pd.concat([data, group_dummies], axis=1)
    dummy_columns = group_dummies.columns.tolist()
    common = fit_hc3(
        common_data, outcome, ["Turbulence_z", *dummy_columns, *COVARIATES]
    )
    common_row = coefficient_row(
        common,
        "Turbulence_z",
        Panel=spec["panel"],
        Outcome=outcome,
        Outcome_Label=spec["label"],
        Model="Group_adjusted_common_slope",
        Group="All",
    )

    interaction_columns = []
    for dummy in dummy_columns:
        name = f"Turbulence_z_x_{dummy}"
        common_data[name] = common_data["Turbulence_z"] * common_data[dummy]
        interaction_columns.append(name)
    interaction = fit_hc3(
        common_data,
        outcome,
        ["Turbulence_z", *dummy_columns, *interaction_columns, *COVARIATES],
    )
    restrictions = np.zeros((len(interaction_columns), len(interaction.params)))
    for row_index, term in enumerate(interaction_columns):
        restrictions[row_index, interaction.params.index.get_loc(term)] = 1
    joint = interaction.wald_test(restrictions, scalar=True)
    interaction_row = {
        "Panel": spec["panel"],
        "Outcome": outcome,
        "Outcome_Label": spec["label"],
        "N": int(interaction.nobs),
        "Test": "Joint turbulence-by-group interaction",
        "Wald_Chi2": float(joint.statistic),
        "DF": int(len(interaction_columns)),
        "P_HC3": float(joint.pvalue),
    }
    return pooled_row, within_rows, common_row, interaction_row, mean, sd


def make_forest_plot(estimates: pd.DataFrame, output_dir: Path) -> None:
    plot_rows = estimates.loc[
        estimates["Model"].isin(["Pooled_original", "Within_stage", "Group_adjusted_common_slope"])
    ].copy()
    figure, axes = plt.subplots(1, 5, figsize=(16, 4.8), sharey=False)
    order = ["Pooled", "HC−", "HC+", "MCI+", "AD+", "Group-adjusted"]
    for axis, spec in zip(axes, OUTCOMES):
        subset = plot_rows.loc[plot_rows["Outcome"] == spec["column"]].copy()
        subset["Display"] = np.where(
            subset["Model"] == "Pooled_original",
            "Pooled",
            np.where(
                subset["Model"] == "Group_adjusted_common_slope",
                "Group-adjusted",
                subset["Group"].map(GROUP_LABELS),
            ),
        )
        subset = subset.set_index("Display").loc[order].reset_index()
        colors = ["#555555", *[GROUP_COLORS[g] for g in GROUP_ORDER], "#111111"]
        y = np.arange(len(order))[::-1]
        axis.axvline(0, color="#AAAAAA", lw=1, zorder=0)
        for index, color in enumerate(colors):
            estimate = subset.loc[index, "Beta_per_turbulence_SD"]
            axis.errorbar(
                estimate,
                y[index],
                xerr=np.array(
                    [[estimate - subset.loc[index, "CI95_Lower"]],
                     [subset.loc[index, "CI95_Upper"] - estimate]]
                ),
                fmt="none",
                ecolor=color,
                elinewidth=1.5,
                capsize=3,
            )
        axis.scatter(subset["Beta_per_turbulence_SD"], y, c=colors, s=32, zorder=2)
        axis.set_yticks(y, order if axis is axes[0] else [])
        axis.set_title(spec["label"], fontsize=10)
        axis.set_xlabel("β per turbulence SD")
        axis.spines[["top", "right"]].set_visible(False)
    figure.suptitle("Turbulence–AT(N) associations: pooled and stage-specific estimates")
    figure.tight_layout()
    for extension in ("pdf", "png"):
        figure.savefig(
            output_dir / f"N145_ATN_turbulence_stage_specific_forest.{extension}",
            dpi=300,
            bbox_inches="tight",
        )
    plt.close(figure)


def main() -> None:
    args = parse_args()
    frame = load_and_merge_inputs(
        args.harmonized_file, args.tau_file, args.vbm_file, args.metadata_file
    )
    pooled_rows, within_rows, common_rows, interaction_rows = [], [], [], []
    scaling_rows = []
    for spec in OUTCOMES:
        pooled, within, common, interaction, mean, sd = analyse_outcome(frame, spec)
        pooled_rows.append(pooled)
        within_rows.extend(within)
        common_rows.append(common)
        interaction_rows.append(interaction)
        scaling_rows.append({"Outcome": spec["column"], "Turbulence_mean": mean, "Turbulence_SD": sd})

    pooled = pd.DataFrame(pooled_rows)
    within = pd.DataFrame(within_rows)
    common = pd.DataFrame(common_rows)
    interactions = pd.DataFrame(interaction_rows)
    add_fdr(pooled, "P_HC3", "P_FDR_BH_Across_5_Outcomes")
    add_fdr(common, "P_HC3", "P_FDR_BH_Across_5_Outcomes")
    add_fdr(interactions, "P_HC3", "P_FDR_BH_Across_5_Outcomes")
    add_fdr(within, "P_HC3", "P_FDR_BH_Across_20_Stage_Outcome_Tests")
    within["P_FDR_BH_Within_Outcome_4_Groups"] = np.nan
    for outcome, indexes in within.groupby("Outcome").groups.items():
        within.loc[indexes, "P_FDR_BH_Within_Outcome_4_Groups"] = multipletests(
            within.loc[indexes, "P_HC3"], method="fdr_bh"
        )[1]

    estimates = pd.concat([pooled, within, common], ignore_index=True)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    estimates.to_csv(args.output_dir / "N145_ATN_turbulence_stage_specific_estimates.csv", index=False)
    interactions.to_csv(args.output_dir / "N145_ATN_turbulence_group_interactions.csv", index=False)
    pd.DataFrame(scaling_rows).to_csv(args.output_dir / "N145_ATN_turbulence_scaling.csv", index=False)
    make_forest_plot(estimates, args.output_dir)

    provenance = {
        "analysis": "N145 turbulence--AT(N) stage-specific sensitivity analysis",
        "models": {
            "pooled": "outcome ~ z(turbulence) + age + sex + education",
            "within_stage": "same model fitted separately in HC-, HC+, MCI+, and AD+",
            "group_adjusted": "outcome ~ z(turbulence) + group + age + sex + education",
            "interaction": "outcome ~ z(turbulence) * group + age + sex + education",
        },
        "covariance": "HC3 heteroscedasticity-consistent",
        "harmonized_input": portable_input_path(args.harmonized_file),
        "tau_input": portable_input_path(args.tau_file),
        "vbm_input": portable_input_path(args.vbm_file),
        "metadata_input": portable_input_path(args.metadata_file),
        "participant_level_merged_data_saved": False,
    }
    (args.output_dir / "N145_ATN_turbulence_stage_sensitivity_provenance.json").write_text(
        json.dumps(provenance, indent=2) + "\n", encoding="utf-8"
    )

    print("\nPooled versus group-adjusted common slopes:")
    print(pd.concat([pooled, common]).loc[:, ["Outcome_Label", "Model", "N", "Beta_per_turbulence_SD", "P_HC3", "P_FDR_BH_Across_5_Outcomes"]].to_string(index=False))
    print("\nWithin-stage slopes:")
    print(within.loc[:, ["Outcome_Label", "Group", "N", "Beta_per_turbulence_SD", "P_HC3", "P_FDR_BH_Across_20_Stage_Outcome_Tests"]].to_string(index=False))
    print("\nJoint interaction tests:")
    print(interactions.loc[:, ["Outcome_Label", "Wald_Chi2", "DF", "P_HC3", "P_FDR_BH_Across_5_Outcomes"]].to_string(index=False))
    print(f"\nSaved aggregate results and figure under: {args.output_dir.resolve()}")


if __name__ == "__main__":
    main()
