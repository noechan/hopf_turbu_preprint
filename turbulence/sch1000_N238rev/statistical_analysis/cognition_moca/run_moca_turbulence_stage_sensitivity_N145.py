#!/usr/bin/env python3
"""Test the turbulence--MoCA association within and across disease stages."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import statsmodels.api as sm
from statsmodels.stats.multitest import multipletests

from run_moca_information_flow_regression_N145 import (
    DEFAULT_CLINICAL_FILE,
    DEFAULT_HARMONIZED_FILE,
    DEFAULT_METADATA_FILE,
    GROUP_ORDER,
    PREDICTOR_SPECS,
    SCRIPT_DIR,
    load_analysis_data,
    portable_path,
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
DEFAULT_OUTPUT_DIR = SCRIPT_DIR / "results" / "N145_MOCA_turbulence_stage_sensitivity"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--harmonized-file", type=Path, default=DEFAULT_HARMONIZED_FILE)
    parser.add_argument("--clinical-file", type=Path, default=DEFAULT_CLINICAL_FILE)
    parser.add_argument("--metadata-file", type=Path, default=DEFAULT_METADATA_FILE)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    return parser.parse_args()


def fit_hc3(data: pd.DataFrame, columns: list[str]):
    design = sm.add_constant(data[columns], has_constant="add")
    return sm.OLS(data["MOCA"], design).fit(cov_type="HC3")


def coefficient_row(model, term: str, **metadata) -> dict[str, object]:
    ci = model.conf_int().loc[term]
    return {
        **metadata,
        "N": int(model.nobs),
        "Beta_MOCA_per_turbulence_SD": float(model.params[term]),
        "HC3_SE": float(model.bse[term]),
        "CI95_Lower": float(ci.iloc[0]),
        "CI95_Upper": float(ci.iloc[1]),
        "T": float(model.tvalues[term]),
        "P_HC3": float(model.pvalues[term]),
        "R2_Full_Model": float(model.rsquared),
        "Adjusted_R2_Full_Model": float(model.rsquared_adj),
    }


def make_forest_plot(estimates: pd.DataFrame, interaction_p: float, output_dir: Path) -> None:
    order = ["Pooled", "HC−", "HC+", "MCI+", "AD+", "Group-adjusted"]
    estimates = estimates.copy()
    estimates["Display"] = np.where(
        estimates["Model"] == "Pooled_original",
        "Pooled",
        np.where(
            estimates["Model"] == "Group_adjusted_common_slope",
            "Group-adjusted",
            estimates["Group"].map(GROUP_LABELS),
        ),
    )
    estimates = estimates.set_index("Display").loc[order].reset_index()
    colors = ["#555555", *[GROUP_COLORS[g] for g in GROUP_ORDER], "#111111"]
    y = np.arange(len(order))[::-1]
    figure, axis = plt.subplots(figsize=(6.5, 4.6))
    axis.axvline(0, color="#AAAAAA", lw=1, zorder=0)
    for index, color in enumerate(colors):
        estimate = estimates.loc[index, "Beta_MOCA_per_turbulence_SD"]
        axis.errorbar(
            estimate,
            y[index],
            xerr=np.array(
                [[estimate - estimates.loc[index, "CI95_Lower"]],
                 [estimates.loc[index, "CI95_Upper"] - estimate]]
            ),
            fmt="none",
            ecolor=color,
            elinewidth=1.6,
            capsize=3,
        )
    axis.scatter(estimates["Beta_MOCA_per_turbulence_SD"], y, c=colors, s=38, zorder=2)
    axis.set_yticks(y, order)
    axis.set_xlabel("Change in MoCA per turbulence SD (β, 95% CI)")
    axis.set_title(f"Turbulence–MoCA associations\ninteraction $p$ = {interaction_p:.3f}")
    axis.spines[["top", "right"]].set_visible(False)
    figure.tight_layout()
    for extension in ("pdf", "png"):
        figure.savefig(
            output_dir / f"N145_MOCA_turbulence_stage_specific_forest.{extension}",
            dpi=300,
            bbox_inches="tight",
        )
    plt.close(figure)


def main() -> None:
    args = parse_args()
    spec = PREDICTOR_SPECS["turbulence"]
    loader_args = argparse.Namespace(
        harmonized_file=args.harmonized_file,
        clinical_file=args.clinical_file,
        metadata_file=args.metadata_file,
    )
    analysis, missing_moca_ptids = load_analysis_data(
        loader_args,
        spec["source"],
        spec["value"],
        spec["z"],
        spec["label"],
        spec["transform"],
    )
    predictor = spec["z"]

    pooled = fit_hc3(analysis, [predictor, *COVARIATES])
    rows = [
        coefficient_row(
            pooled, predictor, Model="Pooled_original", Group="All"
        )
    ]
    within_rows = []
    for group in GROUP_ORDER:
        subset = analysis.loc[analysis["Group"] == group].copy()
        model = fit_hc3(subset, [predictor, *COVARIATES])
        within_rows.append(
            coefficient_row(model, predictor, Model="Within_stage", Group=group)
        )

    dummies = pd.get_dummies(
        analysis["Group"].astype(pd.CategoricalDtype(GROUP_ORDER)),
        prefix="Group",
        drop_first=True,
        dtype=float,
    )
    model_data = pd.concat([analysis, dummies], axis=1)
    dummy_columns = dummies.columns.tolist()
    common = fit_hc3(model_data, [predictor, *dummy_columns, *COVARIATES])
    rows.extend(within_rows)
    rows.append(
        coefficient_row(
            common, predictor, Model="Group_adjusted_common_slope", Group="All"
        )
    )

    interaction_columns = []
    for dummy in dummy_columns:
        name = f"{predictor}_x_{dummy}"
        model_data[name] = model_data[predictor] * model_data[dummy]
        interaction_columns.append(name)
    interaction = fit_hc3(
        model_data, [predictor, *dummy_columns, *interaction_columns, *COVARIATES]
    )
    restrictions = np.zeros((len(interaction_columns), len(interaction.params)))
    for row_index, term in enumerate(interaction_columns):
        restrictions[row_index, interaction.params.index.get_loc(term)] = 1
    joint = interaction.wald_test(restrictions, scalar=True)
    interaction_result = pd.DataFrame(
        [
            {
                "Test": "Joint turbulence-by-group interaction",
                "N": int(interaction.nobs),
                "Wald_Chi2": float(joint.statistic),
                "DF": len(interaction_columns),
                "P_HC3": float(joint.pvalue),
            }
        ]
    )

    estimates = pd.DataFrame(rows)
    within_mask = estimates["Model"] == "Within_stage"
    estimates["P_FDR_BH_Across_4_Within_Stage_Tests"] = np.nan
    estimates.loc[within_mask, "P_FDR_BH_Across_4_Within_Stage_Tests"] = multipletests(
        estimates.loc[within_mask, "P_HC3"], method="fdr_bh"
    )[1]

    args.output_dir.mkdir(parents=True, exist_ok=True)
    estimates.to_csv(args.output_dir / "N145_MOCA_turbulence_stage_specific_estimates.csv", index=False)
    interaction_result.to_csv(args.output_dir / "N145_MOCA_turbulence_group_interaction.csv", index=False)
    make_forest_plot(estimates, float(joint.pvalue), args.output_dir)
    provenance = {
        "analysis": "N145 turbulence--MoCA stage-specific sensitivity analysis",
        "models": {
            "pooled": "MOCA ~ z(turbulence) + age + sex + education",
            "within_stage": "same model fitted separately in HC-, HC+, MCI+, and AD+",
            "group_adjusted": "MOCA ~ z(turbulence) + group + age + sex + education",
            "interaction": "MOCA ~ z(turbulence) * group + age + sex + education",
        },
        "covariance": "HC3 heteroscedasticity-consistent",
        "harmonized_input": portable_path(args.harmonized_file),
        "clinical_input": portable_path(args.clinical_file),
        "metadata_input": portable_path(args.metadata_file),
        "participants_missing_MOCA": len(missing_moca_ptids),
        "participant_identifiers_saved": False,
    }
    (args.output_dir / "N145_MOCA_turbulence_stage_sensitivity_provenance.json").write_text(
        json.dumps(provenance, indent=2) + "\n", encoding="utf-8"
    )

    print("\nTurbulence--MoCA estimates:")
    print(estimates.loc[:, ["Model", "Group", "N", "Beta_MOCA_per_turbulence_SD", "P_HC3", "P_FDR_BH_Across_4_Within_Stage_Tests"]].to_string(index=False))
    print(f"\nJoint turbulence-by-group interaction: chi2({len(interaction_columns)})={float(joint.statistic):.4f}, p={float(joint.pvalue):.6g}")
    print(f"Saved aggregate results and figure under: {args.output_dir.resolve()}")


if __name__ == "__main__":
    main()
