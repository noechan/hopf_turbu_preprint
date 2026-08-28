#!/usr/bin/env python3
"""Combine the two N145 AT(N)-style MOCA models and apply BH-FDR."""

from __future__ import annotations

from pathlib import Path

import pandas as pd
from statsmodels.stats.multitest import multipletests


SCRIPT_DIR = Path(__file__).resolve().parent
RESULTS_ROOT = SCRIPT_DIR / "results"
INPUTS = [
    RESULTS_ROOT
    / "N145_MOCA_information_flow_regression"
    / "N145_MOCA_information_flow_HC3_age_sex_education.csv",
    RESULTS_ROOT
    / "N145_MOCA_turbulence_regression"
    / "N145_MOCA_turbulence_HC3_age_sex_education.csv",
    RESULTS_ROOT
    / "N145_MOCA_one_minus_information_transfer_regression"
    / "N145_MOCA_one_minus_information_transfer_HC3_age_sex_education.csv",
]
OUTPUT = RESULTS_ROOT / "N145_MOCA_single_predictor_regression_comparison.csv"


def main() -> None:
    missing = [path for path in INPUTS if not path.is_file()]
    if missing:
        raise FileNotFoundError(f"Missing single-predictor results: {missing}")
    frames = [pd.read_csv(path) for path in INPUTS]
    if any(len(frame) != 1 for frame in frames):
        raise ValueError("Each single-predictor result must contain exactly one row")
    comparison = pd.concat(frames, ignore_index=True)
    if comparison["Predictor"].duplicated().any() or len(comparison) != 3:
        raise ValueError("Expected three distinct single-predictor models")
    rejected, adjusted, _, _ = multipletests(
        comparison["P_HC3"].to_numpy(), alpha=0.05, method="fdr_bh"
    )
    comparison["P_FDR_BH_Across_3_Dynamics"] = adjusted
    comparison["Significant_FDR_0_05"] = rejected
    keep = [
        "Predictor",
        "N",
        "Beta_MOCA_per_predictor_SD",
        "HC3_SE",
        "CI95_Lower",
        "CI95_Upper",
        "P_HC3",
        "P_FDR_BH_Across_3_Dynamics",
        "Significant_FDR_0_05",
        "R2_Full_Model",
        "Adjusted_R2_Full_Model",
        "Delta_R2_Predictor",
        "Partial_R2_Predictor",
    ]
    comparison[keep].to_csv(OUTPUT, index=False)
    print("\nMOCA three-predictor comparison with BH-FDR:")
    print(comparison[keep].to_string(index=False))
    print(f"Saved: {OUTPUT.resolve()}")


if __name__ == "__main__":
    main()
