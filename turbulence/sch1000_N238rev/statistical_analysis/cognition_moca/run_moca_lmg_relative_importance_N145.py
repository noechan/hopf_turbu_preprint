#!/usr/bin/env python3
"""Bootstrap LMG relative importance for the three N145 MOCA predictors.

The primary AT(N)-style model adjusts for age, sex, and education and compares
turbulence, information flow, and 1-information transfer at lambda=0.01. LMG
importance is the average incremental R-squared across all six predictor entry
orders. A paired participant bootstrap provides confidence intervals and tests
whether information flow contributes more than each alternative measure.
"""

from __future__ import annotations

import argparse
import itertools
import json
from pathlib import Path
from types import SimpleNamespace
from typing import Final

import numpy as np
import pandas as pd
from statsmodels.stats.multitest import multipletests

import run_moca_information_flow_regression_N145 as single_model


SCRIPT_DIR: Final = Path(__file__).resolve().parent
DEFAULT_OUTPUT_DIR: Final = SCRIPT_DIR / "results" / "N145_MOCA_LMG_relative_importance"
PREDICTOR_KEYS: Final = (
    "turbulence",
    "information-flow",
    "one-minus-information-transfer",
)
PREDICTOR_COLUMNS: Final = {
    "turbulence": "Turbulence",
    "information-flow": "InformationFlow",
    "one-minus-information-transfer": "OneMinusInformationTransfer",
}
PREDICTOR_LABELS: Final = {
    "Turbulence": "Turbulence",
    "InformationFlow": "Information flow",
    "OneMinusInformationTransfer": "1 - Information transfer",
}
COMPARISONS: Final = (
    ("InformationFlow", "Turbulence"),
    ("InformationFlow", "OneMinusInformationTransfer"),
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--harmonized-file",
        type=Path,
        default=single_model.DEFAULT_HARMONIZED_FILE,
    )
    parser.add_argument(
        "--clinical-file",
        type=Path,
        default=single_model.DEFAULT_CLINICAL_FILE,
    )
    parser.add_argument(
        "--metadata-file",
        type=Path,
        default=single_model.DEFAULT_METADATA_FILE,
    )
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument(
        "--n-bootstrap",
        type=int,
        default=10000,
        help="Paired participant bootstrap draws (default: 10000)",
    )
    parser.add_argument("--seed", type=int, default=2025)
    parser.add_argument(
        "--include-group",
        action="store_true",
        help="Sensitivity model additionally adjusted for amyloid-status Group",
    )
    parser.add_argument(
        "--validate-only",
        action="store_true",
        help="Validate and print point estimates without bootstrapping or writing",
    )
    return parser.parse_args()


def load_analysis_frame(args: argparse.Namespace) -> tuple[pd.DataFrame, list[str]]:
    loader_args = SimpleNamespace(
        harmonized_file=args.harmonized_file,
        clinical_file=args.clinical_file,
        metadata_file=args.metadata_file,
    )
    loaded: dict[str, pd.DataFrame] = {}
    missing_ptids: list[str] | None = None
    for key in PREDICTOR_KEYS:
        spec = single_model.PREDICTOR_SPECS[key]
        frame, missing = single_model.load_analysis_data(
            loader_args,
            spec["source"],
            spec["value"],
            spec["z"],
            spec["label"],
            spec["transform"],
        )
        loaded[key] = frame
        if missing_ptids is None:
            missing_ptids = missing
        elif missing != missing_ptids:
            raise ValueError("The three predictor samples have different missing PTIDs")

    base = loaded["information-flow"][
        ["PTID", "Group", "MOCA", "AGE", "SEX_NUM", "EDUCATION"]
    ].copy()
    for key in PREDICTOR_KEYS:
        spec = single_model.PREDICTOR_SPECS[key]
        column = PREDICTOR_COLUMNS[key]
        source = loaded[key][["PTID", spec["value"]]].rename(
            columns={spec["value"]: column}
        )
        base = base.merge(source, on="PTID", how="left", validate="one_to_one")

    required = [
        "PTID",
        "Group",
        "MOCA",
        "AGE",
        "SEX_NUM",
        "EDUCATION",
        *PREDICTOR_COLUMNS.values(),
    ]
    if base[required].isna().any().any():
        raise ValueError("The joined LMG analysis frame contains missing values")
    if len(base) != 144 or base["PTID"].nunique() != 144:
        raise ValueError("Expected 144 unique participants in the LMG analysis")
    return base, missing_ptids or []


def design_components(
    frame: pd.DataFrame,
    include_group: bool,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    y = frame["MOCA"].to_numpy(dtype=float)
    covariates = frame[["AGE", "SEX_NUM", "EDUCATION"]].to_numpy(dtype=float)
    if include_group:
        group_dummies = pd.get_dummies(
            frame["Group"], drop_first=True, dtype=float
        ).to_numpy(dtype=float)
        covariates = np.column_stack([covariates, group_dummies])
    predictors = frame[list(PREDICTOR_LABELS)].to_numpy(dtype=float)
    return y, covariates, predictors


def model_r2(y: np.ndarray, covariates: np.ndarray, predictors: np.ndarray) -> float:
    columns = [np.ones(len(y)), covariates]
    if predictors.shape[1]:
        columns.append(predictors)
    design = np.column_stack(columns)
    coefficients, _, _, _ = np.linalg.lstsq(design, y, rcond=None)
    residuals = y - design @ coefficients
    total = np.sum((y - y.mean()) ** 2)
    return float(1 - np.sum(residuals**2) / total)


def subset_r2(
    y: np.ndarray,
    covariates: np.ndarray,
    predictors: np.ndarray,
) -> dict[int, float]:
    results: dict[int, float] = {}
    for mask in range(1 << predictors.shape[1]):
        indices = [index for index in range(predictors.shape[1]) if mask & (1 << index)]
        selected = predictors[:, indices] if indices else np.empty((len(y), 0))
        results[mask] = model_r2(y, covariates, selected)
    return results


def lmg_from_subset_r2(values: dict[int, float]) -> np.ndarray:
    n_predictors = len(PREDICTOR_LABELS)
    contributions = np.zeros(n_predictors, dtype=float)
    permutations = list(itertools.permutations(range(n_predictors)))
    for order in permutations:
        mask = 0
        previous = values[mask]
        for predictor in order:
            mask |= 1 << predictor
            current = values[mask]
            contributions[predictor] += current - previous
            previous = current
    return contributions / len(permutations)


def point_statistics(
    y: np.ndarray,
    covariates: np.ndarray,
    predictors: np.ndarray,
) -> tuple[dict[int, float], np.ndarray, np.ndarray, np.ndarray]:
    values = subset_r2(y, covariates, predictors)
    lmg = lmg_from_subset_r2(values)
    full_mask = (1 << predictors.shape[1]) - 1
    unique = np.array(
        [values[full_mask] - values[full_mask ^ (1 << index)] for index in range(3)]
    )
    standalone = np.array([values[1 << index] - values[0] for index in range(3)])
    return values, lmg, unique, standalone


def bootstrap_statistics(
    y: np.ndarray,
    covariates: np.ndarray,
    predictors: np.ndarray,
    n_bootstrap: int,
    seed: int,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    rng = np.random.default_rng(seed)
    n = len(y)
    lmg_draws = np.empty((n_bootstrap, 3), dtype=float)
    unique_draws = np.empty((n_bootstrap, 3), dtype=float)
    standalone_draws = np.empty((n_bootstrap, 3), dtype=float)
    progress_every = max(1, n_bootstrap // 10)
    for draw in range(n_bootstrap):
        indices = rng.integers(0, n, size=n)
        _, lmg, unique, standalone = point_statistics(
            y[indices], covariates[indices], predictors[indices]
        )
        lmg_draws[draw] = lmg
        unique_draws[draw] = unique
        standalone_draws[draw] = standalone
        if (draw + 1) % progress_every == 0:
            print(f"  completed {draw + 1}/{n_bootstrap} bootstrap draws", flush=True)
    return lmg_draws, unique_draws, standalone_draws


def percentile_interval(draws: np.ndarray) -> tuple[float, float]:
    lower, upper = np.quantile(draws, [0.025, 0.975])
    return float(lower), float(upper)


def paired_bootstrap_sign_p(draws: np.ndarray) -> float:
    n = len(draws)
    lower_tail = (np.count_nonzero(draws <= 0) + 1) / (n + 1)
    upper_tail = (np.count_nonzero(draws >= 0) + 1) / (n + 1)
    return float(min(1.0, 2 * min(lower_tail, upper_tail)))


def write_results(
    args: argparse.Namespace,
    frame: pd.DataFrame,
    missing_ptids: list[str],
    values: dict[int, float],
    lmg: np.ndarray,
    unique: np.ndarray,
    standalone: np.ndarray,
    lmg_draws: np.ndarray,
    unique_draws: np.ndarray,
    standalone_draws: np.ndarray,
) -> None:
    predictor_names = list(PREDICTOR_LABELS)
    full_mask = (1 << len(predictor_names)) - 1
    total_dynamics = values[full_mask] - values[0]
    importance_rows = []
    for index, predictor in enumerate(predictor_names):
        lmg_low, lmg_high = percentile_interval(lmg_draws[:, index])
        unique_low, unique_high = percentile_interval(unique_draws[:, index])
        importance_rows.append(
            {
                "Predictor": predictor,
                "Predictor_label": PREDICTOR_LABELS[predictor],
                "N": len(frame),
                "LMG_R2": lmg[index],
                "LMG_percent_of_total_dynamics_R2": 100 * lmg[index] / total_dynamics,
                "LMG_bootstrap_SE": lmg_draws[:, index].std(ddof=1),
                "LMG_CI95_low": lmg_low,
                "LMG_CI95_high": lmg_high,
                "Standalone_delta_R2": standalone[index],
                "Unique_delta_R2_given_other_two": unique[index],
                "Unique_CI95_low": unique_low,
                "Unique_CI95_high": unique_high,
                "Total_dynamics_delta_R2": total_dynamics,
                "Full_model_R2": values[full_mask],
                "Covariate_model_R2": values[0],
            }
        )
    importance = pd.DataFrame(importance_rows)

    comparison_rows = []
    for left, right in COMPARISONS:
        left_index = predictor_names.index(left)
        right_index = predictor_names.index(right)
        difference_draws = lmg_draws[:, left_index] - lmg_draws[:, right_index]
        low, high = percentile_interval(difference_draws)
        comparison_rows.append(
            {
                "Comparison": f"{PREDICTOR_LABELS[left]} > {PREDICTOR_LABELS[right]}",
                "LMG_difference": lmg[left_index] - lmg[right_index],
                "Bootstrap_SE": difference_draws.std(ddof=1),
                "CI95_low": low,
                "CI95_high": high,
                "P_bootstrap_sign_two_sided": paired_bootstrap_sign_p(difference_draws),
                "Bootstrap_probability_difference_gt_0": float(
                    np.mean(difference_draws > 0)
                ),
                "Bootstrap_draws": args.n_bootstrap,
            }
        )
    comparisons = pd.DataFrame(comparison_rows)
    rejected, adjusted, _, _ = multipletests(
        comparisons["P_bootstrap_sign_two_sided"].to_numpy(),
        alpha=0.05,
        method="fdr_bh",
    )
    comparisons["P_FDR_BH_across_2_comparisons"] = adjusted
    comparisons["Significant_FDR_0_05"] = rejected

    subset_rows = []
    for mask, r2 in values.items():
        included = [
            PREDICTOR_LABELS[name]
            for index, name in enumerate(predictor_names)
            if mask & (1 << index)
        ]
        subset_rows.append(
            {
                "Mask": mask,
                "Dynamical_predictors": " + ".join(included) if included else "None",
                "R2": r2,
                "Delta_R2_vs_covariates": r2 - values[0],
            }
        )
    subsets = pd.DataFrame(subset_rows)
    bootstrap_frame = pd.DataFrame(
        {
            "Bootstrap_draw": np.arange(1, args.n_bootstrap + 1),
            **{
                f"LMG_{name}": lmg_draws[:, index]
                for index, name in enumerate(predictor_names)
            },
            **{
                f"Unique_{name}": unique_draws[:, index]
                for index, name in enumerate(predictor_names)
            },
            **{
                f"Standalone_{name}": standalone_draws[:, index]
                for index, name in enumerate(predictor_names)
            },
        }
    )
    bootstrap_frame["LMG_difference_InfoFlow_minus_Turbulence"] = (
        bootstrap_frame["LMG_InformationFlow"] - bootstrap_frame["LMG_Turbulence"]
    )
    bootstrap_frame["LMG_difference_InfoFlow_minus_OneMinusTransfer"] = (
        bootstrap_frame["LMG_InformationFlow"]
        - bootstrap_frame["LMG_OneMinusInformationTransfer"]
    )

    suffix = "group_adjusted" if args.include_group else "age_sex_education"
    args.output_dir.mkdir(parents=True, exist_ok=True)
    importance_path = args.output_dir / f"N145_MOCA_LMG_importance_{suffix}.csv"
    comparison_path = args.output_dir / f"N145_MOCA_LMG_pairwise_tests_{suffix}.csv"
    subset_path = args.output_dir / f"N145_MOCA_LMG_subset_R2_{suffix}.csv"
    bootstrap_path = args.output_dir / f"N145_MOCA_LMG_bootstrap_draws_{suffix}.csv"
    provenance_path = args.output_dir / f"N145_MOCA_LMG_provenance_{suffix}.json"
    importance.to_csv(importance_path, index=False)
    comparisons.to_csv(comparison_path, index=False)
    subsets.to_csv(subset_path, index=False)
    bootstrap_frame.to_csv(bootstrap_path, index=False)
    provenance_path.write_text(
        json.dumps(
            {
                "analysis": "N145 MOCA LMG relative importance",
                "outcome": "MOCA",
                "covariates": [
                    "Age",
                    "Sex",
                    "Education",
                    *(["Group"] if args.include_group else []),
                ],
                "predictors": [PREDICTOR_LABELS[name] for name in predictor_names],
                "scale": "lambda=0.01",
                "bootstrap": "paired participant resampling with replacement",
                "bootstrap_draws": args.n_bootstrap,
                "seed": args.seed,
                "confidence_intervals": "nonparametric percentile 95%",
                "pairwise_test": "two-sided paired bootstrap sign probability",
                "multiplicity": "BH-FDR across information-flow-vs-alternative comparisons",
                "harmonized_input": single_model.portable_path(args.harmonized_file),
                "clinical_input": single_model.portable_path(args.clinical_file),
                "metadata_input": single_model.portable_path(args.metadata_file),
                "missing_moca_ptids": missing_ptids,
            },
            indent=2,
        )
        + "\n",
        encoding="utf-8",
    )

    print("\nLMG relative importance:")
    print(importance.to_string(index=False))
    print("\nPaired bootstrap comparisons:")
    print(comparisons.to_string(index=False))
    print(f"\nSaved results under: {args.output_dir.resolve()}")


def main() -> None:
    args = parse_args()
    if args.n_bootstrap < 100:
        raise ValueError("--n-bootstrap must be at least 100")
    frame, missing_ptids = load_analysis_frame(args)
    y, covariates, predictors = design_components(frame, args.include_group)
    values, lmg, unique, standalone = point_statistics(y, covariates, predictors)
    total_dynamics = values[(1 << len(PREDICTOR_LABELS)) - 1] - values[0]
    if not np.isclose(lmg.sum(), total_dynamics, rtol=1e-10, atol=1e-12):
        raise RuntimeError("LMG contributions do not sum to the total dynamics delta R2")
    print(
        f"Validated N={len(frame)} LMG sample; covariates=age, sex, education"
        f"{', Group' if args.include_group else ''}; predictors=3"
    )
    print(
        "Point LMG contributions: "
        + " | ".join(
            f"{PREDICTOR_LABELS[name]}={lmg[index]:.6f}"
            for index, name in enumerate(PREDICTOR_LABELS)
        )
    )
    if args.validate_only:
        print("Validation completed; no bootstrap was run and no files were written.")
        return
    print(f"Running {args.n_bootstrap} paired participant bootstrap draws...")
    lmg_draws, unique_draws, standalone_draws = bootstrap_statistics(
        y, covariates, predictors, args.n_bootstrap, args.seed
    )
    write_results(
        args,
        frame,
        missing_ptids,
        values,
        lmg,
        unique,
        standalone,
        lmg_draws,
        unique_draws,
        standalone_draws,
    )


if __name__ == "__main__":
    main()
