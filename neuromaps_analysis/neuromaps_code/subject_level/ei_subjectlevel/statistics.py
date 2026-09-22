"""Participant-level regression, permutation inference, and descriptives."""

from __future__ import annotations

import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests


def _hc3_fit(y: np.ndarray, design: np.ndarray) -> dict[str, np.ndarray | int]:
    """Fit ordinary least squares with HC3 heteroscedasticity-robust errors."""

    xtx = np.einsum("ni,nj->ij", design, design, optimize=True)
    xtx_inv = np.linalg.pinv(xtx)
    beta = np.linalg.lstsq(design, y, rcond=None)[0]
    residual = y - np.einsum("ij,j->i", design, beta, optimize=True)
    design_xtx_inv = np.einsum("ij,jk->ik", design, xtx_inv, optimize=True)
    leverage = np.einsum("ij,ij->i", design_xtx_inv, design, optimize=True)
    adjusted = residual / np.clip(1 - leverage, 1e-8, None)
    meat = np.einsum("ni,n,nj->ij", design, adjusted**2, design, optimize=True)
    covariance = np.einsum(
        "ij,jk,kl->il", xtx_inv, meat, xtx_inv, optimize=True
    )
    standard_error = np.sqrt(np.clip(np.diag(covariance), 0, None))
    df = design.shape[0] - np.linalg.matrix_rank(design)
    t_values = beta / standard_error
    p_values = 2 * stats.t.sf(np.abs(t_values), df)
    return {
        "beta": beta,
        "se": standard_error,
        "t": t_values,
        "p": p_values,
        "df": int(df),
    }


def _freedman_lane_coefficient_test(
    y: np.ndarray,
    full_design: np.ndarray,
    tested_column: int,
    n_perm: int,
    rng: np.random.Generator,
) -> tuple[float, float]:
    """Test one regression coefficient by permuting reduced-model residuals."""

    reduced = np.delete(full_design, tested_column, axis=1)
    reduced_beta = np.linalg.pinv(reduced) @ y
    fitted = reduced @ reduced_beta
    residual = y - fitted
    full_pinv = np.linalg.pinv(full_design)
    observed = float((full_pinv @ y)[tested_column])
    permuted = np.empty(n_perm, dtype=float)
    for permutation in range(n_perm):
        permuted_y = fitted + residual[rng.permutation(len(y))]
        permuted[permutation] = (full_pinv @ permuted_y)[tested_column]
    two_sided = (np.sum(np.abs(permuted) >= abs(observed)) + 1) / (n_perm + 1)
    lower_tailed = (np.sum(permuted <= observed) + 1) / (n_perm + 1)
    return float(two_sided), float(lower_tailed)


def run_models(
    participant_table: pd.DataFrame,
    outcome: str,
    group_order: list[str],
    n_perm: int,
    seed: int,
) -> tuple[list[dict[str, object]], dict[str, object]]:
    """Run the ordered-stage model and three planned HC-negative contrasts."""

    y = participant_table[outcome].to_numpy(dtype=float)
    age = stats.zscore(participant_table["age"].to_numpy(dtype=float), ddof=0)
    edu = stats.zscore(participant_table["edu"].to_numpy(dtype=float), ddof=0)
    gender = participant_table["gender"].to_numpy(dtype=float)
    stage = participant_table["stage"].to_numpy(dtype=float)
    rng = np.random.default_rng(seed)

    # Confirmatory model: does coupling decrease monotonically from stage 0 to 3?
    trend_design = np.column_stack([np.ones(len(y)), stage, age, gender, edu])
    trend_fit = _hc3_fit(y, trend_design)
    trend_perm_two, trend_perm_lower = _freedman_lane_coefficient_test(
        y, trend_design, tested_column=1, n_perm=n_perm, rng=rng
    )
    critical = stats.t.ppf(0.975, trend_fit["df"])
    trend = {
        "outcome": outcome,
        "test": "ordered_linear_trend",
        "coefficient": float(trend_fit["beta"][1]),
        "standard_error_hc3": float(trend_fit["se"][1]),
        "ci95_low_hc3": float(trend_fit["beta"][1] - critical * trend_fit["se"][1]),
        "ci95_high_hc3": float(trend_fit["beta"][1] + critical * trend_fit["se"][1]),
        "t_hc3": float(trend_fit["t"][1]),
        "df_residual": int(trend_fit["df"]),
        "p_hc3_two_sided": float(trend_fit["p"][1]),
        "p_permutation_two_sided": trend_perm_two,
        "p_permutation_directional_decrease": trend_perm_lower,
        "n_permutations": n_perm,
    }

    # Planned categorical models: compare each later group with HC amyloid-negative.
    group = pd.Categorical(participant_table["Group"], categories=group_order)
    dummies = pd.get_dummies(group, dtype=float).drop(columns=group_order[0])
    categorical_design = np.column_stack(
        [np.ones(len(y)), dummies.to_numpy(), age, gender, edu]
    )
    categorical_fit = _hc3_fit(y, categorical_design)
    comparisons: list[dict[str, object]] = []
    for index, comparison_group in enumerate(group_order[1:], start=1):
        p_perm_two, p_perm_lower = _freedman_lane_coefficient_test(
            y,
            categorical_design,
            tested_column=index,
            n_perm=n_perm,
            rng=rng,
        )
        comparisons.append(
            {
                "outcome": outcome,
                "test": "planned_group_comparison",
                "reference_group": group_order[0],
                "comparison_group": comparison_group,
                "coefficient": float(categorical_fit["beta"][index]),
                "standard_error_hc3": float(categorical_fit["se"][index]),
                "ci95_low_hc3": float(
                    categorical_fit["beta"][index]
                    - critical * categorical_fit["se"][index]
                ),
                "ci95_high_hc3": float(
                    categorical_fit["beta"][index]
                    + critical * categorical_fit["se"][index]
                ),
                "t_hc3": float(categorical_fit["t"][index]),
                "df_residual": int(categorical_fit["df"]),
                "p_hc3_two_sided": float(categorical_fit["p"][index]),
                "p_permutation_two_sided": p_perm_two,
                "p_permutation_directional_decrease": p_perm_lower,
                "n_permutations": n_perm,
            }
        )

    rejected, adjusted, _, _ = multipletests(
        [row["p_permutation_two_sided"] for row in comparisons],
        alpha=0.05,
        method="fdr_bh",
    )
    for row, reject, p_fdr in zip(comparisons, rejected, adjusted):
        row["p_permutation_fdr_bh"] = float(p_fdr)
        row["significant_fdr_0.05"] = bool(reject)
    return comparisons, trend


def run_all_models(
    participant_table: pd.DataFrame,
    outcomes: list[str],
    group_order: list[str],
    n_perm: int,
    seed: int,
) -> tuple[list[dict[str, object]], list[dict[str, object]]]:
    """Run identical inference for each configured coupling outcome."""

    all_comparisons: list[dict[str, object]] = []
    trends: list[dict[str, object]] = []
    for outcome_index, outcome in enumerate(outcomes):
        comparisons, trend = run_models(
            participant_table,
            outcome,
            group_order,
            n_perm,
            seed=seed + outcome_index * 100000,
        )
        all_comparisons.extend(comparisons)
        trends.append(trend)
    return all_comparisons, trends


def group_descriptives(
    participant_table: pd.DataFrame,
    group_order: list[str],
    display_names: dict[str, str],
    seed: int,
) -> pd.DataFrame:
    """Summarize participant Pearson correlations within each group."""

    rng = np.random.default_rng(seed)
    rows: list[dict[str, object]] = []
    for group in group_order:
        subset = participant_table.loc[participant_table["Group"] == group]
        r_values = subset["pearson_r"].to_numpy(dtype=float)
        bootstrap = np.mean(
            rng.choice(r_values, size=(10000, r_values.size), replace=True), axis=1
        )
        rows.append(
            {
                "group": group,
                "display_name": display_names[group],
                "sample_size": len(subset),
                "mean_r": float(r_values.mean()),
                "sd_r": float(r_values.std(ddof=1)),
                "median_r": float(np.median(r_values)),
                "q1_r": float(np.quantile(r_values, 0.25)),
                "q3_r": float(np.quantile(r_values, 0.75)),
                "bootstrap_mean_r_ci95_low": float(np.quantile(bootstrap, 0.025)),
                "bootstrap_mean_r_ci95_high": float(np.quantile(bootstrap, 0.975)),
                "mean_fisher_z": float(subset["fisher_z"].mean()),
                "mean_spin_z": float(subset["spin_z"].mean()),
            }
        )
    return pd.DataFrame(rows)
