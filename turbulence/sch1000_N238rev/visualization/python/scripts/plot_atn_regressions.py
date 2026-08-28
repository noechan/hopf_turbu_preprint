#!/usr/bin/env python3
"""Plot manuscript AT(N) panels from saved N145 regression outputs.

Statistical models are fitted upstream by
``statistical_analysis/atn_biomarkers/run_atn_regressions_N145.py``. This
script loads aggregate coefficients and prediction lines; it does not fit or
modify any statistical model.
"""

from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
import tempfile
from typing import Final

os.environ.setdefault(
    "MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "turbu_n145_matplotlib")
)
os.environ.setdefault("MPLBACKEND", "Agg")

import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
import pandas as pd


PYTHON_DIR: Final = Path(__file__).resolve().parent.parent
SCH1000_ROOT: Final = PYTHON_DIR.parent.parent
DEFAULT_HARMONIZED_FILE: Final = (
    SCH1000_ROOT
    / "harmonization_allfeat"
    / "recomputed"
    / "Turbu_ComBat_ADNI3_allfeatures_N145.xlsx"
)
DEFAULT_RESULTS_DIR: Final = (
    SCH1000_ROOT
    / "statistical_analysis"
    / "atn_biomarkers"
    / "results"
    / "N145_atn_regressions"
)
DEFAULT_OUTPUT_ROOT: Final = (
    SCH1000_ROOT
    / "figures_N145"
    / "sch1000"
    / "Abeta_Status"
    / "harmonized_allfeat"
    / "python"
    / "atn_regressions"
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

PREDICTOR: Final = "Turbu_lam_0_01"
GROUP_ORDER: Final = ["HC_ABneg", "HC_ABpos", "MCI_ABpos", "AD_ABpos"]
EXPECTED_GROUP_COUNTS: Final = {
    "HC_ABneg": 51,
    "HC_ABpos": 37,
    "MCI_ABpos": 31,
    "AD_ABpos": 26,
}
GROUP_COLORS: Final = {
    "HC_ABneg": "#0173B2",
    "HC_ABpos": "#DE8F05",
    "MCI_ABpos": "#029E73",
    "AD_ABpos": "#D55E00",
}
OUTCOMES: Final = [
    {
        "panel": "Fig4c",
        "column": "CL_pvc",
        "label": "Amyloid burden (CL-PVC)",
        "filename": "Fig4c_Amyloid_CL_vs_Turbulence_lambda_0_01",
    },
    {
        "panel": "Fig4d",
        "column": "Mean_GMV_HIP_BI",
        "label": "Bilateral hippocampal volume",
        "filename": "Fig4d_Hippocampal_GMV_vs_Turbulence_lambda_0_01",
    },
    {
        "panel": "Fig4j",
        "column": "tau_mesial_pvc",
        "label": "Mesial tau (PVC)",
        "filename": "Fig4j_Tau_mesial_vs_Turbulence_lambda_0_01",
    },
    {
        "panel": "Fig4k",
        "column": "tau_metatemporal_pvc",
        "label": "Metatemporal tau (PVC)",
        "filename": "Fig4k_Tau_metatemporal_vs_Turbulence_lambda_0_01",
    },
    {
        "panel": "Fig4l",
        "column": "tau_temporoparietal_pvc",
        "label": "Temporoparietal tau (PVC)",
        "filename": "Fig4l_Tau_temporoparietal_vs_Turbulence_lambda_0_01",
    },
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--model",
        choices=("age_sex", "age_sex_education"),
        default="age_sex_education",
        help="Saved statistical model to visualize (default: age_sex_education)",
    )
    parser.add_argument(
        "--harmonized-file", type=Path, default=DEFAULT_HARMONIZED_FILE
    )
    parser.add_argument("--tau-file", type=Path, default=DEFAULT_TAU_FILE)
    parser.add_argument("--vbm-file", type=Path, default=DEFAULT_VBM_FILE)
    parser.add_argument("--results-dir", type=Path, default=DEFAULT_RESULTS_DIR)
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=None,
        help="Figure destination (default: atn_regressions/<model>)",
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


def load_scatter_inputs(
    harmonized_file: Path,
    tau_file: Path,
    vbm_file: Path,
) -> pd.DataFrame:
    for path, label in (
        (harmonized_file, "harmonized N145 input"),
        (tau_file, "amyloid/demographic/tau input"),
        (vbm_file, "VBM input"),
    ):
        require_file(path, label)
    harmonized = normalize_ptid(pd.read_excel(harmonized_file), "harmonized input")
    tau = normalize_ptid(pd.read_excel(tau_file), "amyloid/demographic/tau input")
    vbm = normalize_ptid(pd.read_csv(vbm_file), "VBM input")
    tau_columns = [
        "PTID",
        "CL_pvc",
        "tau_mesial_pvc",
        "tau_metatemporal_pvc",
        "tau_temporoparietal_pvc",
    ]
    vbm_columns = ["PTID", "Mean_GMV_HIP_LH", "Mean_GMV_HIP_RH"]
    require_columns(harmonized, ["PTID", "Group", PREDICTOR], "harmonized input")
    require_columns(tau, tau_columns, "amyloid/demographic/tau input")
    require_columns(vbm, vbm_columns, "VBM input")
    group_counts = harmonized["Group"].value_counts().to_dict()
    if group_counts != EXPECTED_GROUP_COUNTS:
        raise ValueError(
            f"Expected canonical groups {EXPECTED_GROUP_COUNTS}; found {group_counts}"
        )
    merged = harmonized[["PTID", "Group", PREDICTOR]].merge(
        tau[tau_columns], on="PTID", how="left", validate="one_to_one"
    )
    merged = merged.merge(
        vbm[vbm_columns], on="PTID", how="left", validate="one_to_one"
    )
    for column in [PREDICTOR, *tau_columns[1:], *vbm_columns[1:]]:
        merged[column] = to_numeric(merged[column], column)
    merged["Mean_GMV_HIP_BI"] = merged[
        ["Mean_GMV_HIP_LH", "Mean_GMV_HIP_RH"]
    ].mean(axis=1, skipna=False)
    return merged


def load_statistical_outputs(
    results_dir: Path,
    model_name: str,
) -> tuple[pd.DataFrame, pd.DataFrame, Path, Path]:
    results_path = results_dir / f"N145_ATN_regression_{model_name}.csv"
    predictions_path = results_dir / f"N145_ATN_prediction_lines_{model_name}.csv"
    require_file(results_path, "AT(N) aggregate regression results")
    require_file(predictions_path, "AT(N) saved prediction lines")
    results = pd.read_csv(results_path)
    predictions = pd.read_csv(predictions_path)
    require_columns(
        results,
        [
            "Model",
            "Panel",
            "Outcome",
            "N",
            "Beta_Unstandardized",
            "P_FDR_BH_Across_5_Outcomes",
            "R2_Full_Model",
        ],
        "regression results",
    )
    require_columns(
        predictions,
        ["Model", "Panel", "Outcome", PREDICTOR, "Predicted_Outcome"],
        "prediction lines",
    )
    expected_panels = {item["panel"] for item in OUTCOMES}
    if len(results) != len(OUTCOMES) or set(results["Panel"]) != expected_panels:
        raise ValueError("Regression result panels do not match the five AT(N) outcomes")
    if set(results["Model"]) != {model_name} or set(predictions["Model"]) != {model_name}:
        raise ValueError(f"Saved outputs do not match requested model {model_name}")
    return results, predictions, results_path, predictions_path


def plot_outcome(
    frame: pd.DataFrame,
    outcome: dict[str, str],
    result: pd.Series,
    prediction: pd.DataFrame,
    output_dir: Path,
) -> None:
    column = outcome["column"]
    plotted = frame[["Group", PREDICTOR, column]].dropna().copy()
    if len(plotted) != int(result["N"]):
        raise ValueError(
            f"Plot sample for {outcome['panel']} is N={len(plotted)}, "
            f"but statistical output records N={int(result['N'])}"
        )
    fig, ax = plt.subplots(figsize=(7, 6))
    for group in GROUP_ORDER:
        group_frame = plotted.loc[plotted["Group"] == group]
        ax.scatter(
            group_frame[PREDICTOR],
            group_frame[column],
            color=GROUP_COLORS[group],
            alpha=0.75,
            s=100,
            edgecolor="none",
        )
    line_label = (
        rf'$\beta$ = {result["Beta_Unstandardized"]:.3f}'
        + "\n"
        + rf'$R^2$ = {result["R2_Full_Model"]:.3f}'
    )
    ax.plot(
        prediction[PREDICTOR],
        prediction["Predicted_Outcome"],
        color="#800080",
        linewidth=5,
        label=line_label,
    )
    ax.set_xlabel(r"Turbulence ($\lambda = 0.01$)", fontsize=18)
    ax.set_ylabel(outcome["label"], fontsize=18)
    ax.xaxis.set_major_locator(MaxNLocator(nbins=5))
    ax.yaxis.set_major_locator(MaxNLocator(nbins=5))
    ax.tick_params(axis="both", direction="out", length=8, width=2, labelsize=25)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_linewidth(1.5)
    ax.spines["bottom"].set_linewidth(1.5)
    ax.legend(fontsize=18, handlelength=0.7, frameon=False, loc="best")
    fig.tight_layout()
    for extension in ("png", "pdf"):
        fig.savefig(
            output_dir / f'{outcome["filename"]}.{extension}',
            dpi=300,
            bbox_inches="tight",
        )
    plt.close(fig)


def main() -> None:
    args = parse_args()
    output_dir = (
        args.output_dir.resolve()
        if args.output_dir is not None
        else (DEFAULT_OUTPUT_ROOT / args.model).resolve()
    )
    output_dir.mkdir(parents=True, exist_ok=True)
    frame = load_scatter_inputs(args.harmonized_file, args.tau_file, args.vbm_file)
    results, predictions, results_path, predictions_path = load_statistical_outputs(
        args.results_dir, args.model
    )
    plt.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "axes.grid": False,
        }
    )
    for outcome in OUTCOMES:
        result = results.loc[results["Panel"] == outcome["panel"]].iloc[0]
        prediction = predictions.loc[predictions["Panel"] == outcome["panel"]]
        if len(prediction) != 200:
            raise ValueError(
                f"Expected 200 prediction points for {outcome['panel']}; found {len(prediction)}"
            )
        plot_outcome(frame, outcome, result, prediction, output_dir)

    provenance = {
        "analysis": "Visualization of saved N145 turbulence--AT(N) regressions",
        "model": args.model,
        "statistical_results_input": str(results_path.resolve()),
        "prediction_lines_input": str(predictions_path.resolve()),
        "statistical_models_refitted_during_visualization": False,
        "group_role": "point colour only",
        "output_directory": str(output_dir),
    }
    summary_path = output_dir / "Fig4_ATN_visualization_provenance.json"
    summary_path.write_text(json.dumps(provenance, indent=2) + "\n", encoding="utf-8")
    print(
        f"Saved five PDF/PNG figure pairs for model {args.model} to: {output_dir}"
    )


if __name__ == "__main__":
    main()
