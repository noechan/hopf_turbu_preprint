#!/usr/bin/env python3
"""Plot saved education-adjusted subject-level Hopf--AT(N) regressions."""

from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
import sys
import tempfile
from typing import Final

os.environ.setdefault(
    "MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "turbu_n145_matplotlib")
)
os.environ.setdefault("MPLBACKEND", "Agg")

import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
import pandas as pd


SCRIPT_DIR: Final = Path(__file__).resolve().parent
SCH1000_ROOT: Final = SCRIPT_DIR.parents[2]
ANALYSIS_DIR: Final = (
    SCH1000_ROOT / "statistical_analysis" / "perturbation_measures"
)
sys.path.insert(0, str(ANALYSIS_DIR))

import run_atn_subjectlevel_perturbation_N145 as atn  # noqa: E402


DEFAULT_RESULTS_DIR: Final = (
    ANALYSIS_DIR / "results" / "N145_ATN_subjectlevel_perturbation"
)
DEFAULT_OUTPUT_DIR: Final = (
    SCH1000_ROOT
    / "figures_N145"
    / "sch1000"
    / "Abeta_Status"
    / "harmonized_allfeat"
    / "python"
    / "atn_regressions"
    / "subjectlevel_age_sex_education"
)
RESULT_FILENAME: Final = (
    "N145_ATN_subjectlevel_perturbation_HC3_age_sex_education.csv"
)
PREDICTION_FILENAME: Final = (
    "N145_ATN_subjectlevel_perturbation_prediction_lines_age_sex_education.csv"
)

GROUP_LABELS: Final = {
    "HC_ABneg": r"HC$^-$",
    "HC_ABpos": r"HC$^+$",
    "MCI_ABpos": r"MCI$^+$",
    "AD_ABpos": r"AD$^+$",
}
GROUP_COLORS: Final = {
    "HC_ABneg": "#0173B2",
    "HC_ABpos": "#E69F00",
    "MCI_ABpos": "#029E73",
    "AD_ABpos": "#D55E00",
}
PREDICTOR_SPECS: Final = {
    "information-capability": {
        "result_label": "InformationCapacity",
        "column": "Info_Cap",
        "display": "Information capability",
        "file_display": "InformationCapability",
    },
    "susceptibility": {
        "result_label": "Susceptibility",
        "column": "Susceptibility",
        "display": "Susceptibility",
        "file_display": "Susceptibility",
    },
}
OUTCOME_SPECS: Final = {
    "CL_pvc": {
        "display": "Amyloid burden (CL-PVC)",
        "file_display": "Amyloid_CL",
    },
    "Mean_GMV_HIP_BI": {
        "display": "Bilateral hippocampal volume",
        "file_display": "Hippocampal_GMV",
    },
    "tau_mesial_pvc": {
        "display": "Mesial tau burden (SUVR)",
        "file_display": "Tau_mesial",
    },
    "tau_metatemporal_pvc": {
        "display": "Metatemporal tau burden (SUVR)",
        "file_display": "Tau_metatemporal",
    },
    "tau_temporoparietal_pvc": {
        "display": "Temporoparietal tau burden (SUVR)",
        "file_display": "Tau_temporoparietal",
    },
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--predictor",
        choices=("both", *PREDICTOR_SPECS),
        default="both",
        help="Hopf measure to plot (default: both)",
    )
    parser.add_argument("--hopf-file", type=Path, default=atn.DEFAULT_HOPF_FILE)
    parser.add_argument("--tau-file", type=Path, default=atn.DEFAULT_TAU_FILE)
    parser.add_argument("--vbm-file", type=Path, default=atn.DEFAULT_VBM_FILE)
    parser.add_argument(
        "--metadata-file", type=Path, default=atn.DEFAULT_METADATA_FILE
    )
    parser.add_argument("--results-dir", type=Path, default=DEFAULT_RESULTS_DIR)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    return parser.parse_args()


def load_saved_outputs(
    results_dir: Path,
) -> tuple[pd.DataFrame, pd.DataFrame, Path, Path]:
    result_path = results_dir / RESULT_FILENAME
    prediction_path = results_dir / PREDICTION_FILENAME
    atn.require_file(result_path, "education-adjusted AT(N) result")
    atn.require_file(prediction_path, "education-adjusted AT(N) prediction line")
    results = pd.read_csv(result_path)
    predictions = pd.read_csv(prediction_path)
    atn.require_columns(
        results,
        [
            "Predictor_Label",
            "Predictor",
            "Formula",
            "Outcome",
            "N",
            "Beta_Per_Raw_Predictor_Unit",
            "R2_Full_Model",
            "P_FDR_BH_Across_5_Outcomes",
            "Significant_FDR_0_05",
        ],
        "AT(N) result",
    )
    atn.require_columns(
        predictions,
        [
            "Predictor_Label",
            "Predictor",
            "Outcome",
            "Predictor_Raw_Value",
            "Predicted_Outcome",
        ],
        "AT(N) prediction line",
    )
    if len(results) != len(atn.PREDICTORS) * len(atn.OUTCOMES):
        raise ValueError(f"Expected 10 saved AT(N) models; found {len(results)}")
    if not results["Formula"].str.contains("Education", regex=False).all():
        raise ValueError("Saved AT(N) results are not all adjusted for education")
    return results, predictions, result_path, prediction_path


def plot_model(
    frame: pd.DataFrame,
    predictor_spec: dict[str, str],
    outcome: dict[str, str],
    result: pd.Series,
    prediction: pd.DataFrame,
    output_dir: Path,
) -> list[Path]:
    predictor = predictor_spec["column"]
    outcome_column = outcome["column"]
    outcome_spec = OUTCOME_SPECS[outcome_column]
    plotted = frame[["Group", predictor, outcome_column]].dropna().copy()
    if len(plotted) != int(result["N"]):
        raise ValueError(
            f"{predictor}/{outcome_column}: plot N={len(plotted)} but result N={int(result['N'])}"
        )
    if len(prediction) != 200:
        raise ValueError(
            f"{predictor}/{outcome_column}: expected 200 prediction points; "
            f"found {len(prediction)}"
        )
    prediction = prediction.sort_values("Predictor_Raw_Value")

    fig, ax = plt.subplots(figsize=(7.0, 6.0))
    for group in atn.GROUP_ORDER:
        group_frame = plotted.loc[plotted["Group"] == group]
        ax.scatter(
            group_frame[predictor],
            group_frame[outcome_column],
            color=GROUP_COLORS[group],
            alpha=0.72,
            s=78,
            edgecolor="white",
            linewidth=0.45,
        )

    significant = bool(result["Significant_FDR_0_05"])
    beta_marker = r"\beta^{*}" if significant else r"\beta"
    line_label = (
        rf"${beta_marker}={result['Beta_Per_Raw_Predictor_Unit']:.3f}$"
        "\n"
        rf"$R^2={result['R2_Full_Model']:.3f}$"
    )
    ax.plot(
        prediction["Predictor_Raw_Value"],
        prediction["Predicted_Outcome"],
        color="#800080",
        linewidth=4.0,
        label=line_label,
        zorder=5,
    )
    ax.set_xlabel(predictor_spec["display"], fontsize=17)
    ax.set_ylabel(outcome_spec["display"], fontsize=17)
    ax.xaxis.set_major_locator(MaxNLocator(nbins=5))
    ax.yaxis.set_major_locator(MaxNLocator(nbins=5))
    ax.tick_params(axis="both", direction="out", length=7, width=1.5, labelsize=15)
    ax.spines[["top", "right"]].set_visible(False)
    ax.spines[["left", "bottom"]].set_linewidth(1.4)
    ax.legend(frameon=False, fontsize=15, handlelength=0.8, loc="best")
    fig.tight_layout()

    stem = (
        f"{outcome_spec['file_display']}_vs_{predictor_spec['file_display']}_"
        "age_sex_education"
    )
    outputs = [output_dir / f"{stem}.pdf", output_dir / f"{stem}.png"]
    fig.savefig(outputs[0], bbox_inches="tight")
    fig.savefig(outputs[1], dpi=300, bbox_inches="tight")
    plt.close(fig)
    return outputs


def main() -> None:
    args = parse_args()
    frame, predictor_correlation = atn.load_analysis_data(args)
    results, predictions, result_path, prediction_path = load_saved_outputs(
        args.results_dir
    )
    args.output_dir.mkdir(parents=True, exist_ok=True)
    plt.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "axes.grid": False,
        }
    )

    selected = (
        PREDICTOR_SPECS
        if args.predictor == "both"
        else {args.predictor: PREDICTOR_SPECS[args.predictor]}
    )
    output_files: list[Path] = []
    for predictor_spec in selected.values():
        result_label = predictor_spec["result_label"]
        for outcome in atn.OUTCOMES:
            outcome_column = outcome["column"]
            matches = results.loc[
                (results["Predictor_Label"] == result_label)
                & (results["Outcome"] == outcome_column)
            ]
            if len(matches) != 1:
                raise ValueError(
                    f"Expected one result for {result_label}/{outcome_column}; "
                    f"found {len(matches)}"
                )
            prediction = predictions.loc[
                (predictions["Predictor_Label"] == result_label)
                & (predictions["Outcome"] == outcome_column)
            ]
            output_files.extend(
                plot_model(
                    frame,
                    predictor_spec,
                    outcome,
                    matches.iloc[0],
                    prediction,
                    args.output_dir,
                )
            )

    provenance = {
        "analysis": "Visualization of saved N145 subject-level Hopf--AT(N) regressions",
        "statistical_results_input": str(result_path.resolve()),
        "prediction_lines_input": str(prediction_path.resolve()),
        "covariates": ["age", "sex", "education"],
        "statistical_models_refitted_during_visualization": False,
        "group_role": "scatter-point colour only; group is not a model covariate",
        "predictor_correlation": predictor_correlation,
        "star_definition": "BH-FDR p < 0.05 across five outcomes within predictor",
        "outputs": [str(path.resolve()) for path in output_files],
    }
    provenance_path = args.output_dir / "subjectlevel_ATN_visualization_provenance.json"
    provenance_path.write_text(
        json.dumps(provenance, indent=2) + "\n", encoding="utf-8"
    )
    print(
        f"Saved {len(output_files) // 2} PDF/PNG figure pairs to: "
        f"{args.output_dir.resolve()}"
    )


if __name__ == "__main__":
    main()
