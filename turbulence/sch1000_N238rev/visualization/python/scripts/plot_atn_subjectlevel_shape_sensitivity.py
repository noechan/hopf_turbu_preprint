#!/usr/bin/env python3
"""Plot Figure 4 group-specific slopes and pooled quadratic sensitivity fits."""

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
from matplotlib.lines import Line2D
from matplotlib.ticker import MaxNLocator
import pandas as pd


SCRIPT_DIR: Final = Path(__file__).resolve().parent
SCH1000_ROOT: Final = SCRIPT_DIR.parents[2]
ANALYSIS_DIR: Final = SCH1000_ROOT / "statistical_analysis" / "perturbation_measures"
sys.path.insert(0, str(ANALYSIS_DIR))

import run_atn_subjectlevel_perturbation_N145 as atn  # noqa: E402
import run_atn_subjectlevel_shape_sensitivity_N145 as sensitivity  # noqa: E402


DEFAULT_RESULTS_DIR: Final = sensitivity.DEFAULT_OUTPUT_DIR
DEFAULT_OUTPUT_DIR: Final = (
    SCH1000_ROOT
    / "figures_N145"
    / "sch1000"
    / "Abeta_Status"
    / "harmonized_allfeat"
    / "python"
    / "atn_regressions"
    / "subjectlevel_shape_sensitivity"
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


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--hopf-file", type=Path, default=atn.DEFAULT_HOPF_FILE)
    parser.add_argument("--tau-file", type=Path, default=atn.DEFAULT_TAU_FILE)
    parser.add_argument("--vbm-file", type=Path, default=atn.DEFAULT_VBM_FILE)
    parser.add_argument("--metadata-file", type=Path, default=atn.DEFAULT_METADATA_FILE)
    parser.add_argument("--results-dir", type=Path, default=DEFAULT_RESULTS_DIR)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    return parser.parse_args()


def require_saved_results(results_dir: Path) -> dict[str, pd.DataFrame]:
    files = {
        "group": "N145_Figure4_information_capability_group_slopes_HC3.csv",
        "group_predictions": "N145_Figure4_information_capability_group_prediction_lines.csv",
        "interactions": "N145_Figure4_information_capability_group_interactions_HC3.csv",
        "shape": "N145_Figure4_information_capability_pooled_quadratic_HC3.csv",
        "shape_predictions": "N145_Figure4_information_capability_pooled_shape_prediction_lines.csv",
    }
    tables: dict[str, pd.DataFrame] = {}
    for name, filename in files.items():
        path = results_dir / filename
        atn.require_file(path, f"shape-sensitivity {name} result")
        tables[name] = pd.read_csv(path)
    return tables


def scatter_groups(axis: plt.Axes, frame: pd.DataFrame, outcome: str) -> None:
    for group in atn.GROUP_ORDER:
        subset = frame.loc[frame["Group"] == group, [sensitivity.PREDICTOR, outcome]].dropna()
        axis.scatter(
            subset[sensitivity.PREDICTOR],
            subset[outcome],
            color=GROUP_COLORS[group],
            alpha=0.70,
            s=34,
            edgecolor="white",
            linewidth=0.35,
            zorder=2,
        )


def format_axis(axis: plt.Axes, outcome_label: str, show_xlabel: bool = True) -> None:
    if show_xlabel:
        axis.set_xlabel("Information capability")
    axis.set_ylabel(outcome_label)
    axis.xaxis.set_major_locator(MaxNLocator(nbins=5))
    axis.yaxis.set_major_locator(MaxNLocator(nbins=5))
    axis.spines[["top", "right"]].set_visible(False)
    axis.grid(False)


def plot_group_panel(
    axis: plt.Axes,
    frame: pd.DataFrame,
    outcome: dict[str, str],
    results: pd.DataFrame,
    predictions: pd.DataFrame,
    title: str | None = None,
) -> None:
    outcome_column = outcome["column"]
    scatter_groups(axis, frame, outcome_column)
    handles: list[Line2D] = []
    for group in atn.GROUP_ORDER:
        row = results.loc[
            (results["Outcome"] == outcome_column) & (results["Group"] == group)
        ]
        line = predictions.loc[
            (predictions["Outcome"] == outcome_column)
            & (predictions["Group"] == group)
        ].sort_values("Predictor_Raw_Value")
        if len(row) != 1 or len(line) != 200:
            raise ValueError(f"Unexpected saved group result for {group}/{outcome_column}")
        row = row.iloc[0]
        axis.plot(
            line["Predictor_Raw_Value"],
            line["Predicted_Outcome"],
            color=GROUP_COLORS[group],
            linewidth=2.0,
            zorder=4,
        )
        handles.append(
            Line2D(
                [0],
                [0],
                color=GROUP_COLORS[group],
                linewidth=2.0,
                label=(
                    rf"{GROUP_LABELS[group]}: "
                    + r"$\beta="
                    + f"{row['Beta_Per_Pooled_SD']:.2f}"
                    + r"$, $p="
                    + f"{row['P_Raw']:.3f}"
                    + "$"
                ),
            )
        )
    format_axis(axis, outcome["label"])
    if title:
        axis.set_title(title, fontweight="bold")
    axis.legend(handles=handles, frameon=False, fontsize=7.2, loc="best")


def plot_shape_panel(
    axis: plt.Axes,
    frame: pd.DataFrame,
    outcome: dict[str, str],
    result: pd.DataFrame,
    predictions: pd.DataFrame,
    title: str | None = None,
) -> None:
    outcome_column = outcome["column"]
    scatter_groups(axis, frame, outcome_column)
    row = result.loc[result["Outcome"] == outcome_column]
    if len(row) != 1:
        raise ValueError(f"Expected one pooled shape result for {outcome_column}")
    row = row.iloc[0]
    for model_name, color, style, width, label in (
        ("linear", "#555555", "--", 1.8, "Pooled linear"),
        ("quadratic", "#800080", "-", 2.5, "Pooled quadratic"),
    ):
        line = predictions.loc[
            (predictions["Outcome"] == outcome_column)
            & (predictions["Model"] == model_name)
        ].sort_values("Predictor_Raw_Value")
        if len(line) != 250:
            raise ValueError(f"Unexpected saved {model_name} prediction for {outcome_column}")
        axis.plot(
            line["Predictor_Raw_Value"],
            line["Predicted_Outcome"],
            color=color,
            linestyle=style,
            linewidth=width,
            label=label,
            zorder=4,
        )
    format_axis(axis, outcome["label"])
    if title:
        axis.set_title(title, fontweight="bold")
    axis.legend(frameon=False, fontsize=8, loc="best")
    axis.text(
        0.02,
        0.02,
        (
            rf"Quadratic $p={row['Quadratic_HC3_P_Raw']:.3f}$; "
            rf"$p_{{FDR}}={row['Quadratic_HC3_P_FDR_BH_Across_4_Outcomes']:.3f}$"
            "\n"
            rf"$Delta$AIC={row['Delta_AIC_Quadratic_minus_Linear']:.2f}"
        ),
        transform=axis.transAxes,
        ha="left",
        va="bottom",
        fontsize=7.5,
    )


def save_figure(figure: plt.Figure, stem: Path) -> list[Path]:
    outputs = [stem.with_suffix(".pdf"), stem.with_suffix(".png")]
    figure.savefig(outputs[0], bbox_inches="tight")
    figure.savefig(outputs[1], dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(figure)
    return outputs


def main() -> None:
    args = parse_args()
    frame, _ = atn.load_analysis_data(args)
    tables = require_saved_results(args.results_dir)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    plt.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
            "font.size": 9,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "axes.grid": False,
        }
    )

    outputs: list[Path] = []
    for outcome in sensitivity.FIGURE4_OUTCOMES:
        group_figure, group_axis = plt.subplots(figsize=(6.4, 5.4), facecolor="white")
        plot_group_panel(
            group_axis,
            frame,
            outcome,
            tables["group"],
            tables["group_predictions"],
            title=f"Within-group linear fits: {outcome['label']}",
        )
        group_figure.tight_layout()
        outputs.extend(
            save_figure(
                group_figure,
                args.output_dir / f"{outcome['panel']}_within_group_linear",
            )
        )

        shape_figure, shape_axis = plt.subplots(figsize=(6.4, 5.4), facecolor="white")
        plot_shape_panel(
            shape_axis,
            frame,
            outcome,
            tables["shape"],
            tables["shape_predictions"],
            title=f"Pooled linear and quadratic fits: {outcome['label']}",
        )
        shape_figure.tight_layout()
        outputs.extend(
            save_figure(
                shape_figure,
                args.output_dir / f"{outcome['panel']}_pooled_linear_quadratic",
            )
        )

    overview, axes = plt.subplots(2, 4, figsize=(18.0, 8.6), facecolor="white")
    for column, outcome in enumerate(sensitivity.FIGURE4_OUTCOMES):
        plot_group_panel(
            axes[0, column],
            frame,
            outcome,
            tables["group"],
            tables["group_predictions"],
            title=outcome["label"],
        )
        plot_shape_panel(
            axes[1, column],
            frame,
            outcome,
            tables["shape"],
            tables["shape_predictions"],
        )
        if column:
            axes[0, column].set_ylabel("")
            axes[1, column].set_ylabel("")
    axes[0, 0].text(
        -0.22,
        1.12,
        "Within-group linear fits",
        transform=axes[0, 0].transAxes,
        fontsize=12,
        fontweight="bold",
    )
    axes[1, 0].text(
        -0.22,
        1.06,
        "Pooled shape comparison",
        transform=axes[1, 0].transAxes,
        fontsize=12,
        fontweight="bold",
    )
    overview.tight_layout(w_pad=2.0, h_pad=2.3)
    outputs.extend(save_figure(overview, args.output_dir / "Figure4_shape_sensitivity_overview"))

    provenance = {
        "analysis": "Figure 4 shape-sensitivity visualization",
        "participant_level_data_saved": False,
        "statistical_results_refitted_during_visualization": False,
        "results_dir": str(args.results_dir.resolve()),
        "output_files": [str(path.resolve()) for path in outputs],
        "group_line_definition": (
            "Separate HC3 models adjusted for age, sex, and education; lines "
            "evaluated at common sample-wide reference covariates"
        ),
        "pooled_shape_definition": (
            "HC3 pooled linear and quadratic models adjusted for age, sex, and education"
        ),
    }
    provenance_path = args.output_dir / "Figure4_shape_sensitivity_visualization_provenance.json"
    provenance_path.write_text(json.dumps(provenance, indent=2) + "\n", encoding="utf-8")
    print(f"Saved {len(outputs) // 2} PDF/PNG figure pairs to: {args.output_dir.resolve()}")


if __name__ == "__main__":
    main()
