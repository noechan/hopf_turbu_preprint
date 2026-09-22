"""Create ROC, confusion-matrix, and SHAP plots for all constrained Hopf scenarios."""

import os
import re
import sys
from pathlib import Path

import joblib
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.metrics import auc


# Allow VS Code's "Run Python File" action to find the top-level src package.
PATH_REPO = (Path(__file__).parent / ".." / "..").resolve()
if str(PATH_REPO) not in sys.path:
    sys.path.insert(0, str(PATH_REPO))

from src.plotting.confusion_matrix import plot_single_cm


mpl.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial"],
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    "svg.fonttype": "none",
    "axes.labelsize": 16,
    "xtick.labelsize": 16,
    "ytick.labelsize": 16,
    "legend.fontsize": 14,
    "axes.titlesize": 16,
})


# ---------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------
SCENARIOS = (
    "HCneg_vs_HCpos",
    "HCneg_vs_MCIpos",
    "HCneg_vs_ADpos",
    "MCIpos_vs_ADpos",
)

FEATURE_SET = "all_features_turbu_hopf_combat"
CLASSIFIERS = ("LinearSVM", "PolySVM")
N_FEATURES_PLOT = 9  # plots the first N_FEATURES_PLOT + 1 features

CM_CMAP = "Blues"
COLOR_MEAN_SHAP = "tab:blue"
COLOR_SEL_FREQ = "darkblue"

GROUP_COMPARISON_TITLES = {
    "HCneg_vs_HCpos": r"HC- vs HC+",
    "HCneg_vs_MCIpos": r"HC- vs MCI+",
    "HCneg_vs_ADpos": r"HC- vs AD+",
    "MCIpos_vs_ADpos": r"MCI+ vs AD+",
}

GROUP_COMPARISON_LABELS = {
    "HCneg_vs_HCpos": [r"HC-", r"HC+"],
    "HCneg_vs_MCIpos": [r"HC-", r"MCI+"],
    "HCneg_vs_ADpos": [r"HC-", r"AD+"],
    "MCIpos_vs_ADpos": [r"MCI+", r"AD+"],
}

FEATURE_SET_NAMES = {
    "all_features_turbu_hopf_combat": "Turbulence + Hopf",
}

PATH_RESULTS = (
    PATH_REPO
    / "Results"
    / "final_3d_gs_classification_turbu_hopf_sch1000_constrained"
)


def load_scenario_data(classifier, scenario):
    """Load the ROC, confusion-matrix, and SHAP data for one scenario."""
    scenario_results = PATH_RESULTS / classifier / scenario / FEATURE_SET
    if not scenario_results.is_dir():
        raise FileNotFoundError(
            f"Scenario results directory not found: {scenario_results}"
        )

    roc_data_test = joblib.load(scenario_results / "avg_roc_data_test.joblib")

    fold_dirs = [
        directory
        for directory in os.listdir(scenario_results)
        if (scenario_results / directory).is_dir()
        and re.fullmatch(r"fold_\d+", directory)
    ]
    fold_dirs.sort(key=lambda directory: int(directory.split("_")[1]))

    if not fold_dirs:
        raise FileNotFoundError(
            f"No fold directories found under: {scenario_results}"
        )

    confusion_matrices = []
    for directory in fold_dirs:
        cm_path = scenario_results / directory / "test_cm.txt"
        if not cm_path.exists():
            raise FileNotFoundError(f"Expected confusion matrix missing: {cm_path}")
        confusion_matrices.append(np.loadtxt(cm_path))

    avg_cm_test = np.mean(np.asarray(confusion_matrices), axis=0)

    importance = pd.read_csv(scenario_results / "summary_importance.csv")
    importance = importance.iloc[:N_FEATURES_PLOT + 1].copy()
    return roc_data_test, avg_cm_test, importance


def save_roc_plot(classifier, scenario, roc_data_test, output_dir):
    """Create and save the ROC plot for one scenario."""
    fig, ax = plt.subplots(figsize=(5, 5))
    mean_auc = auc(roc_data_test["fpr"], roc_data_test["tpr"])

    ax.plot(
        roc_data_test["fpr"],
        roc_data_test["tpr"],
        label=f"{FEATURE_SET_NAMES[FEATURE_SET]} (AUC = {mean_auc:.3f})",
    )

    tpr_upper = np.minimum(
        roc_data_test["tpr"] + roc_data_test["std_tpr"], 1
    )
    tpr_lower = np.maximum(
        roc_data_test["tpr"] - roc_data_test["std_tpr"], 0
    )
    ax.fill_between(roc_data_test["fpr"], tpr_lower, tpr_upper, alpha=0.2)

    ax.plot([0, 1], [0, 1], linestyle="--", color="grey")
    ax.set_xlabel("False Positive Rate")
    ax.set_ylabel("True Positive Rate")
    ax.set_title(f"ROC Curve ({GROUP_COMPARISON_TITLES[scenario]})")
    ax.legend()
    ax.set_box_aspect(1)

    fig.tight_layout()
    output_stem = (
        output_dir
        / f"classification_roc_{scenario.lower()}_hopf_IC_S_{classifier.lower()}"
    )
    fig.savefig(output_stem.with_suffix(".pdf"), dpi=600)
    fig.savefig(output_stem.with_suffix(".svg"), dpi=600)
    plt.close(fig)


def save_confusion_matrix_plot(classifier, scenario, avg_cm_test, output_dir):
    """Create and save the average confusion-matrix plot for one scenario."""
    fig, ax = plt.subplots(figsize=(5, 5))

    plot_single_cm(
        avg_cm_test,
        GROUP_COMPARISON_LABELS[scenario],
        ax,
        vmin=0,
        vmax=1,
        cmap=CM_CMAP,
        fontsize=18,
    )
    ax.set_title(f"{FEATURE_SET_NAMES[FEATURE_SET]}\nConfusion Matrix")
    ax.set_box_aspect(1)

    fig.tight_layout()
    output_stem = (
        output_dir
        / f"classification_cm_{scenario.lower()}_hopf_IC_S_{classifier.lower()}"
    )
    fig.savefig(output_stem.with_suffix(".pdf"), dpi=600)
    fig.savefig(output_stem.with_suffix(".svg"), dpi=600)
    plt.close(fig)


def save_shap_plot(classifier, scenario, importance, output_dir):
    """Create and save the SHAP and selection-frequency plot."""
    features = importance["feature"].astype(str).tolist()
    mean_abs_shap = importance["mean_abs_SHAP"].tolist()
    selection_frequency = importance["selection_freq"].tolist()
    x_values = np.arange(len(features))

    single_bar_width = 0.8 / 2
    fig, ax_bar = plt.subplots(figsize=(9.5, 4.5))
    ax_bar_twin = ax_bar.twinx()

    ax_bar.bar(
        x_values - single_bar_width / 2,
        mean_abs_shap,
        width=single_bar_width,
        label="Mean SHAP value",
        color=COLOR_MEAN_SHAP,
    )
    ax_bar_twin.bar(
        x_values + single_bar_width / 2,
        selection_frequency,
        width=single_bar_width,
        label="Selection Frequency",
        color=COLOR_SEL_FREQ,
    )

    ax_bar_twin.set_ylim(0, 1)
    ax_bar.set_xticks(x_values)
    ax_bar.set_xticklabels(features)
    plt.setp(ax_bar.get_xticklabels(), rotation=45, ha="right")

    ax_bar.set_title(f"Feature Importance in {FEATURE_SET_NAMES[FEATURE_SET]}")
    ax_bar.set_ylabel("Mean SHAP value")
    ax_bar_twin.set_ylabel("Feature Selection Frequency")

    handles1, labels1 = ax_bar.get_legend_handles_labels()
    handles2, labels2 = ax_bar_twin.get_legend_handles_labels()
    ax_bar.legend(
        handles1 + handles2,
        labels1 + labels2,
        loc="upper left",
        bbox_to_anchor=(1.18, 1.0),
        borderaxespad=0,
        frameon=False,
    )

    ax_bar.spines["left"].set_color(COLOR_MEAN_SHAP)
    ax_bar.tick_params(axis="y", colors=COLOR_MEAN_SHAP)
    ax_bar_twin.spines["right"].set_color(COLOR_SEL_FREQ)
    ax_bar_twin.tick_params(axis="y", colors=COLOR_SEL_FREQ)
    ax_bar_twin.set_xlim(
        x_values[0] - 1.25 * single_bar_width,
        x_values[-1] + 1.25 * single_bar_width,
    )

    # Reserve space for the external legend and the right-hand twin-axis label.
    fig.subplots_adjust(left=0.10, right=0.72, bottom=0.30, top=0.88)
    output_stem = (
        output_dir
        / f"classification_shap_{scenario.lower()}_hopf_IC_S_{classifier.lower()}"
    )
    fig.savefig(output_stem.with_suffix(".pdf"), dpi=600, bbox_inches="tight")
    fig.savefig(output_stem.with_suffix(".svg"), dpi=600, bbox_inches="tight")
    fig.savefig(output_stem.with_suffix(".png"), dpi=300, bbox_inches="tight")
    plt.close(fig)


def generate_plots(classifier, scenario, output_dir):
    """Generate all plots for one classification scenario."""
    print(f"Generating {classifier} plots for {scenario}...")
    roc_data_test, avg_cm_test, importance = load_scenario_data(
        classifier, scenario
    )
    save_roc_plot(classifier, scenario, roc_data_test, output_dir)
    save_confusion_matrix_plot(classifier, scenario, avg_cm_test, output_dir)
    save_shap_plot(classifier, scenario, importance, output_dir)


def main():
    for classifier in CLASSIFIERS:
        output_dir = PATH_RESULTS / "figures" / classifier
        output_dir.mkdir(parents=True, exist_ok=True)
        for scenario in SCENARIOS:
            generate_plots(classifier, scenario, output_dir)
        print(f"{classifier} plots saved under: {output_dir}")


if __name__ == "__main__":
    main()
