from pathlib import Path
import sys

import joblib
import matplotlib as mpl
import matplotlib.pyplot as plt


# Allow the script to be launched from any working directory.
PATH_REPO = Path(__file__).resolve().parents[2]
if str(PATH_REPO) not in sys.path:
    sys.path.insert(0, str(PATH_REPO))

# ---------------------------------------------------------------------
# Matplotlib defaults
# ---------------------------------------------------------------------
mpl.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial"],
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    "svg.fonttype": "none",
    "axes.labelsize": 15,
    "xtick.labelsize": 13,
    "ytick.labelsize": 13,
    "legend.fontsize": 12,
    "axes.titlesize": 15,
    "axes.spines.top": False,
    "axes.spines.right": False,
})

# ---------------------------------------------------------------------
# Labels and configuration
# ---------------------------------------------------------------------
SCENARIOS = [
    "HCneg_vs_ADpos",
    "HCneg_vs_MCIpos",
    "HCneg_vs_HCpos",
    "MCIpos_vs_ADpos",
]

SCENARIO_TITLES = {
    "HCneg_vs_HCpos": r"HC$^-$ vs HC$^+$",
    "HCneg_vs_MCIpos": r"HC$^-$ vs MCI$^+$",
    "HCneg_vs_ADpos": r"HC$^-$ vs AD$^+$",
    "MCIpos_vs_ADpos": r"MCI$^+$ vs AD$^+$",
}

# Existing result folders in your repo
TURBU_RESULTS_ROOT = "final_3d_gs_classification_turbu_sch1000"
HOPF_RESULTS_ROOT = "final_3d_gs_classification_turbu_hopf_sch1000_constrained"

# Existing feature-set folder names in your repo
TURBU_FEATURE_SET = "all_features_turbu_combat"
HOPF_FEATURE_SET = "all_features_turbu_hopf_combat"

CLASSIFIER = "LogReg"

# Plot appearance
COLOR_TURBU = "0.45"   # grey
COLOR_HOPF = "navy"
DIAG_COLOR = "0.7"
LW = 2.5
BAND_ALPHA = 0.18


def format_auc_label(roc_data):
    """Format the stored mean +/- SD AUC summary for a plot legend."""
    stored_auc = str(roc_data["auc"])
    if "+-" in stored_auc:
        mean_auc, sd_auc = (value.strip() for value in stored_auc.split("+-", 1))
        return rf"AUC = {mean_auc} $\pm$ {sd_auc}"
    return f"AUC = {stored_auc}"


def load_avg_roc_data(results_root: Path, scenario: str, feature_set: str):
    """
    Load precomputed average ROC data from a scenario/feature_set folder.
    """
    roc_path = (
        results_root
        / CLASSIFIER
        / scenario
        / feature_set
        / "avg_roc_data_test.joblib"
    )

    if not roc_path.exists():
        raise FileNotFoundError(f"Missing ROC data file:\n{roc_path}")

    return joblib.load(roc_path)


def style_axis(ax):
    """
    Apply clean journal-style axis formatting.
    """
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_linewidth(1)
    ax.spines["bottom"].set_linewidth(1)
    ax.tick_params(width=1)
    ax.grid(False)


def add_single_comparison(ax, roc_turbu, roc_hopf, title: str):
    """
    Overlay Turbulence and Turbulence+Perturbation ROC curves on one axis.
    """
    for roc_data, label, color in (
        (roc_turbu, "Turbulence", COLOR_TURBU),
        (roc_hopf, "Turbulence + Hopf", COLOR_HOPF),
    ):
        fpr = roc_data["fpr"]
        tpr = roc_data["tpr"]
        ax.plot(
            fpr,
            tpr,
            color=color,
            lw=LW,
            label=f"{label} ({format_auc_label(roc_data)})",
        )
        tpr_std = roc_data.get("std_tpr")
        if tpr_std is not None:
            ax.fill_between(
                fpr,
                (tpr - tpr_std).clip(min=0),
                (tpr + tpr_std).clip(max=1),
                color=color,
                alpha=BAND_ALPHA,
                linewidth=0,
            )

    # Draw the chance line only once after adding both model curves.
    ax.plot([0, 1], [0, 1], color=DIAG_COLOR, lw=1.5, linestyle="--")

    ax.set_title(title)
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.set_xlabel("False Positive Rate")
    ax.set_ylabel("True Positive Rate")
    ax.legend(loc="lower right", frameon=False)
    ax.set_box_aspect(1)
    style_axis(ax)


def save_single_panel_figure(path_repo: Path, output_dir: Path, scenario: str):
    """
    Save one standalone ROC comparison figure for a single contrast.
    """
    turbu_root = path_repo / "Results" / TURBU_RESULTS_ROOT
    hopf_root = path_repo / "Results" / HOPF_RESULTS_ROOT

    roc_turbu = load_avg_roc_data(turbu_root, scenario, TURBU_FEATURE_SET)
    roc_hopf = load_avg_roc_data(hopf_root, scenario, HOPF_FEATURE_SET)

    fig, ax = plt.subplots(figsize=(5.2, 5.2))
    add_single_comparison(ax, roc_turbu, roc_hopf, SCENARIO_TITLES[scenario])
    fig.tight_layout()

    output_stem = output_dir / f"roc_compare_{scenario.lower()}"
    fig.savefig(output_stem.with_suffix(".pdf"), dpi=600, bbox_inches="tight")
    fig.savefig(output_stem.with_suffix(".svg"), dpi=600, bbox_inches="tight")
    fig.savefig(output_stem.with_suffix(".png"), dpi=300, bbox_inches="tight")
    plt.close(fig)


def save_four_panel_figure(path_repo: Path, output_dir: Path):
    """
    Save a 2x2 panel ROC comparison figure covering all contrasts.
    """
    turbu_root = path_repo / "Results" / TURBU_RESULTS_ROOT
    hopf_root = path_repo / "Results" / HOPF_RESULTS_ROOT

    fig, axes = plt.subplots(2, 2, figsize=(10.5, 10.0))
    axes = axes.ravel()

    for ax, scenario in zip(axes, SCENARIOS):
        roc_turbu = load_avg_roc_data(turbu_root, scenario, TURBU_FEATURE_SET)
        roc_hopf = load_avg_roc_data(hopf_root, scenario, HOPF_FEATURE_SET)
        add_single_comparison(ax, roc_turbu, roc_hopf, SCENARIO_TITLES[scenario])

    fig.tight_layout()

    output_stem = output_dir / "roc_compare_all_contrasts"
    fig.savefig(output_stem.with_suffix(".pdf"), dpi=600, bbox_inches="tight")
    fig.savefig(output_stem.with_suffix(".svg"), dpi=600, bbox_inches="tight")
    fig.savefig(output_stem.with_suffix(".png"), dpi=300, bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    output_dir = PATH_REPO / "Results" / "final_3d_gs_classification_compare_roc"
    output_dir.mkdir(parents=True, exist_ok=True)

    # Save one figure per contrast
    for scenario in SCENARIOS:
        save_single_panel_figure(PATH_REPO, output_dir, scenario)

    # Save one combined 2x2 figure
    save_four_panel_figure(PATH_REPO, output_dir)

    print(f"ROC comparison plots saved to:\n{output_dir}")
