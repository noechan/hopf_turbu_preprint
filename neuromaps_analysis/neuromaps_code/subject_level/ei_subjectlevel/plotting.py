"""Publication figure for participant-level E:I--turbulence coupling."""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


# Match the sequential ``Blues`` palette used for the turbulence boxplots in
# manuscript Figure 2.  The light-to-dark ordering encodes progression along
# the amyloid-clinical continuum without assigning an unrelated hue to each
# diagnostic group.
STAGE_COLORS = ["#D0E1F2", "#94C4DF", "#4A98C9", "#1764AB"]
SHORT_LABELS = ["HC-", "HC+", "MCI+", "AD+"]
BOX_EDGE_COLOR = "#6B6B6B"


def draw_stage_boxplot(
    axis: plt.Axes,
    values: list[np.ndarray],
    positions: np.ndarray | None = None,
) -> dict[str, list]:
    """Draw a Figure-2-style boxplot ordered from HC- through AD+."""

    if len(values) != len(STAGE_COLORS):
        raise ValueError(f"Expected four staging groups; received {len(values)}.")
    if positions is None:
        positions = np.arange(1, len(values) + 1)
    box = axis.boxplot(
        values,
        positions=positions,
        widths=0.55,
        showfliers=True,
        patch_artist=True,
        boxprops={"edgecolor": BOX_EDGE_COLOR, "linewidth": 0.8},
        medianprops={"color": BOX_EDGE_COLOR, "linewidth": 1.0},
        whiskerprops={"color": BOX_EDGE_COLOR, "linewidth": 0.8},
        capprops={"color": BOX_EDGE_COLOR, "linewidth": 0.8},
        flierprops={
            "marker": "o",
            "markersize": 3.0,
            "markerfacecolor": "white",
            "markeredgecolor": BOX_EDGE_COLOR,
            "markeredgewidth": 0.6,
            "alpha": 1.0,
        },
    )
    for patch, color in zip(box["boxes"], STAGE_COLORS, strict=True):
        patch.set_facecolor(color)
        patch.set_alpha(1.0)
    return box


def plot_participant_results(
    participant_table: pd.DataFrame,
    group_order: list[str],
    display_names: dict[str, str],
    primary_trend: dict[str, object],
    sensitivity_trend: dict[str, object],
    pdf_path: Path,
    png_path: Path,
    seed: int,
) -> None:
    """Plot raw Pearson coupling and spin-normalized coupling by stage."""

    # Kept in the signature because the configuration supplies publication names.
    del display_names
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
    # Retain the argument for reproducibility compatibility with earlier runs.
    del seed
    figure, axes = plt.subplots(1, 2, figsize=(10.2, 4.7))
    outcomes = [
        ("pearson_r", "Participant E:I-turbulence correlation", primary_trend),
        ("spin_z", "Spatially normalized coupling (z)", sensitivity_trend),
    ]
    for panel, (axis, (column, ylabel, trend)) in enumerate(zip(axes, outcomes)):
        values = [
            participant_table.loc[
                participant_table["Group"] == group, column
            ].to_numpy()
            for group in group_order
        ]
        draw_stage_boxplot(axis, values)
        axis.set_xticks(np.arange(1, 5), SHORT_LABELS)
        axis.set_xlabel("Amyloid-staging group")
        axis.set_ylabel(ylabel)
        outcome_name = "Fisher-z" if column == "pearson_r" else "Spin-z"
        axis.set_title(
            f"{outcome_name} ordered trend: beta={float(trend['coefficient']):.3f}, "
            f"p$_{{perm}}$={float(trend['p_permutation_two_sided']):.3f}"
        )
        axis.text(
            -0.13,
            1.06,
            chr(ord("A") + panel),
            transform=axis.transAxes,
            fontsize=13,
            fontweight="bold",
            va="top",
        )
        axis.spines[["top", "right"]].set_visible(False)
    figure.suptitle(
        "Participant-level coupling between cortical E:I and turbulence at lambda=0.01",
        fontsize=13,
        fontweight="bold",
    )
    figure.tight_layout(rect=(0.01, 0.01, 0.99, 0.94))
    pdf_path.parent.mkdir(parents=True, exist_ok=True)
    png_path.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(pdf_path, bbox_inches="tight")
    figure.savefig(png_path, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(figure)
