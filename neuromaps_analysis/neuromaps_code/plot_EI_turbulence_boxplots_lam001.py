#!/usr/bin/env python3
"""Plot spatial-null r distributions for HC- and HC- baseline contrasts."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import tempfile


ROOT = Path(__file__).resolve().parents[1]
os.environ.setdefault(
    "MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "neuromaps_matplotlib")
)

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
import numpy as np
import pandas as pd


PDF_NAME = "Figure_EI_turbulence_boxplots_lam001.pdf"
PNG_NAME = "Figure_EI_turbulence_boxplots_lam001.png"
METADATA_NAME = "Figure_EI_turbulence_boxplots_lam001_metadata.json"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--pdf-dir", type=Path, default=ROOT / "output" / "pdf"
    )
    parser.add_argument(
        "--figure-dir", type=Path, default=ROOT / "output" / "figures"
    )
    return parser.parse_args()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> None:
    args = parse_args()
    hc_results = ROOT / "results" / "ei_hc_turbulence_lam001"
    contrast_results = ROOT / "results" / "ei_hc_baseline_contrasts_lam001"
    hc_summary_path = hc_results / "EI_HC_ABetaNeg_turbulence_lam001_summary.csv"
    contrast_summary_path = (
        contrast_results / "EI_HC_ABetaNeg_baseline_contrasts_lam001_summary.csv"
    )
    for path in (hc_summary_path, contrast_summary_path):
        if not path.is_file():
            raise FileNotFoundError(f"Missing analysis summary: {path}")

    hc_summary = pd.read_csv(hc_summary_path)
    contrast_summary = pd.read_csv(contrast_summary_path)
    if len(hc_summary) != 1 or len(contrast_summary) != 3:
        raise ValueError("Expected one HC-only row and three HC-contrast rows.")

    categories = [
        {
            "id": "HC_ABetaNeg",
            "label": "HC- only",
            "definition": "HC amyloid-negative group mean",
            "r": float(hc_summary.iloc[0]["empirical_r"]),
            "p": float(hc_summary.iloc[0]["p_spin"]),
            "n_perm": int(hc_summary.iloc[0]["n_perm"]),
            "null_path": hc_results / hc_summary.iloc[0]["null_distribution_file"],
        }
    ]
    labels = {
        "HCneg_minus_HCpos": "HC- minus HC+",
        "HCneg_minus_MCIpos": "HC- minus MCI+",
        "HCneg_minus_ADpos": "HC- minus AD+",
    }
    for _, row in contrast_summary.iterrows():
        contrast = row["contrast"]
        categories.append(
            {
                "id": contrast,
                "label": labels[contrast],
                "definition": row["operation"],
                "r": float(row["empirical_r"]),
                "p": float(row["p_spin"]),
                "n_perm": int(row["n_perm"]),
                "null_path": contrast_results / row["null_distribution_file"],
            }
        )

    null_distributions: list[np.ndarray] = []
    for category in categories:
        null_path = category["null_path"]
        if not null_path.is_file():
            raise FileNotFoundError(f"Missing null distribution: {null_path}")
        null_values = np.load(null_path)
        if (
            null_values.shape != (category["n_perm"],)
            or not np.isfinite(null_values).all()
        ):
            raise ValueError(f"Invalid null distribution: {null_path}")
        recomputed_p = (
            np.sum(np.abs(null_values) >= abs(category["r"])) + 1
        ) / (category["n_perm"] + 1)
        if not np.isclose(recomputed_p, category["p"]):
            raise ValueError(
                f"Stored p for {category['id']} does not match its null distribution."
            )
        null_distributions.append(null_values)

    plt.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
            "font.size": 10,
            "axes.labelsize": 11,
            "axes.titlesize": 13,
            "xtick.labelsize": 9,
            "ytick.labelsize": 9,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "axes.grid": False,
        }
    )

    colors = ["#4D4D4D", "#0072B2", "#E69F00", "#D55E00"]
    positions = np.arange(1, 5)
    figure, axis = plt.subplots(figsize=(8.5, 5.6))
    boxplot = axis.boxplot(
        null_distributions,
        positions=positions,
        widths=0.58,
        whis=(2.5, 97.5),
        showfliers=False,
        patch_artist=True,
        medianprops={"color": "#333333", "linewidth": 1.2},
        whiskerprops={"color": "#666666", "linewidth": 1.0},
        capprops={"color": "#666666", "linewidth": 1.0},
    )
    for box, color in zip(boxplot["boxes"], colors):
        box.set_facecolor(color)
        box.set_edgecolor(color)
        box.set_alpha(0.20)

    empirical = np.asarray([category["r"] for category in categories])
    axis.scatter(
        positions,
        empirical,
        c=colors,
        s=78,
        edgecolor="white",
        linewidth=0.9,
        zorder=4,
    )
    axis.plot(
        positions[1:],
        empirical[1:],
        color="#555555",
        linewidth=1.4,
        zorder=3,
    )
    axis.axvline(1.5, color="#AAAAAA", linestyle=":", linewidth=1.0)
    axis.axhline(0, color="#555555", linewidth=0.8)

    for x, category in zip(positions, categories):
        axis.annotate(
            f"r={category['r']:.3f}\np$_{{spin}}$={category['p']:.3f}",
            (x, category["r"]),
            xytext=(0, 11),
            textcoords="offset points",
            ha="center",
            va="bottom",
            fontsize=9,
        )

    axis.set_xticks(
        positions,
        ["HC-\nonly", "HC- minus\nHC+", "HC- minus\nMCI+", "HC- minus\nAD+"],
    )
    axis.set_ylabel("Pearson's r")
    axis.set_title("E:I association with turbulence at lambda=0.01", pad=13)
    lower = min(float(values.min()) for values in null_distributions) - 0.05
    upper = max(float(empirical.max()), max(float(v.max()) for v in null_distributions)) + 0.14
    axis.set_ylim(lower, upper)
    axis.legend(
        handles=[
            Patch(
                facecolor="#8FAABD",
                edgecolor="#667788",
                alpha=0.3,
                label="Spatial-null r distribution",
            ),
            Line2D(
                [0],
                [0],
                marker="o",
                color="none",
                markerfacecolor="#D55E00",
                markeredgecolor="white",
                markersize=8,
                label="Observed r",
            ),
        ],
        loc="lower left",
        frameon=False,
        fontsize=9,
    )
    axis.spines[["top", "right"]].set_visible(False)
    figure.text(
        0.5,
        0.01,
        "Boxes show the interquartile range; whiskers show the 2.5th-97.5th percentiles of 1,000 spatial rotations.",
        ha="center",
        va="bottom",
        fontsize=8,
        color="#444444",
    )
    figure.tight_layout(rect=(0.02, 0.05, 0.98, 0.98))

    pdf_dir = args.pdf_dir.expanduser().resolve()
    figure_dir = args.figure_dir.expanduser().resolve()
    pdf_dir.mkdir(parents=True, exist_ok=True)
    figure_dir.mkdir(parents=True, exist_ok=True)
    pdf_path = pdf_dir / PDF_NAME
    png_path = figure_dir / PNG_NAME
    metadata_path = figure_dir / METADATA_NAME
    figure.savefig(pdf_path, bbox_inches="tight")
    figure.savefig(png_path, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(figure)

    metadata = {
        "figure_definition": (
            "Boxplots of 1,000 spatial-null Pearson correlations with observed "
            "correlations overlaid for the HC- mean and three signed HC- contrasts."
        ),
        "box_definition": "IQR with whiskers at the 2.5th and 97.5th percentiles",
        "observed_r_is_a_distribution": False,
        "inputs": {
            "hc_summary": {
                "path": str(hc_summary_path),
                "sha256": sha256(hc_summary_path),
            },
            "contrast_summary": {
                "path": str(contrast_summary_path),
                "sha256": sha256(contrast_summary_path),
            },
            "null_distributions": [
                {"path": str(category["null_path"]), "sha256": sha256(category["null_path"])}
                for category in categories
            ],
        },
        "categories": [
            {key: value for key, value in category.items() if key != "null_path"}
            for category in categories
        ],
        "outputs": {"pdf": str(pdf_path), "png": str(png_path)},
    }
    with metadata_path.open("w", encoding="utf-8") as stream:
        json.dump(metadata, stream, indent=2)
        stream.write("\n")

    print(f"Saved boxplot PDF: {pdf_path}")
    print(f"Saved boxplot PNG: {png_path}")
    print(f"Saved figure metadata: {metadata_path}")


if __name__ == "__main__":
    main()
