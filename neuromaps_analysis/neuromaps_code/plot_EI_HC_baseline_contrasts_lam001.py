#!/usr/bin/env python3
"""Create the manuscript figure for E:I and HC- baseline turbulence contrasts."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import tempfile


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_CONFIG = ROOT / "config" / "ei_hc_baseline_contrasts_lam001.json"
DEFAULT_CACHE = Path(os.environ.get("XDG_CACHE_HOME", Path.home() / ".cache")) / "neuromaps"
os.environ.setdefault("NEUROMAPS_DATA", str(DEFAULT_CACHE))
os.environ.setdefault(
    "MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "neuromaps_matplotlib")
)

import matplotlib

matplotlib.use("Agg")
from matplotlib.colors import LinearSegmentedColormap
import matplotlib.pyplot as plt
from neuromaps import transforms
import numpy as np
import pandas as pd


PDF_NAME = "Figure_EI_HC_baseline_turbulence_contrasts_lam001.pdf"
PNG_NAME = "Figure_EI_HC_baseline_turbulence_contrasts_lam001.png"
METADATA_NAME = "Figure_EI_HC_baseline_turbulence_contrasts_lam001_metadata.json"
COLORS = ["#0072B2", "#E69F00", "#D55E00"]
SHORT_LABELS = ["HC- minus HC+", "HC- minus MCI+", "HC- minus AD+"]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    parser.add_argument(
        "--pdf-dir", type=Path, default=ROOT / "output" / "pdf"
    )
    parser.add_argument(
        "--figure-dir", type=Path, default=ROOT / "output" / "figures"
    )
    return parser.parse_args()


def resolve_from_config(config_path: Path, value: str) -> Path:
    path = Path(value).expanduser()
    if not path.is_absolute():
        path = config_path.parent / path
    return path.resolve()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def surface_to_array(surface: tuple) -> np.ndarray:
    data = np.hstack(
        [
            np.asarray(hemisphere.agg_data(), dtype=float).squeeze()
            for hemisphere in surface
        ]
    )
    if data.ndim != 1 or not np.isfinite(data).all():
        raise ValueError("Surface data must be a finite one-dimensional array.")
    return data


def standardize(values: np.ndarray) -> np.ndarray:
    scale = values.std(ddof=0)
    if not np.isfinite(scale) or scale == 0:
        raise ValueError("Cannot standardize a constant or non-finite map.")
    return (values - values.mean()) / scale


def density_cmap(color: str) -> LinearSegmentedColormap:
    return LinearSegmentedColormap.from_list("density", ["#F7F7F7", color])


def main() -> None:
    args = parse_args()
    config_path = args.config.expanduser().resolve()
    with config_path.open(encoding="utf-8") as stream:
        config = json.load(stream)

    results_dir = resolve_from_config(config_path, config["output_dir"])
    contrast_dir = resolve_from_config(config_path, config["contrast_output_dir"])
    summary_path = results_dir / "EI_HC_ABetaNeg_baseline_contrasts_lam001_summary.csv"
    if not summary_path.is_file():
        raise FileNotFoundError(f"Missing correlation summary: {summary_path}")
    summary = pd.read_csv(summary_path)
    configured_contrasts = list(config["contrasts"])
    if summary["contrast"].tolist() != configured_contrasts:
        raise ValueError("Summary contrast order does not match the configuration.")

    ei_path = resolve_from_config(config_path, config["ei_map"])
    atlas = config["surface"]["atlas"]
    density = config["surface"]["density"]
    print(f"Transforming E:I map to {atlas} {density}")
    ei = surface_to_array(transforms.mni152_to_fsaverage(str(ei_path), density))

    ei_z_by_contrast: list[np.ndarray] = []
    contrast_z: list[np.ndarray] = []
    nulls: list[np.ndarray] = []
    contrast_paths: list[Path] = []
    for _, row in summary.iterrows():
        contrast = row["contrast"]
        contrast_path = contrast_dir / (
            f"turbu_{contrast}_lam001_sch1000_raw_2mm.nii.gz"
        )
        if not contrast_path.is_file():
            raise FileNotFoundError(f"Missing contrast map: {contrast_path}")
        print(f"Transforming {contrast}")
        values = surface_to_array(
            transforms.mni152_to_fsaverage(str(contrast_path), density)
        )
        if values.shape != ei.shape:
            raise ValueError(f"Surface size mismatch for {contrast}.")
        valid = np.logical_and(~np.isclose(values, 0), ~np.isclose(ei, 0))
        if valid.sum() < 3:
            raise ValueError(f"Too few jointly nonzero vertices for {contrast}.")
        observed_r = float(np.corrcoef(values[valid], ei[valid])[0, 1])
        if not np.isclose(observed_r, float(row["empirical_r"]), atol=1e-7):
            raise ValueError(
                f"Recomputed r for {contrast} ({observed_r}) differs from summary."
            )
        ei_z_by_contrast.append(standardize(ei[valid]))
        contrast_z.append(standardize(values[valid]))
        contrast_paths.append(contrast_path)

        null_path = results_dir / row["null_distribution_file"]
        null = np.load(null_path)
        if null.shape != (int(row["n_perm"]),) or not np.isfinite(null).all():
            raise ValueError(f"Invalid null distribution for {contrast}.")
        nulls.append(null)

    plt.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
            "font.size": 9,
            "axes.labelsize": 9,
            "axes.titlesize": 10,
            "xtick.labelsize": 8,
            "ytick.labelsize": 8,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "axes.grid": False,
        }
    )

    figure = plt.figure(figsize=(11.4, 7.2), constrained_layout=True)
    grid = figure.add_gridspec(2, 3, height_ratios=[1.05, 0.90])
    top_axes = [figure.add_subplot(grid[0, index]) for index in range(3)]
    summary_axis = figure.add_subplot(grid[1, :])

    for index, (axis, x, y, color, label) in enumerate(
        zip(top_axes, ei_z_by_contrast, contrast_z, COLORS, SHORT_LABELS)
    ):
        row = summary.iloc[index]
        axis.hexbin(
            x,
            y,
            gridsize=44,
            bins="log",
            mincnt=1,
            cmap=density_cmap(color),
            linewidths=0,
            rasterized=True,
        )
        x_line = np.linspace(np.quantile(x, 0.005), np.quantile(x, 0.995), 100)
        slope, intercept = np.polyfit(x, y, 1)
        axis.plot(x_line, intercept + slope * x_line, color=color, linewidth=2)
        axis.axhline(0, color="#777777", linewidth=0.7, linestyle="--")
        axis.axvline(0, color="#777777", linewidth=0.7, linestyle="--")
        axis.set_xlabel("E:I expression (z score)")
        if index == 0:
            axis.set_ylabel("Turbulence contrast (z score)")
        axis.set_title(label)
        axis.text(
            0.04,
            0.96,
            f"r = {float(row['empirical_r']):.3f}\n"
            f"p$_{{spin}}$ = {float(row['p_spin']):.3f}\n"
            f"q$_{{FDR}}$ = {float(row['p_fdr_bh']):.3f}",
            transform=axis.transAxes,
            ha="left",
            va="top",
            fontsize=8.5,
            bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.82, "pad": 2},
        )
        axis.text(
            -0.15,
            1.08,
            chr(ord("A") + index),
            transform=axis.transAxes,
            fontsize=13,
            fontweight="bold",
            va="top",
        )
        axis.spines[["top", "right"]].set_visible(False)

    positions = np.arange(1, 4)
    violin = summary_axis.violinplot(
        nulls,
        positions=positions,
        widths=0.58,
        showmeans=False,
        showmedians=True,
        showextrema=False,
    )
    for body, color in zip(violin["bodies"], COLORS):
        body.set_facecolor(color)
        body.set_edgecolor("none")
        body.set_alpha(0.20)
    violin["cmedians"].set_color("#555555")
    violin["cmedians"].set_linewidth(1)

    empirical = summary["empirical_r"].to_numpy(dtype=float)
    summary_axis.plot(positions, empirical, color="#333333", linewidth=1.5, zorder=3)
    for x, value, color in zip(positions, empirical, COLORS):
        summary_axis.scatter(
            x,
            value,
            s=65,
            color=color,
            edgecolor="white",
            linewidth=0.8,
            zorder=4,
        )
    summary_axis.axhline(0, color="#777777", linewidth=0.8)
    summary_axis.set_xticks(positions, SHORT_LABELS)
    summary_axis.set_ylabel("Pearson's r")
    summary_axis.set_title(
        "Observed spatial correlations and Alexander-Bloch null distributions"
    )
    summary_axis.set_ylim(
        min(float(np.min(null)) for null in nulls) - 0.05,
        max(float(empirical.max()), max(float(np.max(null)) for null in nulls)) + 0.12,
    )
    for x, row in zip(positions, summary.itertuples(index=False)):
        summary_axis.annotate(
            f"r={row.empirical_r:.3f}",
            (x, row.empirical_r),
            xytext=(0, 9),
            textcoords="offset points",
            ha="center",
            fontsize=8.5,
        )
    summary_axis.text(
        -0.05,
        1.08,
        "D",
        transform=summary_axis.transAxes,
        fontsize=13,
        fontweight="bold",
        va="top",
    )
    summary_axis.spines[["top", "right"]].set_visible(False)

    figure.suptitle(
        "Cortical E:I organization tracks turbulence changes from the HC- baseline",
        fontsize=14,
        fontweight="bold",
    )
    figure.text(
        0.5,
        -0.012,
        "Physical scale: lambda=0.01. Hexagons show log vertex density; p values use 1,000 spatial rotations.",
        ha="center",
        va="bottom",
        fontsize=8,
        color="#444444",
    )

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
            "Manuscript visualization of spatial E:I associations with three signed "
            "HC amyloid-negative baseline turbulence contrasts at lambda=0.01."
        ),
        "statistics_recomputed": False,
        "surface_values_recomputed_for_display": True,
        "configuration": {"path": str(config_path), "sha256": sha256(config_path)},
        "summary": {"path": str(summary_path), "sha256": sha256(summary_path)},
        "ei_map": {"path": str(ei_path), "sha256": sha256(ei_path)},
        "contrast_maps": [
            {"path": str(path), "sha256": sha256(path)} for path in contrast_paths
        ],
        "outputs": {"pdf": str(pdf_path), "png": str(png_path)},
    }
    with metadata_path.open("w", encoding="utf-8") as stream:
        json.dump(metadata, stream, indent=2)
        stream.write("\n")

    print(f"Saved manuscript PDF: {pdf_path}")
    print(f"Saved manuscript PNG: {png_path}")
    print(f"Saved figure metadata: {metadata_path}")


if __name__ == "__main__":
    main()
