#!/usr/bin/env python3
"""Create the manuscript E:I surface and participant spin-coupling figure."""

from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path
import tempfile


NEUROMAPS_ROOT = Path(__file__).resolve().parents[2]
CONFIG_PATH = NEUROMAPS_ROOT / "config" / "ei_subjectlevel_coupling_lam001.json"
RESULTS_DIR = NEUROMAPS_ROOT / "results" / "ei_subjectlevel_coupling_lam001"
PARTICIPANT_PATH = (
    RESULTS_DIR / "participant_EI_turbulence_coupling_lam001_N145.csv"
)
TREND_PATH = RESULTS_DIR / "ordered_trend_EI_turbulence_lam001_N145.csv"
COMPARISON_PATH = RESULTS_DIR / "confirmatory_tests_EI_turbulence_lam001_N145.csv"
PDF_PATH = (
    NEUROMAPS_ROOT
    / "output"
    / "pdf"
    / "Figure_EI_surface_spin_coupling_lam001.pdf"
)
PNG_PATH = (
    NEUROMAPS_ROOT
    / "output"
    / "figures"
    / "Figure_EI_surface_spin_coupling_lam001.png"
)
METADATA_PATH = (
    NEUROMAPS_ROOT
    / "output"
    / "figures"
    / "Figure_EI_surface_spin_coupling_lam001_metadata.json"
)

DEFAULT_CACHE = NEUROMAPS_ROOT / "neuromaps-data" / "cache"
os.environ.setdefault("NEUROMAPS_DATA", str(DEFAULT_CACHE))
os.environ.setdefault(
    "MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "neuromaps_matplotlib")
)

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.cm import ScalarMappable
from matplotlib.colors import Normalize
from neuromaps import transforms
from neuromaps.datasets import fetch_atlas
from nilearn import plotting
from nilearn import surface as nilearn_surface
import numpy as np
import pandas as pd
from scipy.spatial import cKDTree

from ei_subjectlevel.plotting import SHORT_LABELS, draw_stage_boxplot


GROUP_ORDER = ["HC_ABneg", "HC_ABpos", "MCI_ABpos", "AD_ABpos"]
GROUP_COUNTS = {"HC_ABneg": 51, "HC_ABpos": 37, "MCI_ABpos": 31, "AD_ABpos": 26}


def resolve_from_config(value: str) -> Path:
    path = Path(value).expanduser()
    if not path.is_absolute():
        path = CONFIG_PATH.parent / path
    return path.resolve()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def gifti_data(image) -> np.ndarray:
    values = np.asarray(image.agg_data(), dtype=float).squeeze()
    if values.ndim != 1 or not np.isfinite(values).all():
        raise ValueError("E:I surface data must be finite and one-dimensional.")
    return values


def load_inputs() -> tuple[dict, pd.DataFrame, dict[str, float], dict[str, float]]:
    with CONFIG_PATH.open(encoding="utf-8") as stream:
        config = json.load(stream)
    participants = pd.read_csv(PARTICIPANT_PATH)
    trends = pd.read_csv(TREND_PATH)
    comparisons = pd.read_csv(COMPARISON_PATH)

    counts = participants["Group"].value_counts().to_dict()
    if counts != GROUP_COUNTS or len(participants) != 145:
        raise ValueError(f"Unexpected participant cohort: N={len(participants)}, counts={counts}")
    if not np.isfinite(participants["spin_z"].to_numpy(dtype=float)).all():
        raise ValueError("Participant spin-z values contain non-finite entries.")

    spin_trend = trends.loc[trends["outcome"] == "spin_z"]
    if len(spin_trend) != 1:
        raise ValueError("Expected exactly one spin-z ordered trend result.")
    trend = spin_trend.iloc[0].to_dict()
    ad_comparison = comparisons.loc[
        (comparisons["outcome"] == "spin_z")
        & (comparisons["comparison_group"] == "AD_ABpos")
    ]
    if len(ad_comparison) != 1:
        raise ValueError("Expected exactly one HC- versus AD+ spin-z comparison.")
    return config, participants, trend, ad_comparison.iloc[0].to_dict()


def render_ei_surfaces(
    figure: plt.Figure,
    grid,
    ei_hemispheres: tuple[np.ndarray, np.ndarray],
    atlas,
    vmin: float,
    vmax: float,
) -> None:
    specifications = [
        (0, 0, "left", "lateral", "Left lateral"),
        (0, 1, "right", "lateral", "Right lateral"),
        (1, 0, "left", "medial", "Left medial"),
        (1, 1, "right", "medial", "Right medial"),
    ]
    for row, column, hemisphere, view, title in specifications:
        index = 0 if hemisphere == "left" else 1
        axis = figure.add_subplot(grid[row, column], projection="3d")
        values = prepare_surface_values(
            ei_hemispheres[index],
            atlas.medial[index],
            atlas.sphere[index],
            vmin,
            vmax,
        )
        plotting.plot_surf_stat_map(
            surf_mesh=atlas.inflated[index],
            stat_map=values,
            bg_map=atlas.sulc[index],
            hemi=hemisphere,
            view=view,
            axes=axis,
            figure=figure,
            cmap="viridis",
            colorbar=False,
            symmetric_cbar=False,
            vmin=vmin,
            vmax=vmax,
            bg_on_data=True,
            darkness=0.55,
            alpha=1.0,
        )
        axis.set_title(title, fontsize=8.5, pad=-2)
        axis.set_axis_off()
        for collection in axis.collections:
            collection.set_rasterized(True)

    color_axis = figure.add_subplot(grid[2, :])
    colorbar = figure.colorbar(
        ScalarMappable(norm=Normalize(vmin=vmin, vmax=vmax), cmap="viridis"),
        cax=color_axis,
        orientation="horizontal",
    )
    colorbar.set_label("E:I expression ratio (raw values)", labelpad=3)
    colorbar.ax.tick_params(labelsize=8, length=2)
    colorbar.outline.set_linewidth(0.6)


def prepare_surface_values(
    values: np.ndarray,
    medial_path: Path,
    sphere_path: Path,
    vmin: float,
    vmax: float,
) -> np.ndarray:
    """Mask the medial wall and fill volume-projection edge gaps for display."""
    cortex = np.asarray(nilearn_surface.load_surf_data(medial_path)).squeeze() > 0
    sphere_coordinates = nilearn_surface.load_surf_mesh(sphere_path).coordinates
    valid = cortex & np.isfinite(values) & (values >= vmin - 1e-6)
    missing = cortex & ~valid
    if valid.sum() < 1000:
        raise ValueError("Too few valid E:I vertices for surface rendering.")
    prepared = values.copy()
    if missing.any():
        nearest = cKDTree(sphere_coordinates[valid]).query(
            sphere_coordinates[missing], k=1
        )[1]
        prepared[missing] = prepared[valid][nearest]
    prepared = np.clip(prepared, vmin, vmax)
    prepared[~cortex] = np.nan
    return prepared


def render_spin_panel(
    axis: plt.Axes,
    participants: pd.DataFrame,
    trend: dict[str, float],
    seed: int,
) -> None:
    # Retain the argument for reproducibility compatibility with earlier runs.
    del seed
    values = [
        participants.loc[participants["Group"] == group, "spin_z"].to_numpy(dtype=float)
        for group in GROUP_ORDER
    ]
    draw_stage_boxplot(axis, values)

    labels = [f"{label}\n(n={GROUP_COUNTS[group]})" for label, group in zip(SHORT_LABELS, GROUP_ORDER)]
    axis.set_xticks(np.arange(1, 5), labels)
    axis.set_xlim(0.5, 4.5)
    axis.set_ylim(bottom=0)
    axis.set_ylabel("Spatially normalized E:I-turbulence coupling (z)")
    axis.set_xlabel("Amyloid-clinical stage")
    axis.set_title(
        "Progressive reduction in correspondence\n"
        "with the normative E:I architecture",
        fontsize=11,
        fontweight="bold",
        pad=30,
    )
    axis.text(
        0.5,
        1.035,
        rf"Ordered trend: $\beta_{{stage}}$ = {trend['coefficient']:.3f}, "
        rf"95% CI [{trend['ci95_low_hc3']:.3f}, {trend['ci95_high_hc3']:.3f}], "
        rf"$p_{{perm}}$={trend['p_permutation_two_sided']:.4f}",
        transform=axis.transAxes,
        ha="center",
        va="bottom",
        fontsize=8.5,
    )
    axis.spines[["top", "right"]].set_visible(False)
    axis.tick_params(axis="both", labelsize=8.5)


def main() -> None:
    config, participants, trend, ad_comparison = load_inputs()
    ei_path = resolve_from_config(config["ei_nifti"])
    density = config["surface"]["density"]
    ei_surface = transforms.mni152_to_fsaverage(str(ei_path), density)
    ei_hemispheres = tuple(gifti_data(image) for image in ei_surface)
    ei_values_path = resolve_from_config(config["ei_values"])
    ei_table = pd.read_csv(ei_values_path)
    ei_values = ei_table[config["ei_column"]].to_numpy(dtype=float)
    if (
        ei_values.size != 500
        or not np.isfinite(ei_values).all()
        or "hemisphere" not in ei_table
        or not ei_table["hemisphere"].eq("L").all()
    ):
        raise ValueError("Expected 500 finite left-hemisphere raw E:I parcel values.")
    vmin, vmax = float(ei_values.min()), float(ei_values.max())
    atlas = fetch_atlas(config["surface"]["atlas"], density)

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
    figure = plt.figure(figsize=(11.7, 5.5), facecolor="white")
    outer = figure.add_gridspec(
        1,
        2,
        width_ratios=[1.42, 1.0],
        left=0.025,
        right=0.985,
        bottom=0.14,
        top=0.86,
        wspace=0.13,
    )
    surface_grid = outer[0].subgridspec(
        3,
        2,
        height_ratios=[1.0, 1.0, 0.09],
        hspace=-0.28,
        wspace=-0.18,
    )
    render_ei_surfaces(figure, surface_grid, ei_hemispheres, atlas, vmin, vmax)
    spin_axis = figure.add_subplot(outer[1])
    render_spin_panel(spin_axis, participants, trend, seed=1234)

    figure.text(0.016, 0.955, "A", fontsize=14, fontweight="bold", va="top")
    figure.text(
        0.285,
        0.94,
        "Normative cortical E:I architecture",
        fontsize=11,
        fontweight="bold",
        ha="center",
        va="top",
    )
    figure.text(0.602, 0.955, "B", fontsize=14, fontweight="bold", va="top")

    PDF_PATH.parent.mkdir(parents=True, exist_ok=True)
    PNG_PATH.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(PDF_PATH, dpi=600, bbox_inches="tight", facecolor="white")
    figure.savefig(PNG_PATH, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(figure)

    metadata = {
        "figure_definition": (
            "Normative cortical E:I surface map and participant-level spatially "
            "normalized E:I-turbulence coupling at physical lambda=0.01"
        ),
        "inputs": {
            "configuration": {"path": str(CONFIG_PATH), "sha256": sha256(CONFIG_PATH)},
            "ei_values": {
                "path": str(ei_values_path),
                "sha256": sha256(ei_values_path),
            },
            "ei_nifti": {"path": str(ei_path), "sha256": sha256(ei_path)},
            "participant_scores": {
                "path": str(PARTICIPANT_PATH),
                "sha256": sha256(PARTICIPANT_PATH),
            },
            "ordered_trends": {"path": str(TREND_PATH), "sha256": sha256(TREND_PATH)},
            "planned_comparisons": {
                "path": str(COMPARISON_PATH),
                "sha256": sha256(COMPARISON_PATH),
            },
        },
        "surface": {
            "atlas": config["surface"]["atlas"],
            "density": density,
            "views": ["left_lateral", "right_lateral", "left_medial", "right_medial"],
            "color_limits": {
                "policy": "full_raw_surface_range",
                "vmin": float(vmin),
                "vmax": float(vmax),
            },
            "display_edge_handling": (
                "Medial-wall vertices are masked; cortical vertices below the raw "
                "parcel minimum after volume-to-surface interpolation are filled "
                "from the nearest valid vertex on the spherical mesh."
            ),
        },
        "panel_b": {
            "outcome": "spin_z",
            "group_counts": GROUP_COUNTS,
            "ordered_trend": trend,
            "hc_negative_vs_ad_positive": ad_comparison,
        },
        "outputs": {"pdf": str(PDF_PATH), "png": str(PNG_PATH)},
    }
    with METADATA_PATH.open("w", encoding="utf-8") as stream:
        json.dump(metadata, stream, indent=2)
        stream.write("\n")

    print(f"Saved PDF: {PDF_PATH}")
    print(f"Saved PNG: {PNG_PATH}")
    print(f"Saved metadata: {METADATA_PATH}")


if __name__ == "__main__":
    main()
