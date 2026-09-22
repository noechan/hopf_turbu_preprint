#!/usr/bin/env python3
"""Correlate primary E:I with the HC mean turbulence map at lambda=0.01."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import tempfile


NEUROMAPS_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_CONFIG = NEUROMAPS_ROOT / "config" / "ei_hc_turbulence_lam001.json"
DEFAULT_CACHE = NEUROMAPS_ROOT / "neuromaps-data" / "cache"
os.environ.setdefault("NEUROMAPS_DATA", str(DEFAULT_CACHE))
os.environ.setdefault(
    "MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "neuromaps_matplotlib")
)

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import neuromaps
from neuromaps import nulls, transforms
from neuromaps.stats import compare_images
import nibabel
import numpy as np
import pandas as pd


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    parser.add_argument("--n-perm", type=int, help="Override configured rotations.")
    parser.add_argument("--seed", type=int, help="Override configured random seed.")
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


def surface_to_array(surface: tuple, hemisphere: str) -> np.ndarray:
    hemispheres = [
        np.asarray(image.agg_data(), dtype=float).squeeze() for image in surface
    ]
    if hemisphere != "left":
        raise ValueError("This paper-style analysis requires surface.hemisphere='left'.")
    data = hemispheres[0]
    if data.ndim != 1 or not np.isfinite(data).all():
        raise ValueError("Transformed surface map must be a finite one-dimensional array.")
    return data


def left_spin_indices(
    atlas: str, density: str, n_vertices: int, n_perm: int, seed: int
) -> np.ndarray:
    full = np.asarray(
        nulls.alexander_bloch(
            None, atlas=atlas, density=density, n_perm=n_perm, seed=seed
        )
    )
    if full.shape[0] != 2 * n_vertices:
        raise ValueError(f"Unexpected bilateral spin-index shape: {full.shape}")
    left = full[:n_vertices]
    if left.min() < 0 or left.max() >= n_vertices:
        raise ValueError("Left-hemisphere spins reference right-hemisphere vertices.")
    return left


def main() -> None:
    args = parse_args()
    config_path = args.config.expanduser().resolve()
    with config_path.open(encoding="utf-8") as stream:
        config = json.load(stream)

    turbulence = config["turbulence"]
    ei_path = resolve_from_config(config_path, config["ei_map"])
    turbulence_path = resolve_from_config(config_path, turbulence["group_mean_map"])
    output_dir = resolve_from_config(config_path, config["output_dir"])
    figure_path = resolve_from_config(config_path, config["figure_file"])
    output_dir.mkdir(parents=True, exist_ok=True)
    figure_path.parent.mkdir(parents=True, exist_ok=True)
    for label, path in (("E:I", ei_path), ("HC turbulence", turbulence_path)):
        if not path.is_file():
            raise FileNotFoundError(f"Missing {label} map: {path}")

    atlas = config["surface"]["atlas"]
    density = config["surface"]["density"]
    hemisphere = config["surface"]["hemisphere"]
    source_parcels = int(config["surface"]["source_parcels"])
    if hemisphere != "left" or source_parcels != 500:
        raise ValueError("Statistical correlations must use the 500 left parcels.")
    metric = config["metric"]
    n_perm = args.n_perm if args.n_perm is not None else int(
        config["null_model"]["n_perm"]
    )
    seed = args.seed if args.seed is not None else int(config["null_model"]["seed"])
    if n_perm < 1:
        raise ValueError("n_perm must be a positive integer.")

    print(f"Transforming E:I map to {atlas} {density}: {ei_path}")
    ei_data = surface_to_array(
        transforms.mni152_to_fsaverage(str(ei_path), density), hemisphere
    )
    print(f"Transforming HC turbulence map to {atlas} {density}: {turbulence_path}")
    turbulence_data = surface_to_array(
        transforms.mni152_to_fsaverage(str(turbulence_path), density), hemisphere
    )
    if turbulence_data.shape != ei_data.shape:
        raise ValueError(
            f"Surface array mismatch: turbulence={turbulence_data.shape}, "
            f"E:I={ei_data.shape}"
        )

    print(f"Generating {n_perm} Alexander-Bloch rotations (seed={seed})")
    spin_indices = left_spin_indices(atlas, density, ei_data.size, n_perm, seed)
    rotated = turbulence_data[spin_indices]
    r_value, p_spin, null_distribution = compare_images(
        turbulence_data,
        ei_data,
        nulls=rotated,
        metric=metric,
        return_nulls=True,
    )
    null_distribution = np.asarray(null_distribution, dtype=float)
    if null_distribution.shape != (n_perm,) or not np.isfinite(null_distribution).all():
        raise ValueError("Spin-test null distribution is incomplete or non-finite.")

    null_path = output_dir / "null_distribution_HC_ABetaNeg_lam001.npy"
    np.save(null_path, null_distribution)
    summary_path = output_dir / "EI_HC_ABetaNeg_turbulence_lam001_summary.csv"
    summary = pd.DataFrame(
        [
            {
                "group": turbulence["group"],
                "data_stage": turbulence["data_stage"],
                "sample_size": int(turbulence["sample_size"]),
                "physical_lambda": float(turbulence["physical_lambda"]),
                "legacy_lambda_suffix": turbulence["legacy_lambda_suffix"],
                "metric": metric,
                "empirical_r": float(r_value),
                "p_spin": float(p_spin),
                "n_perm": n_perm,
                "seed": seed,
                "surface_atlas": atlas,
                "surface_density": density,
                "hemisphere": hemisphere,
                "n_source_parcels": source_parcels,
                "n_surface_vertices": int(ei_data.size),
                "null_distribution_file": null_path.name,
            }
        ]
    )
    summary.to_csv(summary_path, index=False)

    figure, axis = plt.subplots(figsize=(6.4, 4.8))
    axis.hist(
        null_distribution,
        bins=35,
        color="#C9D4DF",
        edgecolor="white",
        linewidth=0.5,
    )
    axis.axvline(
        r_value,
        color="#D55E00",
        linewidth=2.2,
        label=f"Empirical r = {r_value:.3f}",
    )
    axis.axvline(0, color="black", linewidth=0.8, alpha=0.55)
    axis.set_xlabel("Pearson's r under spatial rotation")
    axis.set_ylabel("Count")
    axis.set_title("E:I vs HC turbulence at lambda=0.01")
    axis.text(
        0.98,
        0.95,
        f"HC amyloid-negative, N={int(turbulence['sample_size'])}\n"
        f"Alexander-Bloch p={p_spin:.3f} ({n_perm:,} rotations)",
        transform=axis.transAxes,
        ha="right",
        va="top",
        fontsize=9,
    )
    axis.legend(loc="upper left", frameon=False)
    figure.tight_layout()
    figure.savefig(figure_path, dpi=300)
    plt.close(figure)

    metadata = {
        "analysis_name": config["analysis_name"],
        "analysis_definition": (
            "Spatial correlation between the primary E:I map and the parcel-wise "
            "HC amyloid-negative group-mean turbulence map; this is not a group contrast."
        ),
        "configuration": config,
        "effective_settings": {"n_perm": n_perm, "seed": seed},
        "inputs": {
            "ei_map": {"path": str(ei_path), "sha256": sha256(ei_path)},
            "hc_turbulence_map": {
                "path": str(turbulence_path),
                "sha256": sha256(turbulence_path),
            },
        },
        "software": {
            "neuromaps": getattr(neuromaps, "__version__", "unknown"),
            "nibabel": getattr(nibabel, "__version__", "unknown"),
            "numpy": np.__version__,
            "pandas": pd.__version__,
        },
        "result": {"empirical_r": float(r_value), "p_spin": float(p_spin)},
        "outputs": {
            "summary": str(summary_path),
            "null_distribution": str(null_path),
            "figure": str(figure_path),
        },
    }
    metadata_path = output_dir / "EI_HC_ABetaNeg_turbulence_lam001_metadata.json"
    with metadata_path.open("w", encoding="utf-8") as stream:
        json.dump(metadata, stream, indent=2)
        stream.write("\n")

    print(f"HC_ABetaNeg: r={r_value:.6f}, p_spin={p_spin:.6f}")
    print(f"Saved summary: {summary_path}")
    print(f"Saved figure: {figure_path}")
    print(f"Saved analysis metadata: {metadata_path}")


if __name__ == "__main__":
    main()
