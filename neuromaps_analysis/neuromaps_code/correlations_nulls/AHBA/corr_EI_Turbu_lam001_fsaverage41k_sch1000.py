#!/usr/bin/env python3
"""Correlate the primary Schaefer-1000 E:I map with lambda=0.01 turbulence."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import tempfile


NEUROMAPS_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_CONFIG = NEUROMAPS_ROOT / "config" / "ei_turbulence_lam001.json"
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
    parser.add_argument(
        "--contrast",
        action="append",
        help="Contrast to run; repeat for multiple contrasts. Defaults to all configured.",
    )
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
    """Return the requested in-memory GIFTI hemisphere."""
    hemispheres = [
        np.asarray(image.agg_data(), dtype=float).squeeze() for image in surface
    ]
    if hemisphere != "left":
        raise ValueError("This paper-style analysis requires surface.hemisphere='left'.")
    data = hemispheres[0]
    if data.ndim != 1 or not np.isfinite(data).all():
        raise ValueError("Transformed surface map must be a finite one-dimensional array.")
    return data


def main() -> None:
    args = parse_args()
    config_path = args.config.expanduser().resolve()
    with config_path.open(encoding="utf-8") as stream:
        config = json.load(stream)

    configured_maps = config["turbulence_maps"]
    contrasts = args.contrast or list(configured_maps)
    unknown = sorted(set(contrasts) - set(configured_maps))
    if unknown:
        raise ValueError(f"Unknown contrast(s): {unknown}")

    ei_path = resolve_from_config(config_path, config["ei_map"])
    output_dir = resolve_from_config(config_path, config["output_dir"])
    output_dir.mkdir(parents=True, exist_ok=True)
    if not ei_path.is_file():
        raise FileNotFoundError(f"Missing E:I map: {ei_path}")

    atlas = config["surface"]["atlas"]
    density = config["surface"]["density"]
    hemisphere = config["surface"]["hemisphere"]
    source_parcels = int(config["surface"]["source_parcels"])
    if hemisphere != "left" or source_parcels != 500:
        raise ValueError("Statistical correlations must use the 500 left parcels.")
    metric = config["metric"]
    n_perm = (
        args.n_perm
        if args.n_perm is not None
        else int(config["null_model"]["n_perm"])
    )
    seed = args.seed if args.seed is not None else int(config["null_model"]["seed"])
    if n_perm < 1:
        raise ValueError("n_perm must be a positive integer.")

    print(f"Transforming E:I map to {atlas} {density}: {ei_path}")
    ei_surface = transforms.mni152_to_fsaverage(str(ei_path), density)
    ei_data = surface_to_array(ei_surface, hemisphere)
    full_spin_indices = np.asarray(
        nulls.alexander_bloch(
            None, atlas=atlas, density=density, n_perm=n_perm, seed=seed
        )
    )
    if full_spin_indices.shape != (2 * ei_data.size, n_perm):
        raise ValueError(f"Unexpected bilateral spin-index shape: {full_spin_indices.shape}")
    spin_indices = full_spin_indices[:ei_data.size]
    if spin_indices.min() < 0 or spin_indices.max() >= ei_data.size:
        raise ValueError("Left-hemisphere spins reference right-hemisphere vertices.")
    results: list[dict[str, object]] = []
    null_distributions: dict[str, np.ndarray] = {}

    for contrast in contrasts:
        turbulence_path = resolve_from_config(config_path, configured_maps[contrast])
        if not turbulence_path.is_file():
            raise FileNotFoundError(f"Missing turbulence map: {turbulence_path}")

        print(f"Transforming {contrast} turbulence map: {turbulence_path}")
        turbulence_surface = transforms.mni152_to_fsaverage(
            str(turbulence_path), density
        )
        turbulence_data = surface_to_array(turbulence_surface, hemisphere)
        if turbulence_data.shape != ei_data.shape:
            raise ValueError(
                f"Surface array mismatch: turbulence={turbulence_data.shape}, "
                f"E:I={ei_data.shape}"
            )

        print(
            f"Generating {n_perm} Alexander-Bloch rotations "
            f"for {contrast} (seed={seed})"
        )
        rotated = turbulence_data[spin_indices]
        r_value, p_spin, null_distribution = compare_images(
            turbulence_data,
            ei_data,
            nulls=rotated,
            metric=metric,
            return_nulls=True,
        )
        null_distribution = np.asarray(null_distribution, dtype=float)
        null_path = output_dir / f"null_distribution_{contrast}_lam001.npy"
        np.save(null_path, null_distribution)
        null_distributions[contrast] = null_distribution
        results.append(
            {
                "contrast": contrast,
                "physical_lambda": float(config["physical_lambda"]),
                "legacy_lambda_suffix": config["legacy_lambda_suffix"],
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
        )
        print(f"{contrast}: r={r_value:.6f}, p_spin={p_spin:.6f}")

    summary_path = output_dir / "EI_turbulence_lam001_summary.csv"
    pd.DataFrame(results).to_csv(summary_path, index=False)

    figure, axis = plt.subplots(figsize=(7, 5))
    positions = np.arange(1, len(contrasts) + 1)
    axis.boxplot(
        [null_distributions[name] for name in contrasts],
        positions=positions,
        patch_artist=True,
        boxprops={"facecolor": "lightgray"},
    )
    axis.scatter(
        positions,
        [row["empirical_r"] for row in results],
        color="darkorange",
        s=65,
        label="Empirical r",
        zorder=3,
    )
    axis.axhline(0, color="black", linewidth=0.8, alpha=0.5)
    axis.set_xticks(positions, [name.upper().replace("_", "-") for name in contrasts])
    axis.set_ylabel("Pearson's r")
    axis.set_title("E:I-turbulence spatial correlation at lambda=0.01")
    axis.legend()
    figure.tight_layout()
    figure_path = output_dir / "EI_turbulence_lam001_nulls.pdf"
    figure.savefig(figure_path, dpi=300)
    plt.close(figure)

    metadata = {
        "analysis_name": config["analysis_name"],
        "configuration": config,
        "effective_settings": {"contrasts": contrasts, "n_perm": n_perm, "seed": seed},
        "inputs": {
            "ei_map": str(ei_path),
            "ei_map_sha256": sha256(ei_path),
            "turbulence_maps": {
                contrast: {
                    "path": str(resolve_from_config(config_path, configured_maps[contrast])),
                    "sha256": sha256(
                        resolve_from_config(config_path, configured_maps[contrast])
                    ),
                }
                for contrast in contrasts
            },
        },
        "software": {
            "neuromaps": getattr(neuromaps, "__version__", "unknown"),
            "nibabel": getattr(nibabel, "__version__", "unknown"),
            "numpy": np.__version__,
            "pandas": pd.__version__,
        },
        "outputs": {
            "summary": summary_path.name,
            "figure": figure_path.name,
        },
    }
    metadata_path = output_dir / "EI_turbulence_lam001_metadata.json"
    with metadata_path.open("w", encoding="utf-8") as stream:
        json.dump(metadata, stream, indent=2)
        stream.write("\n")

    print(f"Saved summary: {summary_path}")
    print(f"Saved figure: {figure_path}")
    print(f"Saved analysis metadata: {metadata_path}")


if __name__ == "__main__":
    main()
