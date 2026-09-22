#!/usr/bin/env python3
"""Correlate E:I with turbulence contrasts relative to HC amyloid-negative."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import tempfile


NEUROMAPS_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_CONFIG = NEUROMAPS_ROOT / "config" / "ei_hc_baseline_contrasts_lam001.json"
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
from statsmodels.stats.multitest import multipletests


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


def contrast_map_path(contrast_dir: Path, contrast: str) -> Path:
    return contrast_dir / f"turbu_{contrast}_lam001_sch1000_raw_2mm.nii.gz"


def main() -> None:
    args = parse_args()
    config_path = args.config.expanduser().resolve()
    with config_path.open(encoding="utf-8") as stream:
        config = json.load(stream)

    ei_path = resolve_from_config(config_path, config["ei_map"])
    contrast_dir = resolve_from_config(config_path, config["contrast_output_dir"])
    output_dir = resolve_from_config(config_path, config["output_dir"])
    figure_path = resolve_from_config(config_path, config["figure_file"])
    output_dir.mkdir(parents=True, exist_ok=True)
    figure_path.parent.mkdir(parents=True, exist_ok=True)
    if not ei_path.is_file():
        raise FileNotFoundError(f"Missing E:I map: {ei_path}")

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

    print(f"Transforming E:I map to {atlas} {density}: {ei_path}")
    ei_data = surface_to_array(
        transforms.mni152_to_fsaverage(str(ei_path), density), hemisphere
    )
    print(f"Generating one shared set of {n_perm} rotations (seed={seed})")
    full_spin_indices = np.asarray(
        nulls.alexander_bloch(
            None,
            atlas=atlas,
            density=density,
            n_perm=n_perm,
            seed=seed,
        )
    )
    if full_spin_indices.shape != (2 * ei_data.size, n_perm):
        raise ValueError(f"Unexpected bilateral spin-index shape: {full_spin_indices.shape}")
    spin_indices = full_spin_indices[:ei_data.size]
    if spin_indices.min() < 0 or spin_indices.max() >= ei_data.size:
        raise ValueError("Left-hemisphere spins reference right-hemisphere vertices.")

    results: list[dict[str, object]] = []
    null_distributions: list[np.ndarray] = []
    input_metadata: dict[str, dict[str, str]] = {}
    for contrast, settings in config["contrasts"].items():
        map_path = contrast_map_path(contrast_dir, contrast)
        if not map_path.is_file():
            raise FileNotFoundError(f"Missing contrast map: {map_path}")
        print(f"Transforming {contrast}: {map_path}")
        contrast_data = surface_to_array(
            transforms.mni152_to_fsaverage(str(map_path), density), hemisphere
        )
        if contrast_data.shape != ei_data.shape:
            raise ValueError(
                f"Surface mismatch for {contrast}: {contrast_data.shape} vs {ei_data.shape}"
            )

        rotated = contrast_data[spin_indices]
        r_value, p_spin, null_distribution = compare_images(
            contrast_data,
            ei_data,
            nulls=rotated,
            metric=metric,
            return_nulls=True,
        )
        null_distribution = np.asarray(null_distribution, dtype=float)
        if null_distribution.shape != (n_perm,) or not np.isfinite(null_distribution).all():
            raise ValueError(f"Invalid null distribution for {contrast}.")
        null_path = output_dir / f"null_distribution_{contrast}_lam001.npy"
        np.save(null_path, null_distribution)
        null_distributions.append(null_distribution)
        input_metadata[contrast] = {"path": str(map_path), "sha256": sha256(map_path)}
        comparison_group = settings["comparison_group"]
        results.append(
            {
                "contrast": contrast,
                "display_name": settings["display_name"],
                "operation": f"{config['baseline_group']} minus {comparison_group}",
                "baseline_group": config["baseline_group"],
                "comparison_group": comparison_group,
                "baseline_sample_size": int(
                    config["groups"][config["baseline_group"]]["sample_size"]
                ),
                "comparison_sample_size": int(
                    config["groups"][comparison_group]["sample_size"]
                ),
                "data_stage": config["data_stage"],
                "physical_lambda": float(config["physical_lambda"]),
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
        del rotated

    p_values = np.asarray([row["p_spin"] for row in results], dtype=float)
    rejected, p_fdr, _, _ = multipletests(
        p_values, alpha=0.05, method=config["multiple_testing"]
    )
    for row, reject, adjusted in zip(results, rejected, p_fdr):
        row["p_fdr_bh"] = float(adjusted)
        row["significant_fdr_0.05"] = bool(reject)

    summary_path = output_dir / "EI_HC_ABetaNeg_baseline_contrasts_lam001_summary.csv"
    pd.DataFrame(results).to_csv(summary_path, index=False)

    # Match the visual language of the GE15 gene-expression null plots:
    # light-grey null boxes with the empirical correlations overlaid in orange.
    figure, axis = plt.subplots(figsize=(10, 6))
    positions = np.arange(1, len(results) + 1)
    axis.boxplot(
        null_distributions,
        positions=positions,
        patch_artist=True,
        boxprops={"facecolor": "lightgray"},
    )
    empirical = [float(row["empirical_r"]) for row in results]
    axis.scatter(
        positions,
        empirical,
        color="orange",
        marker="o",
        s=60,
        label="Empirical r",
        zorder=3,
    )
    all_values = np.concatenate([*null_distributions, np.asarray(empirical)])
    axis.set_ylim(float(all_values.min()) - 0.06, float(all_values.max()) + 0.14)
    axis.set_xticks(
        positions,
        ["HC- minus HC+", "HC- minus MCI+", "HC- minus AD+"],
    )
    axis.set_ylabel("Pearson's r")
    axis.set_title(
        "Correlation between E:I expression and HC- baseline turbulence "
        "contrasts (fsaverage 41k, lambda=0.01)"
    )
    for x, row in zip(positions, results):
        axis.annotate(
            f"r={float(row['empirical_r']):.3f}\n"
            f"p={float(row['p_spin']):.3f}",
            (x, float(row["empirical_r"])),
            xytext=(0, 10),
            textcoords="offset points",
            ha="center",
            fontsize=8.5,
        )
    axis.legend(loc="best", frameon=False)
    figure.tight_layout()
    figure.savefig(figure_path, dpi=300)
    figure_png_path = figure_path.with_suffix(".png")
    figure.savefig(figure_png_path, dpi=300)
    plt.close(figure)

    metadata = {
        "analysis_name": config["analysis_name"],
        "analysis_definition": config["analysis_definition"],
        "configuration": config,
        "effective_settings": {"n_perm": n_perm, "seed": seed},
        "inputs": {
            "ei_map": {"path": str(ei_path), "sha256": sha256(ei_path)},
            "contrast_maps": input_metadata,
        },
        "software": {
            "neuromaps": getattr(neuromaps, "__version__", "unknown"),
            "nibabel": getattr(nibabel, "__version__", "unknown"),
            "numpy": np.__version__,
            "pandas": pd.__version__,
        },
        "results": results,
        "outputs": {
            "summary": str(summary_path),
            "figure_pdf": str(figure_path),
            "figure_png": str(figure_png_path),
        },
    }
    metadata_path = output_dir / "EI_HC_ABetaNeg_baseline_contrasts_lam001_metadata.json"
    with metadata_path.open("w", encoding="utf-8") as stream:
        json.dump(metadata, stream, indent=2)
        stream.write("\n")

    print(f"Saved summary: {summary_path}")
    print(f"Saved figure PDF: {figure_path}")
    print(f"Saved figure PNG: {figure_png_path}")
    print(f"Saved metadata: {metadata_path}")


if __name__ == "__main__":
    main()
