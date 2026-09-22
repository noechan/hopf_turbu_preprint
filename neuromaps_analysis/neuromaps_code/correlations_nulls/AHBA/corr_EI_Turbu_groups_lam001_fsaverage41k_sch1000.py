#!/usr/bin/env python3
"""Compare E:I correlations across raw group-mean turbulence maps."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import tempfile


NEUROMAPS_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_CONFIG = NEUROMAPS_ROOT / "config" / "ei_group_turbulence_lam001.json"
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


def group_map_path(map_dir: Path, group: str, sample_size: int) -> Path:
    return map_dir / f"turbu_mean_{group}_lam001_sch1000_N{sample_size}_2mm.nii.gz"


def main() -> None:
    args = parse_args()
    config_path = args.config.expanduser().resolve()
    with config_path.open(encoding="utf-8") as stream:
        config = json.load(stream)

    ei_path = resolve_from_config(config_path, config["ei_map"])
    map_dir = resolve_from_config(config_path, config["map_output_dir"])
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
    if n_perm < 1:
        raise ValueError("n_perm must be a positive integer.")

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
    for group, settings in config["groups"].items():
        sample_size = int(settings["sample_size"])
        map_path = group_map_path(map_dir, group, sample_size)
        if not map_path.is_file():
            raise FileNotFoundError(f"Missing {group} turbulence map: {map_path}")
        print(f"Transforming {group}: {map_path}")
        turbulence_data = surface_to_array(
            transforms.mni152_to_fsaverage(str(map_path), density), hemisphere
        )
        if turbulence_data.shape != ei_data.shape:
            raise ValueError(
                f"Surface mismatch for {group}: {turbulence_data.shape} vs {ei_data.shape}"
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
        if null_distribution.shape != (n_perm,) or not np.isfinite(null_distribution).all():
            raise ValueError(f"Invalid null distribution for {group}.")
        null_path = output_dir / f"null_distribution_{group}_lam001.npy"
        np.save(null_path, null_distribution)
        null_distributions.append(null_distribution)
        input_metadata[group] = {"path": str(map_path), "sha256": sha256(map_path)}
        results.append(
            {
                "group": group,
                "display_name": settings["display_name"],
                "data_stage": config["data_stage"],
                "sample_size": sample_size,
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
        print(f"{group}: r={r_value:.6f}, p_spin={p_spin:.6f}")
        del rotated

    p_values = np.asarray([row["p_spin"] for row in results], dtype=float)
    rejected, p_fdr, _, _ = multipletests(p_values, alpha=0.05, method=config["multiple_testing"])
    reference_r = float(results[0]["empirical_r"])
    for row, reject, adjusted in zip(results, rejected, p_fdr):
        row["delta_r_vs_HC_ABetaNeg"] = float(row["empirical_r"]) - reference_r
        row["p_fdr_bh"] = float(adjusted)
        row["significant_fdr_0.05"] = bool(reject)

    summary_path = output_dir / "EI_group_turbulence_lam001_summary.csv"
    pd.DataFrame(results).to_csv(summary_path, index=False)

    figure, axis = plt.subplots(figsize=(8.2, 5.2))
    positions = np.arange(1, len(results) + 1)
    axis.boxplot(
        null_distributions,
        positions=positions,
        widths=0.58,
        showfliers=False,
        patch_artist=True,
        boxprops={"facecolor": "#D4DEE8", "edgecolor": "#667788"},
        medianprops={"color": "#334455"},
        whiskerprops={"color": "#667788"},
        capprops={"color": "#667788"},
    )
    empirical = [float(row["empirical_r"]) for row in results]
    axis.plot(positions, empirical, color="#D55E00", linewidth=1.6, alpha=0.75)
    axis.scatter(
        positions,
        empirical,
        color="#D55E00",
        edgecolor="white",
        linewidth=0.8,
        s=70,
        label="Empirical r",
        zorder=3,
    )
    axis.axhline(0, color="black", linewidth=0.8, alpha=0.55)
    axis.set_ylim(-0.32, 0.47)
    axis.set_xticks(
        positions,
        [
            "HC amyloid-\nnegative",
            "HC amyloid-\npositive",
            "MCI amyloid-\npositive",
            "AD amyloid-\npositive",
        ],
    )
    axis.set_ylabel("Pearson's r")
    axis.set_title(
        "E:I correlation with group-mean turbulence at lambda=0.01", pad=12
    )
    for x, row in zip(positions, results):
        axis.annotate(
            f"r={float(row['empirical_r']):.3f}\np={float(row['p_spin']):.3f}",
            (x, float(row["empirical_r"])),
            xytext=(0, 10),
            textcoords="offset points",
            ha="center",
            fontsize=8.5,
        )
    axis.legend(loc="lower left", frameon=False)
    figure.tight_layout()
    figure.savefig(figure_path, dpi=300)
    plt.close(figure)

    metadata = {
        "analysis_name": config["analysis_name"],
        "analysis_definition": (
            "Separate spatial correlations between the primary E:I map and each "
            "raw group-mean turbulence map; these are not group contrasts."
        ),
        "configuration": config,
        "effective_settings": {"n_perm": n_perm, "seed": seed},
        "inputs": {
            "ei_map": {"path": str(ei_path), "sha256": sha256(ei_path)},
            "group_maps": input_metadata,
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
            "figure": str(figure_path),
        },
    }
    metadata_path = output_dir / "EI_group_turbulence_lam001_metadata.json"
    with metadata_path.open("w", encoding="utf-8") as stream:
        json.dump(metadata, stream, indent=2)
        stream.write("\n")

    print(f"Saved summary: {summary_path}")
    print(f"Saved figure: {figure_path}")
    print(f"Saved metadata: {metadata_path}")


if __name__ == "__main__":
    main()
