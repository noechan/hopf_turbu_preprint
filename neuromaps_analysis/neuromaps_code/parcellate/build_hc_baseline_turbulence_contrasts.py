#!/usr/bin/env python3
"""Build signed turbulence contrasts relative to HC amyloid-negative controls."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import nibabel as nib
import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_CONFIG = ROOT / "config" / "ei_hc_baseline_contrasts_lam001.json"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
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


def group_csv_path(group_dir: Path, group: str, sample_size: int) -> Path:
    return group_dir / f"turbu_mean_{group}_lam001_sch1000_N{sample_size}.csv"


def load_group_values(path: Path) -> np.ndarray:
    if not path.is_file():
        raise FileNotFoundError(f"Missing group-mean table: {path}")
    table = pd.read_csv(path)
    if list(table.columns) != ["parcel_id", "parcel_name", "turbulence_mean"]:
        raise ValueError(f"Unexpected columns in {path}.")
    if not np.array_equal(table["parcel_id"].to_numpy(), np.arange(1, 1001)):
        raise ValueError(f"Parcel order is invalid in {path}.")
    values = table["turbulence_mean"].to_numpy(dtype=float)
    if values.shape != (1000,) or not np.isfinite(values).all():
        raise ValueError(f"Expected 1,000 finite turbulence values in {path}.")
    return values


def main() -> None:
    args = parse_args()
    config_path = args.config.expanduser().resolve()
    with config_path.open(encoding="utf-8") as stream:
        config = json.load(stream)

    atlas_path = resolve_from_config(config_path, config["atlas_map"])
    group_dir = resolve_from_config(config_path, config["group_map_dir"])
    output_dir = resolve_from_config(config_path, config["contrast_output_dir"])
    output_dir.mkdir(parents=True, exist_ok=True)
    if not atlas_path.is_file():
        raise FileNotFoundError(f"Missing Schaefer atlas: {atlas_path}")

    atlas_img = nib.load(atlas_path)
    atlas = np.asarray(atlas_img.dataobj)
    labels = np.unique(atlas[atlas > 0]).astype(int)
    if labels.size != 1000 or not np.array_equal(labels, np.arange(1, 1001)):
        raise ValueError("Atlas must contain integer labels 1 through 1000.")

    baseline_group = config["baseline_group"]
    baseline_settings = config["groups"][baseline_group]
    baseline_csv = group_csv_path(
        group_dir, baseline_group, int(baseline_settings["sample_size"])
    )
    baseline = load_group_values(baseline_csv)

    metadata_rows: list[dict[str, object]] = []
    for contrast, settings in config["contrasts"].items():
        comparison_group = settings["comparison_group"]
        comparison_settings = config["groups"][comparison_group]
        comparison_csv = group_csv_path(
            group_dir, comparison_group, int(comparison_settings["sample_size"])
        )
        comparison = load_group_values(comparison_csv)
        values = baseline - comparison

        csv_path = output_dir / f"turbu_{contrast}_lam001_sch1000_raw.csv"
        map_path = output_dir / f"turbu_{contrast}_lam001_sch1000_raw_2mm.nii.gz"
        pd.DataFrame(
            {
                "parcel_id": np.arange(1, 1001),
                "parcel_name": [f"Schaefer_{i}" for i in range(1, 1001)],
                "turbulence_difference": values,
            }
        ).to_csv(csv_path, index=False)

        mapped = np.zeros(atlas.shape, dtype=np.float32)
        for label, value in enumerate(values, start=1):
            mapped[atlas == label] = value
        header = atlas_img.header.copy()
        header.set_data_dtype(np.float32)
        nib.save(nib.Nifti1Image(mapped, atlas_img.affine, header), map_path)

        metadata_rows.append(
            {
                "contrast": contrast,
                "display_name": settings["display_name"],
                "operation": f"{baseline_group} minus {comparison_group}",
                "baseline_sample_size": int(baseline_settings["sample_size"]),
                "comparison_sample_size": int(comparison_settings["sample_size"]),
                "physical_lambda": float(config["physical_lambda"]),
                "data_stage": config["data_stage"],
                "baseline_csv": {"path": str(baseline_csv), "sha256": sha256(baseline_csv)},
                "comparison_csv": {
                    "path": str(comparison_csv),
                    "sha256": sha256(comparison_csv),
                },
                "difference_range": [float(values.min()), float(values.max())],
                "csv": str(csv_path),
                "nifti": str(map_path),
                "nifti_sha256": sha256(map_path),
            }
        )
        print(
            f"{contrast}: range={values.min():.8f} to {values.max():.8f}, "
            f"mean={values.mean():.8f}"
        )

    metadata_path = output_dir / "turbu_HC_ABetaNeg_baseline_contrasts_metadata.json"
    with metadata_path.open("w", encoding="utf-8") as stream:
        json.dump(
            {
                "analysis_definition": config["analysis_definition"],
                "configuration": str(config_path),
                "atlas": {"path": str(atlas_path), "sha256": sha256(atlas_path)},
                "contrasts": metadata_rows,
            },
            stream,
            indent=2,
        )
        stream.write("\n")
    print(f"Saved contrast metadata: {metadata_path}")


if __name__ == "__main__":
    main()
