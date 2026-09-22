#!/usr/bin/env python3
"""Build raw group-mean Schaefer-1000 turbulence maps at lambda=0.01."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path

import h5py
import nibabel as nib
import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_CONFIG = ROOT / "config" / "ei_group_turbulence_lam001.json"


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


def main() -> None:
    args = parse_args()
    config_path = args.config.expanduser().resolve()
    with config_path.open(encoding="utf-8") as stream:
        config = json.load(stream)

    atlas_path = resolve_from_config(config_path, config["atlas_map"])
    output_dir = resolve_from_config(config_path, config["map_output_dir"])
    output_dir.mkdir(parents=True, exist_ok=True)
    if not atlas_path.is_file():
        raise FileNotFoundError(f"Missing Schaefer atlas: {atlas_path}")

    atlas_img = nib.load(atlas_path)
    atlas = np.asarray(atlas_img.dataobj)
    labels = np.unique(atlas[atlas > 0]).astype(int)
    if labels.size != 1000 or not np.array_equal(labels, np.arange(1, 1001)):
        raise ValueError("Atlas must contain integer labels 1 through 1000.")

    lambda_values = np.asarray(config["lambda_values"], dtype=float)
    matlab_index = int(config["matlab_lambda_index"])
    python_index = matlab_index - 1
    if not np.isclose(lambda_values[python_index], float(config["physical_lambda"])):
        raise ValueError("Configured MATLAB index does not match physical_lambda.")

    manifest: list[dict[str, object]] = []
    for group, settings in config["groups"].items():
        override = os.environ.get(settings["source_env_var"])
        source = Path(override or settings["source_mat"]).expanduser().resolve()
        if not source.is_file():
            raise FileNotFoundError(
                f"Missing {group} source: {source}. Set "
                f"{settings['source_env_var']} to override it."
            )

        with h5py.File(source, "r") as mat:
            variable = config["source_variable"]
            if variable not in mat:
                raise KeyError(f"{variable!r} is absent from {source}")
            data = np.asarray(mat[variable], dtype=float)

        expected_shape = (int(settings["sample_size"]), 1000, lambda_values.size)
        if data.shape != expected_shape:
            raise ValueError(
                f"{group} has HDF5 shape {data.shape}; expected {expected_shape}."
            )
        if not np.isfinite(data).all():
            raise ValueError(f"{group} contains non-finite turbulence values.")

        values = data[:, :, python_index].mean(axis=0)
        csv_path = output_dir / f"turbu_mean_{group}_lam001_sch1000_N{data.shape[0]}.csv"
        map_path = output_dir / (
            f"turbu_mean_{group}_lam001_sch1000_N{data.shape[0]}_2mm.nii.gz"
        )
        pd.DataFrame(
            {
                "parcel_id": np.arange(1, 1001),
                "parcel_name": [f"Schaefer_{i}" for i in range(1, 1001)],
                "turbulence_mean": values,
            }
        ).to_csv(csv_path, index=False)

        mapped = np.zeros(atlas.shape, dtype=np.float32)
        for label, value in enumerate(values, start=1):
            mapped[atlas == label] = value
        header = atlas_img.header.copy()
        header.set_data_dtype(np.float32)
        nib.save(nib.Nifti1Image(mapped, atlas_img.affine, header), map_path)

        manifest.append(
            {
                "group": group,
                "display_name": settings["display_name"],
                "sample_size": data.shape[0],
                "data_stage": config["data_stage"],
                "physical_lambda": float(config["physical_lambda"]),
                "matlab_lambda_index": matlab_index,
                "source_path": str(source),
                "source_sha256": sha256(source),
                "source_hdf5_shape": list(data.shape),
                "parcel_value_min": float(values.min()),
                "parcel_value_max": float(values.max()),
                "csv": str(csv_path),
                "nifti": str(map_path),
                "nifti_sha256": sha256(map_path),
            }
        )
        print(
            f"{group}: N={data.shape[0]}, range={values.min():.8f} "
            f"to {values.max():.8f}"
        )

    manifest_path = output_dir / "turbu_group_means_lam001_build_metadata.json"
    with manifest_path.open("w", encoding="utf-8") as stream:
        json.dump(
            {
                "map_definition": "raw parcel-wise group-mean turbulence",
                "configuration": str(config_path),
                "atlas": {"path": str(atlas_path), "sha256": sha256(atlas_path)},
                "maps": manifest,
            },
            stream,
            indent=2,
        )
        stream.write("\n")
    print(f"Saved build metadata: {manifest_path}")


if __name__ == "__main__":
    main()
