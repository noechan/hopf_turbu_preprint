#!/usr/bin/env python3
"""Build the HC amyloid-negative mean turbulence map at physical lambda=0.01."""

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
DEFAULT_CONFIG = ROOT / "config" / "ei_hc_turbulence_lam001.json"


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


def load_matlab_v73_turbulence(
    source: Path,
    variable: str,
    n_parcels: int,
    n_lambdas: int,
) -> np.ndarray:
    """Return subjects x parcels x lambdas from a MATLAB v7.3 dataset.

    MATLAB v7.3 reverses array dimensions in HDF5. The source MATLAB array is
    lambda x parcel x subject, so h5py exposes subject x parcel x lambda.
    """
    with h5py.File(source, "r") as mat:
        if variable not in mat:
            raise KeyError(f"{variable!r} is absent from {source}")
        data = np.asarray(mat[variable], dtype=float)

    if data.ndim != 3:
        raise ValueError(f"{variable} must be three-dimensional; got {data.shape}.")
    if data.shape[1:] != (n_parcels, n_lambdas):
        raise ValueError(
            "Unexpected MATLAB v7.3 storage shape. Expected "
            f"(subjects, {n_parcels}, {n_lambdas}); got {data.shape}."
        )
    if not np.isfinite(data).all():
        raise ValueError(f"{variable} contains non-finite values.")
    return data


def main() -> None:
    args = parse_args()
    config_path = args.config.expanduser().resolve()
    with config_path.open(encoding="utf-8") as stream:
        config = json.load(stream)

    turbulence = config["turbulence"]
    source_override = os.environ.get(turbulence["source_env_var"])
    source = Path(source_override).expanduser().resolve() if source_override else Path(
        turbulence["source_mat"]
    ).expanduser().resolve()
    atlas_path = resolve_from_config(config_path, config["atlas_map"])
    csv_path = resolve_from_config(config_path, turbulence["group_mean_csv"])
    map_path = resolve_from_config(config_path, turbulence["group_mean_map"])
    metadata_path = resolve_from_config(config_path, turbulence["build_metadata"])

    if not source.is_file():
        raise FileNotFoundError(
            f"Missing HC turbulence source: {source}. Set "
            f"{turbulence['source_env_var']} to override it."
        )
    if not atlas_path.is_file():
        raise FileNotFoundError(f"Missing Schaefer atlas: {atlas_path}")

    lambda_values = np.asarray(turbulence["lambda_values"], dtype=float)
    matlab_index = int(turbulence["matlab_lambda_index"])
    python_index = matlab_index - 1
    physical_lambda = float(turbulence["physical_lambda"])
    if not 0 <= python_index < lambda_values.size:
        raise ValueError("matlab_lambda_index is outside lambda_values.")
    if not np.isclose(lambda_values[python_index], physical_lambda):
        raise ValueError(
            f"Lambda index {matlab_index} maps to {lambda_values[python_index]}, "
            f"not {physical_lambda}."
        )

    atlas_img = nib.load(atlas_path)
    atlas = np.asarray(atlas_img.dataobj)
    labels = np.unique(atlas[atlas > 0]).astype(int)
    n_parcels = labels.size
    if n_parcels != 1000 or not np.array_equal(labels, np.arange(1, 1001)):
        raise ValueError("Atlas must contain integer labels 1 through 1000.")

    data = load_matlab_v73_turbulence(
        source,
        turbulence["source_variable"],
        n_parcels=n_parcels,
        n_lambdas=lambda_values.size,
    )
    n_subjects = data.shape[0]
    expected_n = int(turbulence["sample_size"])
    if n_subjects != expected_n:
        raise ValueError(
            f"Source contains {n_subjects} subjects; configuration expects {expected_n}."
        )

    subject_by_parcel = data[:, :, python_index]
    parcel_mean = subject_by_parcel.mean(axis=0)
    if parcel_mean.shape != (n_parcels,) or not np.isfinite(parcel_mean).all():
        raise ValueError("HC group mean must contain 1,000 finite parcel values.")

    csv_path.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(
        {
            "parcel_id": np.arange(1, n_parcels + 1),
            "parcel_name": [f"Schaefer_{i}" for i in range(1, n_parcels + 1)],
            "turbulence_mean": parcel_mean,
        }
    ).to_csv(csv_path, index=False)

    mapped = np.zeros(atlas.shape, dtype=np.float32)
    for label, value in enumerate(parcel_mean, start=1):
        mapped[atlas == label] = value
    header = atlas_img.header.copy()
    header.set_data_dtype(np.float32)
    nib.save(nib.Nifti1Image(mapped, atlas_img.affine, header), map_path)

    metadata = {
        "map_definition": "parcel-wise mean turbulence across HC amyloid-negative participants",
        "group": turbulence["group"],
        "data_stage": turbulence["data_stage"],
        "sample_size": n_subjects,
        "physical_lambda": physical_lambda,
        "matlab_lambda_index": matlab_index,
        "legacy_lambda_suffix": turbulence["legacy_lambda_suffix"],
        "source": {
            "path": str(source),
            "variable": turbulence["source_variable"],
            "hdf5_shape": list(data.shape),
            "sha256": sha256(source),
        },
        "atlas": {"path": str(atlas_path), "sha256": sha256(atlas_path)},
        "parcel_value_range": [float(parcel_mean.min()), float(parcel_mean.max())],
        "outputs": {"csv": str(csv_path), "nifti": str(map_path)},
    }
    with metadata_path.open("w", encoding="utf-8") as stream:
        json.dump(metadata, stream, indent=2)
        stream.write("\n")

    print(f"HC subjects: {n_subjects}")
    print(f"Physical lambda: {physical_lambda} (MATLAB index {matlab_index})")
    print(f"Parcel range: {parcel_mean.min():.8f} to {parcel_mean.max():.8f}")
    print(f"Saved parcel table: {csv_path}")
    print(f"Saved NIfTI map: {map_path}")
    print(f"Saved build metadata: {metadata_path}")


if __name__ == "__main__":
    main()
