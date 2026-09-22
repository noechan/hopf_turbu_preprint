"""Load and validate the inputs for the participant-level analysis.

This module does not perform harmonization. It reads the canonical N145
Schaefer-1000 table produced by ``run_harmonization.py`` and confirms that it
matches the manuscript cohort and the retained HC-minus-AD group map.
"""

from __future__ import annotations

from dataclasses import dataclass
import hashlib
import json
import os
from pathlib import Path

import nibabel as nib
import numpy as np
import pandas as pd
from scipy.io import loadmat


@dataclass(frozen=True)
class InputPaths:
    """Resolved paths recorded in the analysis provenance."""

    harmonized_node_table: Path
    harmonization_script: Path
    demographics: Path
    n145_cohort_table: Path
    ei_values: Path
    ei_nifti: Path
    atlas_map: Path
    validation_reference_hc_ad: Path


@dataclass
class AnalysisInputs:
    """Validated data needed by the spatial and statistical stages."""

    harmonized_nodes: pd.DataFrame
    participants: pd.DataFrame
    parcel_columns: list[str]
    atlas_img: nib.Nifti1Image
    atlas_volume: np.ndarray
    paths: InputPaths
    validation_max_error: float


def read_config(config_path: Path) -> dict:
    """Read the JSON analysis configuration."""

    with config_path.open(encoding="utf-8") as stream:
        return json.load(stream)


def resolve_from_config(config_path: Path, value: str) -> Path:
    """Resolve a path relative to the directory containing the configuration."""

    path = Path(value).expanduser()
    if not path.is_absolute():
        path = config_path.parent / path
    return path.resolve()


def resolve_external(config_path: Path, section: dict[str, str]) -> Path:
    """Resolve a configurable path, giving its environment variable priority."""

    override = os.environ.get(section["env_var"])
    if override:
        return Path(override).expanduser().resolve()
    return resolve_from_config(config_path, section["path"])


def file_sha256(path: Path) -> str:
    """Return the SHA-256 checksum used in the provenance record."""

    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def resolve_input_paths(config_path: Path, config: dict) -> InputPaths:
    """Resolve all configured input paths and confirm that they exist."""

    paths = InputPaths(
        harmonized_node_table=resolve_external(
            config_path, config["harmonized_node_table"]
        ),
        harmonization_script=resolve_from_config(
            config_path, config["harmonization_script"]
        ),
        demographics=resolve_external(config_path, config["demographics"]),
        n145_cohort_table=resolve_from_config(
            config_path, config["n145_cohort_table"]
        ),
        ei_values=resolve_from_config(config_path, config["ei_values"]),
        ei_nifti=resolve_from_config(config_path, config["ei_nifti"]),
        atlas_map=resolve_from_config(config_path, config["atlas_map"]),
        validation_reference_hc_ad=resolve_from_config(
            config_path, config["validation_reference_hc_ad"]
        ),
    )
    missing = [str(path) for path in vars(paths).values() if not path.is_file()]
    if missing:
        raise FileNotFoundError("Missing required inputs:\n" + "\n".join(missing))
    return paths


def validate_group_counts(
    frame: pd.DataFrame, expected: dict[str, int], label: str
) -> None:
    counts = frame["Group"].value_counts().to_dict()
    if counts != expected:
        raise ValueError(f"Unexpected {label} group counts: {counts}; expected {expected}")


def _validate_cohort(
    harmonized: pd.DataFrame,
    cohort: pd.DataFrame,
    expected_counts: dict[str, int],
) -> None:
    """Confirm that the harmonized table is exactly the manuscript N145 cohort."""

    for frame in (harmonized, cohort):
        frame["PTID"] = frame["PTID"].astype(str)
        if frame["PTID"].duplicated().any():
            raise ValueError("A participant input contains duplicate PTIDs.")

    validate_group_counts(harmonized, expected_counts, "harmonized N145")
    validate_group_counts(cohort, expected_counts, "N145 cohort")
    harmonized_labels = harmonized.set_index("PTID")["Group"].sort_index()
    cohort_labels = cohort.set_index("PTID")["Group"].sort_index()
    if harmonized_labels.equals(cohort_labels):
        return

    missing = sorted(set(cohort_labels.index) - set(harmonized_labels.index))
    unexpected = sorted(set(harmonized_labels.index) - set(cohort_labels.index))
    raise ValueError(
        "The canonical harmonized node table does not match the manuscript "
        "N145 cohort and group labels. "
        f"Missing PTIDs: {missing}; unexpected PTIDs: {unexpected}."
    )

def _validate_parcels(harmonized: pd.DataFrame) -> list[str]:
    """Check the expected ordered Schaefer-1000 columns and finite values."""

    parcel_columns = [f"Schaefer_{index}" for index in range(1, 1001)]
    observed = [
        column for column in harmonized.columns if column.startswith("Schaefer_")
    ]
    if observed != parcel_columns:
        raise ValueError(
            "Node table must contain ordered Schaefer_1 through Schaefer_1000."
        )
    matrix = harmonized[parcel_columns].to_numpy(dtype=float)
    if not np.isfinite(matrix).all():
        raise ValueError("Harmonized node table contains non-finite values.")
    return parcel_columns


def _validate_retained_group_map(
    harmonized: pd.DataFrame,
    parcel_columns: list[str],
    reference_path: Path,
) -> float:
    """Compare the HC-minus-AD mean map with the retained manuscript map."""

    reference = np.asarray(
        loadmat(reference_path)["diff_hc_ad_lam1"], dtype=float
    ).squeeze()
    observed = (
        harmonized.loc[harmonized["Group"] == "HC_ABneg", parcel_columns]
        .mean(axis=0)
        .to_numpy()
        - harmonized.loc[harmonized["Group"] == "AD_ABpos", parcel_columns]
        .mean(axis=0)
        .to_numpy()
    )
    max_error = float(np.max(np.abs(reference - observed)))
    if max_error > 1e-10:
        raise ValueError(
            "The canonical harmonized values do not reproduce the retained "
            f"HC-AD aggregate map (max absolute error={max_error:.3g})."
        )
    return max_error


def load_and_validate_inputs(config_path: Path, config: dict) -> AnalysisInputs:
    """Load every analysis input and assemble the participant metadata table."""

    paths = resolve_input_paths(config_path, config)
    print(
        "Loading canonical harmonized lambda=0.01 node table: "
        f"{paths.harmonized_node_table}"
    )
    harmonized = pd.read_excel(paths.harmonized_node_table, engine="openpyxl")
    cohort = pd.read_excel(
        paths.n145_cohort_table,
        usecols=["PTID", "Group"],
        engine="openpyxl",
    )
    _validate_cohort(harmonized, cohort, config["expected_final_group_counts"])
    parcel_columns = _validate_parcels(harmonized)
    max_error = _validate_retained_group_map(
        harmonized, parcel_columns, paths.validation_reference_hc_ad
    )
    print(f"Harmonized aggregate validation max error: {max_error:.3g}")

    demographics = pd.read_csv(paths.demographics).assign(
        PTID=lambda frame: frame["PTID"].astype(str)
    )
    if demographics["PTID"].duplicated().any():
        raise ValueError("The demographics input contains duplicate PTIDs.")
    participants = harmonized[["PTID", "Group"]].merge(
        demographics[["PTID", "age", "gender", "edu"]],
        on="PTID",
        how="left",
    )
    group_order = list(config["group_order"])
    participants["stage"] = participants["Group"].map(
        {group: index for index, group in enumerate(group_order)}
    )
    if participants.isna().any().any():
        raise ValueError("Final participant metadata contains missing values.")

    atlas_img = nib.load(paths.atlas_map)
    atlas_volume = np.asarray(atlas_img.dataobj)
    atlas_labels = np.unique(atlas_volume[atlas_volume > 0]).astype(int)
    if not np.array_equal(atlas_labels, np.arange(1, 1001)):
        raise ValueError("Atlas must contain labels 1 through 1000.")

    return AnalysisInputs(
        harmonized_nodes=harmonized,
        participants=participants,
        parcel_columns=parcel_columns,
        atlas_img=atlas_img,
        atlas_volume=atlas_volume,
        paths=paths,
        validation_max_error=max_error,
    )
