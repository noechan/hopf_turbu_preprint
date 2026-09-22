#!/usr/bin/env python3
"""Build bilateral Schaefer-1000 and left-only transcriptomic E:I maps.

The AHBA expression matrix is rebuilt with abagen. Expression is estimated
independently for all 1,000 parcels (no hemisphere reflection), preserving the
existing bilateral workflow. A separate 500-parcel left-hemisphere table is
written for the paper-style statistical correlations.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import nibabel as nib
import numpy as np
import pandas as pd


ABAGEN_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_CONFIG = ABAGEN_ROOT / "config" / "ei_sch1000.json"


def enable_abagen_pandas_compatibility() -> None:
    """Supply pandas APIs used by abagen 0.1.3 and removed in pandas 2.x."""
    if not hasattr(pd.DataFrame, "append"):
        def dataframe_append(
            frame: pd.DataFrame,
            other: pd.DataFrame,
            ignore_index: bool = False,
            verify_integrity: bool = False,
            sort: bool = False,
        ) -> pd.DataFrame:
            return pd.concat(
                [frame, other],
                ignore_index=ignore_index,
                verify_integrity=verify_integrity,
                sort=sort,
            )

        pd.DataFrame.append = dataframe_append  # type: ignore[attr-defined]

    if not hasattr(pd.Series, "append"):
        def series_append(
            series: pd.Series,
            other: pd.Series,
            ignore_index: bool = False,
            verify_integrity: bool = False,
        ) -> pd.Series:
            return pd.concat(
                [series, other],
                ignore_index=ignore_index,
                verify_integrity=verify_integrity,
            )

        pd.Series.append = series_append  # type: ignore[attr-defined]

    original_set_axis = pd.DataFrame.set_axis
    if "inplace" not in original_set_axis.__code__.co_varnames:
        def set_axis_compat(
            frame: pd.DataFrame,
            labels: object,
            axis: object = 0,
            copy: bool | None = None,
            inplace: bool = False,
        ) -> pd.DataFrame | None:
            result = original_set_axis(frame, labels, axis=axis, copy=copy)
            if inplace:
                frame.axes[frame._get_axis_number(axis)] = result.axes[
                    result._get_axis_number(axis)
                ]
                return None
            return result

        pd.DataFrame.set_axis = set_axis_compat  # type: ignore[method-assign]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    parser.add_argument(
        "--force",
        action="store_true",
        help="Overwrite existing E:I outputs.",
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


def ensure_available(outputs: list[Path], force: bool) -> None:
    existing = [path for path in outputs if path.exists()]
    if existing and not force:
        names = ", ".join(str(path) for path in existing)
        raise FileExistsError(f"Output already exists; use --force to replace: {names}")


def map_values_to_atlas(
    atlas: np.ndarray,
    labels: np.ndarray,
    values: np.ndarray,
) -> np.ndarray:
    mapped = np.zeros(atlas.shape, dtype=np.float32)
    for label, value in zip(labels, values, strict=True):
        mapped[atlas == label] = value
    return mapped


def parcel_names(path: Path, expected_labels: np.ndarray) -> list[str]:
    rows = [line.split() for line in path.read_text(encoding="utf-8").splitlines()]
    labels = np.asarray([int(row[0]) for row in rows], dtype=int)
    if not np.array_equal(labels, expected_labels):
        raise ValueError("Schaefer label file must contain ordered labels 1 through 1000.")
    return [row[1] for row in rows]


def minmax(values: pd.Series, description: str) -> pd.Series:
    span = float(values.max() - values.min())
    if not np.isfinite(span) or span <= 0:
        raise ValueError(f"The {description} E:I ratio has no finite positive range.")
    return (values - values.min()) / span


def main() -> None:
    args = parse_args()
    config_path = args.config.expanduser().resolve()
    with config_path.open(encoding="utf-8") as stream:
        config = json.load(stream)

    try:
        import abagen
    except ImportError as exc:
        raise RuntimeError(
            "abagen is required. Activate .venv-neuromaps and install "
            "neuromaps_analysis/requirements.txt."
        ) from exc

    enable_abagen_pandas_compatibility()

    data_dir = resolve_from_config(config_path, config["data_dir"])
    atlas_path = resolve_from_config(config_path, config["atlas"])
    labels_path = resolve_from_config(config_path, config["atlas_labels"])
    output_dir = resolve_from_config(config_path, config["output_dir"])
    for label, path in (
        ("abagen data directory", data_dir),
        ("Schaefer atlas", atlas_path),
        ("Schaefer label file", labels_path),
    ):
        if not path.exists():
            raise FileNotFoundError(f"Missing {label}: {path}")

    output_dir.mkdir(parents=True, exist_ok=True)
    expression_path = output_dir / "expression_sch1000_bilateral_paper_settings.csv"
    table_path = output_dir / "EI_values_sch1000.csv"
    left_table_path = output_dir / "EI_values_sch500_left.csv"
    raw_map_path = output_dir / "EI_expression_raw_sch1000_2mm.nii.gz"
    normalized_map_path = output_dir / "EI_expression_minmax_sch1000_2mm.nii.gz"
    report_path = output_dir / "abagen_report_sch1000_bilateral.txt"
    metadata_path = output_dir / "EI_build_metadata.json"
    outputs = [
        expression_path,
        table_path,
        left_table_path,
        raw_map_path,
        normalized_map_path,
        report_path,
        metadata_path,
    ]
    ensure_available(outputs, args.force)

    settings = config["abagen"]
    get_expression_args = {
        key: settings[key]
        for key in (
            "donors",
            "ibf_threshold",
            "probe_selection",
            "donor_probes",
            "tolerance",
            "sample_norm",
            "gene_norm",
            "norm_matched",
            "lr_mirror",
            "missing",
            "region_agg",
            "agg_metric",
            "corrected_mni",
            "reannotated",
        )
    }
    requested_settings = {
        "probe_selection": "rnaseq",
        "sample_norm": None,
        "gene_norm": "srs",
    }
    for key, expected in requested_settings.items():
        if get_expression_args[key] != expected:
            raise ValueError(
                f"Configuration must set abagen.{key}={expected!r}; "
                f"found {get_expression_args[key]!r}."
            )
    if get_expression_args["lr_mirror"] is not None:
        raise ValueError("Bilateral estimation requires abagen.lr_mirror=null.")

    # abagen 0.1.3's RNA-seq selector makes an internal fetch call without
    # forwarding data_dir. Scope that call to the repository-local cache.
    original_fetch_rnaseq = abagen.datasets.fetch_rnaseq

    def fetch_rnaseq_scoped(*fetch_args: object, **fetch_kwargs: object) -> object:
        if fetch_kwargs.get("data_dir") is None:
            fetch_kwargs["data_dir"] = str(data_dir)
        return original_fetch_rnaseq(*fetch_args, **fetch_kwargs)

    abagen.datasets.fetch_rnaseq = fetch_rnaseq_scoped
    print(f"Using abagen data directory: {data_dir}")
    print(f"Estimating expression independently in both hemispheres: {atlas_path}")
    try:
        expression, report = abagen.get_expression_data(
            str(atlas_path),
            data_dir=str(data_dir),
            return_report=True,
            **get_expression_args,
        )
    finally:
        abagen.datasets.fetch_rnaseq = original_fetch_rnaseq
    if not isinstance(expression, pd.DataFrame):
        raise TypeError("Expected abagen to return one aggregated expression table.")

    expected_labels = np.arange(1, int(config["n_parcels"]) + 1)
    try:
        expression.index = pd.Index(expression.index).astype(int)
    except (TypeError, ValueError) as exc:
        raise ValueError("abagen returned non-numeric parcel identifiers.") from exc
    expression = expression.sort_index()
    if not expression.index.is_unique or not np.array_equal(
        expression.index.to_numpy(), expected_labels
    ):
        raise ValueError("abagen must return unique parcel labels 1 through 1000.")

    excitatory = list(config["excitatory_genes"])
    inhibitory = list(config["inhibitory_genes"])
    requested_genes = [*excitatory, *inhibitory]
    missing_genes = [gene for gene in requested_genes if gene not in expression.columns]
    if missing_genes:
        raise KeyError(f"Processed expression matrix is missing genes: {missing_genes}")
    selected = expression.loc[:, requested_genes].astype(float)
    if not np.isfinite(selected.to_numpy()).all():
        bad = selected.index[~np.isfinite(selected).all(axis=1)].tolist()
        raise ValueError(f"Non-finite selected-gene expression in parcels: {bad}")

    excitation_sum = selected[excitatory].sum(axis=1)
    inhibition_sum = selected[inhibitory].sum(axis=1)
    if (inhibition_sum <= 0).any():
        bad = inhibition_sum.index[inhibition_sum <= 0].tolist()
        raise ValueError(f"Non-positive inhibitory expression sums in parcels: {bad}")
    ei_raw = excitation_sum / inhibition_sum
    ei_minmax_bilateral = minmax(ei_raw, "bilateral")

    left_count = int(config["n_left_parcels"])
    left_labels = np.arange(1, left_count + 1)
    left_ei_raw = ei_raw.loc[left_labels]
    ei_minmax_left = minmax(left_ei_raw, "left-hemisphere")
    names = parcel_names(labels_path, expected_labels)
    hemispheres = np.where(expected_labels <= left_count, "L", "R")

    values = pd.DataFrame(
        {
            "label": expected_labels,
            "parcel_name": names,
            "hemisphere": hemispheres,
            "excitation_sum": excitation_sum.to_numpy(),
            "inhibition_sum": inhibition_sum.to_numpy(),
            "ei_raw": ei_raw.to_numpy(),
            "ei_minmax": ei_minmax_bilateral.to_numpy(),
        }
    )
    left_values = values.iloc[:left_count].copy()
    left_values.insert(0, "left_order", np.arange(1, left_count + 1))
    left_values["ei_minmax"] = ei_minmax_left.to_numpy()

    atlas_img = nib.load(atlas_path)
    atlas = np.asarray(atlas_img.dataobj)
    atlas_labels = np.unique(atlas[atlas > 0]).astype(int)
    if not np.array_equal(atlas_labels, expected_labels):
        raise ValueError("Atlas must contain integer labels 1 through 1000 exactly.")
    raw_map = map_values_to_atlas(atlas, expected_labels, ei_raw.to_numpy())
    normalized_map = map_values_to_atlas(
        atlas, expected_labels, ei_minmax_bilateral.to_numpy()
    )
    header = atlas_img.header.copy()
    header.set_data_dtype(np.float32)
    nib.save(nib.Nifti1Image(raw_map, atlas_img.affine, header), raw_map_path)
    nib.save(
        nib.Nifti1Image(normalized_map, atlas_img.affine, header), normalized_map_path
    )

    expression.to_csv(expression_path, index_label="label", float_format="%.17g")
    values.to_csv(table_path, index=False, float_format="%.17g")
    left_values.to_csv(left_table_path, index=False, float_format="%.17g")
    report_path.write_text(str(report).strip() + "\n", encoding="utf-8")

    metadata = {
        "map_name": "primary_transcriptomic_ei_sch1000",
        "formula": "sum(excitatory genes) / sum(inhibitory genes)",
        "normalization": config["normalization"],
        "hemisphere_strategy": config["hemisphere_strategy"],
        "statistical_subset": config["statistical_subset"],
        "n_parcels": int(config["n_parcels"]),
        "n_left_parcels": left_count,
        "excitatory_genes": excitatory,
        "inhibitory_genes": inhibitory,
        "abagen_version": abagen.__version__,
        "abagen_settings": get_expression_args,
        "inputs": {
            "data_dir": str(data_dir),
            "atlas": str(atlas_path),
            "atlas_sha256": sha256(atlas_path),
            "atlas_labels": str(labels_path),
            "atlas_labels_sha256": sha256(labels_path),
        },
        "outputs": {
            "expression_bilateral": expression_path.name,
            "parcel_values_bilateral": table_path.name,
            "parcel_values_left_for_statistics": left_table_path.name,
            "raw_nifti_bilateral": raw_map_path.name,
            "minmax_nifti_bilateral": normalized_map_path.name,
            "abagen_report": report_path.name,
        },
        "summary": {
            "bilateral_ei_raw_min": float(ei_raw.min()),
            "bilateral_ei_raw_max": float(ei_raw.max()),
            "bilateral_ei_minmax_min": float(ei_minmax_bilateral.min()),
            "bilateral_ei_minmax_max": float(ei_minmax_bilateral.max()),
            "left_ei_raw_min": float(left_ei_raw.min()),
            "left_ei_raw_max": float(left_ei_raw.max()),
            "left_ei_minmax_min": float(ei_minmax_left.min()),
            "left_ei_minmax_max": float(ei_minmax_left.max()),
        },
        "method_note": config["method_note"],
    }
    with metadata_path.open("w", encoding="utf-8") as stream:
        json.dump(metadata, stream, indent=2)
        stream.write("\n")

    print(f"Saved bilateral expression matrix: {expression_path}")
    print(f"Saved bilateral 1,000-parcel E:I values: {table_path}")
    print(f"Saved left-only 500-parcel statistical values: {left_table_path}")
    print(f"Saved bilateral raw E:I map: {raw_map_path}")
    print(f"Saved bilateral normalized E:I map: {normalized_map_path}")
    print(f"Saved abagen report: {report_path}")
    print(f"Saved build metadata: {metadata_path}")


if __name__ == "__main__":
    main()
