#!/usr/bin/env python3
"""Build a paper-matched Desikan--Killiany E:I sanity-check map.

The AHBA data are parcellated with abagen in the Desikan--Killiany atlas.
Only the 34 left cortical parcels are used to calculate E:I; those values are
then reflected to the homologous right-hemisphere parcels for DK68 rendering.
"""

from __future__ import annotations

import argparse
import json
import tempfile
from pathlib import Path

import nibabel as nib
import numpy as np
import pandas as pd


ABAGEN_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_CONFIG = ABAGEN_ROOT / "config" / "ei_dk68.json"


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

    # abagen 0.1.3 passes ``inplace=False`` to DataFrame.set_axis(); pandas
    # 2.x removed that keyword. Keep the compatibility shim local to this
    # command-line process instead of altering the installed dependency.
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
        help="Overwrite an existing DK68 E:I build.",
    )
    return parser.parse_args()


def resolve_from_config(config_path: Path, value: str) -> Path:
    path = Path(value).expanduser()
    if not path.is_absolute():
        path = config_path.parent / path
    return path.resolve()


def require_new_outputs(paths: list[Path], force: bool) -> None:
    existing = [path for path in paths if path.exists()]
    if existing and not force:
        joined = ", ".join(str(path) for path in existing)
        raise FileExistsError(
            f"Output already exists; use --force to replace: {joined}"
        )


def cortical_atlas_info(atlas_info_path: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    info = pd.read_csv(atlas_info_path)
    required = {"id", "label", "hemisphere", "structure"}
    missing = required.difference(info.columns)
    if missing:
        raise KeyError(f"DK atlas information is missing columns: {sorted(missing)}")

    cortex = info.loc[info["structure"].eq("cortex")].copy()
    left = cortex.loc[cortex["hemisphere"].eq("L")].sort_values("id")
    right = cortex.loc[cortex["hemisphere"].eq("R")].sort_values("id")
    if len(left) != 34 or len(right) != 34:
        raise ValueError(
            "Expected 34 cortical Desikan--Killiany parcels per hemisphere; "
            f"found {len(left)} left and {len(right)} right."
        )
    if left["label"].tolist() != right["label"].tolist():
        raise ValueError("Left and right DK cortical parcel orders do not match.")
    return left.reset_index(drop=True), right.reset_index(drop=True)


def select_left_expression(
    expression: pd.DataFrame,
    left_info: pd.DataFrame,
) -> pd.DataFrame:
    """Return expression rows in renderer order, with strict ID checking."""
    try:
        numeric_index = pd.Index(expression.index).astype(int)
    except (TypeError, ValueError) as exc:
        raise ValueError("abagen returned non-numeric parcel identifiers.") from exc

    expression = expression.copy()
    expression.index = numeric_index
    if not expression.index.is_unique:
        raise ValueError("abagen returned duplicate parcel identifiers.")

    ids = left_info["id"].astype(int).to_numpy()
    absent = sorted(set(ids).difference(expression.index))
    if absent:
        raise ValueError(f"No expression data were returned for left DK IDs: {absent}")
    return expression.loc[ids]


def make_left_cortex_atlases(
    atlas_images: dict[str, str],
    left_ids: np.ndarray,
    output_dir: Path,
) -> dict[str, str]:
    """Mask donor-native DK atlases before matching and normalization."""
    masked_images: dict[str, str] = {}
    for donor, atlas_path in atlas_images.items():
        image = nib.load(atlas_path)
        atlas = np.asanyarray(image.dataobj)
        masked = np.where(np.isin(atlas, left_ids), atlas, 0)
        header = image.header.copy()
        header.set_data_dtype(atlas.dtype)
        masked_path = output_dir / f"atlas-dk34-left-donor{donor}.nii.gz"
        nib.save(nib.Nifti1Image(masked, image.affine, header), masked_path)
        masked_images[str(donor)] = str(masked_path)
    return masked_images


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
    output_dir = resolve_from_config(config_path, config["output_dir"])
    output_dir.mkdir(parents=True, exist_ok=True)

    expression_path = output_dir / "expression_dk34_left_cortex.csv"
    left_path = output_dir / "EI_values_dk34_left.csv"
    render_path = output_dir / "EI_values_dk68_reflected.csv"
    raw_vector_path = output_dir / "EI_raw_dk68.txt"
    minmax_vector_path = output_dir / "EI_minmax_dk68.txt"
    report_path = output_dir / "abagen_report_dk68.txt"
    metadata_path = output_dir / "EI_build_metadata.json"
    outputs = [
        expression_path,
        left_path,
        render_path,
        raw_vector_path,
        minmax_vector_path,
        report_path,
        metadata_path,
    ]
    require_new_outputs(outputs, args.force)

    settings = config["abagen"]
    if settings["atlas"] != "desikan_killiany":
        raise ValueError("This builder supports only the Desikan--Killiany atlas.")

    atlas = abagen.fetch_desikan_killiany(native=bool(settings["native"]))
    atlas_info_path = Path(atlas["info"]).resolve()
    left_info, right_info = cortical_atlas_info(atlas_info_path)

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
    # abagen 0.1.3's RNA-seq selector makes a second fetch call without
    # forwarding data_dir. Scope that internal call to the configured project
    # cache so results do not depend on a user's home-directory cache.
    original_fetch_rnaseq = abagen.datasets.fetch_rnaseq

    def fetch_rnaseq_scoped(*fetch_args: object, **fetch_kwargs: object) -> object:
        if fetch_kwargs.get("data_dir") is None:
            fetch_kwargs["data_dir"] = str(data_dir)
        return original_fetch_rnaseq(*fetch_args, **fetch_kwargs)

    abagen.datasets.fetch_rnaseq = fetch_rnaseq_scoped
    print(f"Using abagen data directory: {data_dir}")
    try:
        restrict_left = bool(
            config.get("restrict_to_left_cortex_before_normalization", False)
        )
        if restrict_left:
            left_ids = left_info["id"].astype(int).to_numpy()
            with tempfile.TemporaryDirectory(prefix="ei_dk34_left_") as tmp:
                masked_atlases = make_left_cortex_atlases(
                    atlas["image"], left_ids, Path(tmp)
                )
                expression, report = abagen.get_expression_data(
                    masked_atlases,
                    left_info,
                    data_dir=str(data_dir),
                    return_report=True,
                    **get_expression_args,
                )
        else:
            expression, report = abagen.get_expression_data(
                atlas["image"],
                atlas["info"],
                data_dir=str(data_dir),
                return_report=True,
                **get_expression_args,
            )
    finally:
        abagen.datasets.fetch_rnaseq = original_fetch_rnaseq
    if not isinstance(expression, pd.DataFrame):
        raise TypeError("Expected abagen to return one aggregated expression table.")

    excitatory = list(config["excitatory_genes"])
    inhibitory = list(config["inhibitory_genes"])
    requested = [*excitatory, *inhibitory]
    missing_genes = [gene for gene in requested if gene not in expression.columns]
    if missing_genes:
        raise KeyError(f"Processed expression matrix is missing genes: {missing_genes}")

    left_expression = select_left_expression(expression, left_info)
    selected = left_expression.loc[:, requested].astype(float)
    if not np.isfinite(selected.to_numpy()).all():
        bad_ids = selected.index[~np.isfinite(selected).all(axis=1)].tolist()
        raise ValueError(f"Non-finite selected-gene expression in DK IDs: {bad_ids}")

    excitation_sum = selected[excitatory].sum(axis=1)
    inhibition_sum = selected[inhibitory].sum(axis=1)
    if (inhibition_sum <= 0).any():
        bad_ids = inhibition_sum.index[inhibition_sum <= 0].tolist()
        raise ValueError(f"Non-positive inhibitory sums in DK IDs: {bad_ids}")

    ei_raw = excitation_sum / inhibition_sum
    span = float(ei_raw.max() - ei_raw.min())
    if not np.isfinite(span) or span <= 0:
        raise ValueError("The left-hemisphere E:I ratio has no finite positive range.")
    ei_minmax = (ei_raw - ei_raw.min()) / span

    left_values = left_info.loc[:, ["id", "label", "hemisphere"]].copy()
    left_values.insert(0, "left_order", np.arange(1, 35))
    left_values["excitation_sum"] = excitation_sum.to_numpy()
    left_values["inhibition_sum"] = inhibition_sum.to_numpy()
    left_values["ei_raw"] = ei_raw.to_numpy()
    left_values["ei_minmax"] = ei_minmax.to_numpy()

    left_render = left_values.copy()
    left_render["source_left_id"] = left_values["id"].to_numpy()
    left_render["source_left_label"] = left_values["label"].to_numpy()

    right_render = left_values.copy()
    right_render["id"] = right_info["id"].to_numpy()
    right_render["label"] = right_info["label"].to_numpy()
    right_render["hemisphere"] = "R"
    right_render["source_left_id"] = left_values["id"].to_numpy()
    right_render["source_left_label"] = left_values["label"].to_numpy()

    render_values = pd.concat([left_render, right_render], ignore_index=True)
    render_values.insert(0, "render_index", np.arange(1, 69))
    if render_values["label"].iloc[:34].tolist() != render_values[
        "label"
    ].iloc[34:].tolist():
        raise AssertionError("Reflected DK68 parcel order is inconsistent.")

    expression.to_csv(expression_path, float_format="%.17g")
    left_values.to_csv(left_path, index=False, float_format="%.17g")
    render_values.to_csv(render_path, index=False, float_format="%.17g")
    np.savetxt(raw_vector_path, render_values["ei_raw"], fmt="%.17g")
    np.savetxt(minmax_vector_path, render_values["ei_minmax"], fmt="%.17g")
    report_path.write_text(str(report).strip() + "\n", encoding="utf-8")

    metadata = {
        "map_name": "paper_matched_transcriptomic_ei_dk68_sanity_check",
        "formula": "sum(excitatory genes) / sum(inhibitory genes)",
        "normalization": config["normalization"],
        "hemisphere_strategy": config["hemisphere_strategy"],
        "restrict_to_left_cortex_before_normalization": restrict_left,
        "n_left_parcels": int(config["n_left_parcels"]),
        "n_render_parcels": int(config["n_render_parcels"]),
        "excitatory_genes": excitatory,
        "inhibitory_genes": inhibitory,
        "abagen_version": abagen.__version__,
        "abagen_settings": get_expression_args,
        "inputs": {
            "data_dir": str(data_dir),
            "atlas_info": str(atlas_info_path),
            "atlas_images": atlas["image"],
        },
        "outputs": {
            "expression_left_cortex": expression_path.name,
            "left_parcel_values": left_path.name,
            "reflected_render_values": render_path.name,
            "raw_render_vector": raw_vector_path.name,
            "minmax_render_vector": minmax_vector_path.name,
            "abagen_report": report_path.name,
        },
        "summary": {
            "ei_raw_min": float(ei_raw.min()),
            "ei_raw_max": float(ei_raw.max()),
            "ei_minmax_min": float(ei_minmax.min()),
            "ei_minmax_max": float(ei_minmax.max()),
        },
        "method_note": config["method_note"],
    }
    with metadata_path.open("w", encoding="utf-8") as stream:
        json.dump(metadata, stream, indent=2)
        stream.write("\n")

    print(f"Saved left-cortex expression: {expression_path}")
    print(f"Saved 34 left-hemisphere values: {left_path}")
    print(f"Saved reflected DK68 values: {render_path}")
    print(f"Saved raw vector: {raw_vector_path}")
    print(f"Saved min-max vector: {minmax_vector_path}")
    print(f"Saved abagen report: {report_path}")
    print(f"Saved metadata: {metadata_path}")


if __name__ == "__main__":
    main()
