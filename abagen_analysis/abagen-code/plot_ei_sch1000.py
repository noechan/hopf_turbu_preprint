#!/usr/bin/env python3
"""Visualize the raw and min-max-normalized Schaefer-1000 E:I maps."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import nibabel as nib
import numpy as np
from nilearn import datasets, plotting


SCRIPT_DIR = Path(__file__).resolve().parent
DEFAULT_MAP_DIR = SCRIPT_DIR / "ei_maps" / "schaefer1000_2mm"
DEFAULT_OUTPUT_DIR = DEFAULT_MAP_DIR / "figures"
RAW_FILENAME = "EI_expression_raw_sch1000_2mm.nii.gz"
MINMAX_FILENAME = "EI_expression_minmax_sch1000_2mm.nii.gz"
CUT_COORDS = (-40, -20, 0, 20, 40, 60)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--map-dir",
        type=Path,
        default=DEFAULT_MAP_DIR,
        help="Directory containing the raw and min-max E:I NIfTI maps.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=DEFAULT_OUTPUT_DIR,
        help="Directory in which PNG and PDF figures will be written.",
    )
    parser.add_argument("--dpi", type=int, default=300)
    parser.add_argument(
        "--force",
        action="store_true",
        help="Overwrite existing visualization outputs.",
    )
    return parser.parse_args()


def load_and_validate_maps(map_dir: Path) -> tuple[nib.Nifti1Image, nib.Nifti1Image, np.ndarray]:
    raw_path = map_dir / RAW_FILENAME
    minmax_path = map_dir / MINMAX_FILENAME
    for path in (raw_path, minmax_path):
        if not path.is_file():
            raise FileNotFoundError(f"Missing E:I map: {path}")

    raw_img = nib.load(raw_path)
    minmax_img = nib.load(minmax_path)
    if raw_img.shape != minmax_img.shape or not np.allclose(
        raw_img.affine, minmax_img.affine
    ):
        raise ValueError("Raw and min-max E:I maps must have matching geometry.")

    raw_data = np.asarray(raw_img.dataobj, dtype=float)
    minmax_data = np.asarray(minmax_img.dataobj, dtype=float)
    if not np.isfinite(raw_data).all() or not np.isfinite(minmax_data).all():
        raise ValueError("E:I maps contain non-finite values.")

    # The raw ratio is strictly positive inside the atlas and zero outside it.
    atlas_mask = raw_data != 0
    if not atlas_mask.any():
        raise ValueError("The raw E:I map contains no nonzero atlas voxels.")

    raw_values = raw_data[atlas_mask]
    minmax_values = minmax_data[atlas_mask]
    expected = (raw_values - raw_values.min()) / (
        raw_values.max() - raw_values.min()
    )
    if not np.allclose(minmax_values, expected, atol=1e-5, rtol=1e-5):
        raise ValueError("The normalized map is not a min-max transform of the raw map.")

    return raw_img, minmax_img, atlas_mask


def plot_one_map(
    image: nib.Nifti1Image,
    background: nib.Nifti1Image,
    title: str,
    output_path: Path,
    vmin: float,
    vmax: float,
    threshold: float,
    dpi: int,
) -> None:
    figure = plt.figure(figsize=(13.5, 3.1), facecolor="white")
    figure.suptitle(title, fontsize=16, fontweight="semibold", y=0.98)
    axes = figure.add_axes((0.02, 0.06, 0.96, 0.82))
    display = plotting.plot_stat_map(
        image,
        display_mode="z",
        cut_coords=CUT_COORDS,
        cmap="viridis",
        colorbar=True,
        symmetric_cbar=False,
        threshold=threshold,
        vmin=vmin,
        vmax=vmax,
        bg_img=background,
        black_bg=False,
        annotate=True,
        draw_cross=False,
        radiological=False,
        title=None,
        figure=figure,
        axes=axes,
    )
    display.savefig(output_path, dpi=dpi)
    display.close()


def plot_combined(
    raw_img: nib.Nifti1Image,
    minmax_img: nib.Nifti1Image,
    background: nib.Nifti1Image,
    raw_limits: tuple[float, float],
    output_path: Path,
    dpi: int,
) -> None:
    figure = plt.figure(figsize=(13.5, 6.4), facecolor="white")
    specifications = (
        (
            raw_img,
            "Raw transcriptomic E:I ratio",
            (0.02, 0.53, 0.96, 0.39),
            raw_limits,
            raw_limits[0] - 1e-6,
        ),
        (
            minmax_img,
            "Min–max-normalized transcriptomic E:I ratio",
            (0.02, 0.06, 0.96, 0.39),
            (0.0, 1.0),
            0.0,
        ),
    )

    displays = []
    for index, (image, title, axes_position, limits, threshold) in enumerate(specifications):
        figure.text(
            0.5,
            0.965 if index == 0 else 0.495,
            title,
            ha="center",
            va="top",
            fontsize=16,
            fontweight="semibold",
        )
        axes = figure.add_axes(axes_position)
        display = plotting.plot_stat_map(
            image,
            display_mode="z",
            cut_coords=CUT_COORDS,
            cmap="viridis",
            colorbar=True,
            symmetric_cbar=False,
            threshold=threshold,
            vmin=limits[0],
            vmax=limits[1],
            bg_img=background,
            black_bg=False,
            annotate=True,
            draw_cross=False,
            radiological=False,
            title=None,
            figure=figure,
            axes=axes,
        )
        displays.append(display)

    figure.savefig(output_path, dpi=dpi, bbox_inches="tight", facecolor="white")
    for display in displays:
        display.close()


def main() -> None:
    args = parse_args()
    map_dir = args.map_dir.expanduser().resolve()
    output_dir = args.output_dir.expanduser().resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    outputs = {
        "raw_png": output_dir / "EI_expression_raw_sch1000_slices.png",
        "minmax_png": output_dir / "EI_expression_minmax_sch1000_slices.png",
        "combined_png": output_dir / "EI_expression_raw_and_minmax_sch1000_slices.png",
        "combined_pdf": output_dir / "EI_expression_raw_and_minmax_sch1000_slices.pdf",
        "metadata": output_dir / "EI_visualization_metadata.json",
    }
    existing = [path for path in outputs.values() if path.exists()]
    if existing and not args.force:
        names = ", ".join(str(path) for path in existing)
        raise FileExistsError(f"Output already exists; use --force to replace: {names}")

    raw_img, minmax_img, atlas_mask = load_and_validate_maps(map_dir)
    background = datasets.load_mni152_template(resolution=2)
    raw_data = np.asarray(raw_img.dataobj, dtype=float)
    minmax_data = np.asarray(minmax_img.dataobj, dtype=float)
    raw_limits = (float(raw_data[atlas_mask].min()), float(raw_data[atlas_mask].max()))
    plot_one_map(
        raw_img,
        background,
        "Raw transcriptomic E:I ratio",
        outputs["raw_png"],
        *raw_limits,
        raw_limits[0] - 1e-6,
        args.dpi,
    )
    plot_one_map(
        minmax_img,
        background,
        "Min–max-normalized transcriptomic E:I ratio",
        outputs["minmax_png"],
        0.0,
        1.0,
        0.0,
        args.dpi,
    )
    plot_combined(
        raw_img,
        minmax_img,
        background,
        raw_limits,
        outputs["combined_png"],
        args.dpi,
    )
    plot_combined(
        raw_img,
        minmax_img,
        background,
        raw_limits,
        outputs["combined_pdf"],
        args.dpi,
    )

    metadata = {
        "raw_map": str((map_dir / RAW_FILENAME).resolve()),
        "minmax_map": str((map_dir / MINMAX_FILENAME).resolve()),
        "atlas_voxels": int(atlas_mask.sum()),
        "raw_range": list(raw_limits),
        "minmax_range": [
            float(minmax_data[atlas_mask].min()),
            float(minmax_data[atlas_mask].max()),
        ],
        "display_mode": "axial",
        "cut_coordinates_mni_z_mm": list(CUT_COORDS),
        "colormap": "viridis",
        "note": (
            "Min-max normalization is a monotonic linear rescaling; therefore "
            "the two maps have identical spatial ordering and differ only in units."
        ),
        "outputs": {name: path.name for name, path in outputs.items() if name != "metadata"},
    }
    with outputs["metadata"].open("w", encoding="utf-8") as stream:
        json.dump(metadata, stream, indent=2)
        stream.write("\n")

    for path in outputs.values():
        print(f"Saved: {path}")


if __name__ == "__main__":
    main()
