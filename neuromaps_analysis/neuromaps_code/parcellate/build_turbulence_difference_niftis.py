#!/usr/bin/env python3
"""Map the six aggregate Schaefer-1000 turbulence differences to MNI152."""

from __future__ import annotations

import argparse
from pathlib import Path

import nibabel as nib
import numpy as np
from scipy.io import loadmat


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_DATA_DIR = (
    ROOT / "neuromaps-data" / "annotations" / "turbu" / "MNI152" / "N145_sch1000"
)
DEFAULT_ATLAS = (
    ROOT / "nilearn_data" / "schaefer_2018"
    / "Schaefer2018_1000Parcels_7Networks_order_FSLMNI152_2mm.nii.gz"
)

INPUTS = {
    ("hc_ad", 1): (
        "brain_nodes_diff_hcneg_adpos_lam1_sch1000_harmonized_N145.mat",
        "diff_hc_ad_lam1",
    ),
    ("hc_ad", 3): (
        "brain_nodes_diff_hcneg_adpos_lam3_sch1000_harmonized_N145.mat",
        "diff_hc_ad_lam3",
    ),
    ("hc_ad", 6): (
        "brain_nodes_diff_hcneg_adpos_lam6_sch1000_harmonized_N145.mat",
        "diff_hc_ad_lam6",
    ),
    ("mci_ad", 1): (
        "brain_nodes_diff_mcipos_adpos_lam1_sch1000_harmonized_N145.mat",
        "diff_mci_ad_lam1",
    ),
    ("mci_ad", 3): (
        "brain_nodes_diff_mcipos_adpos_lam3_sch1000_harmonized_N145.mat",
        "diff_mci_ad_lam3",
    ),
    ("mci_ad", 6): (
        "brain_nodes_diff_mcipos_adpos_lam6_sch1000_harmonized_N145.mat",
        "diff_mci_ad_lam6",
    ),
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-dir", type=Path, default=DEFAULT_DATA_DIR)
    parser.add_argument("--atlas", type=Path, default=DEFAULT_ATLAS)
    parser.add_argument(
        "--contrast", choices=("all", "hc_ad", "mci_ad"), default="all"
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    data_dir = args.data_dir.expanduser().resolve()
    atlas_path = args.atlas.expanduser().resolve()
    if not atlas_path.is_file():
        raise FileNotFoundError(f"Missing Schaefer atlas: {atlas_path}")

    atlas_img = nib.load(atlas_path)
    atlas = np.asarray(atlas_img.dataobj)
    labels = np.unique(atlas[atlas > 0]).astype(int)
    if labels.size != 1000 or not np.array_equal(labels, np.arange(1, 1001)):
        raise ValueError("Atlas must contain integer labels 1 through 1000.")

    selected = {
        key: value
        for key, value in INPUTS.items()
        if args.contrast == "all" or key[0] == args.contrast
    }
    for (contrast, suffix), (mat_name, variable) in selected.items():
        mat_path = data_dir / mat_name
        if not mat_path.is_file():
            raise FileNotFoundError(f"Missing aggregate input: {mat_path}")
        mat = loadmat(mat_path)
        if variable not in mat:
            raise KeyError(f"{variable!r} is absent from {mat_path.name}")
        values = np.asarray(mat[variable], dtype=float).squeeze()
        if values.shape != (1000,) or not np.all(np.isfinite(values)):
            raise ValueError(f"{variable} must be 1,000 finite parcel values.")

        mapped = np.zeros(atlas.shape, dtype=np.float32)
        for label, value in enumerate(values, start=1):
            mapped[atlas == label] = value

        output = data_dir / (
            f"turbu_diff_{contrast}_lam{suffix}_sch1000_N145_2mm.nii.gz"
        )
        nib.save(nib.Nifti1Image(mapped, atlas_img.affine, atlas_img.header), output)
        print(f"Saved {output}")


if __name__ == "__main__":
    main()
