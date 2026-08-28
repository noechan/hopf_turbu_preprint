#!/usr/bin/env python3
"""Run the planned N145 turbulence--Neurosynth spatial correlations."""

from __future__ import annotations

import argparse
import os
from pathlib import Path
import tempfile

ROOT = Path(__file__).resolve().parents[2]
DEFAULT_CACHE = Path(os.environ.get("XDG_CACHE_HOME", Path.home() / ".cache")) / "neuromaps"
os.environ.setdefault("NEUROMAPS_DATA", str(DEFAULT_CACHE))
os.environ.setdefault(
    "MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "neuromaps_matplotlib")
)

import numpy as np
import pandas as pd
from neuromaps import nulls, transforms
from neuromaps.stats import compare_images
from statsmodels.stats.multitest import multipletests


TURBULENCE_DIR = (
    ROOT / "neuromaps-data" / "annotations" / "turbu" / "MNI152" / "N145_sch1000"
)
NEUROSYNTH_DIR = (
    ROOT / "neuromaps-data" / "annotations" / "neurosynth_memo" / "MNI152"
)
TERMS = {
    "memory": "parcellated_memory_association-test_z_FDR_0.01_sch1000_2mm.nii.gz",
    "episodic_memory": (
        "parcellated_episodic memory_association-test_z_FDR_0.01_sch1000_2mm.nii.gz"
    ),
    "autobiographical_memory": (
        "parcellated_autobiographical memory_association-test_z_FDR_0.01_sch1000_2mm.nii.gz"
    ),
}
LAMBDAS = {1: 0.01, 3: 0.03, 6: 0.06}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--contrast", choices=("all", "hc_ad", "mci_ad"), default="all"
    )
    parser.add_argument("--n-perm", type=int, default=1000)
    parser.add_argument("--seed", type=int, default=1234)
    parser.add_argument(
        "--output-dir", type=Path, default=ROOT / "results" / "neurosynth"
    )
    return parser.parse_args()


def surface_to_array(surface: tuple) -> np.ndarray:
    """Concatenate in-memory GIFTI hemispheres (neuromaps 0.0.5)."""
    return np.hstack(
        [np.asarray(hemisphere.agg_data(), dtype=float).squeeze() for hemisphere in surface]
    )


def project(path: Path) -> np.ndarray:
    if not path.is_file():
        raise FileNotFoundError(f"Missing spatial map: {path}")
    return surface_to_array(transforms.mni152_to_fsaverage(str(path), "41k"))


def run_contrast(contrast: str, n_perm: int, seed: int) -> pd.DataFrame:
    cognition = {term: project(NEUROSYNTH_DIR / name) for term, name in TERMS.items()}
    cognitive_nulls = {
        term: nulls.alexander_bloch(
            values,
            atlas="fsaverage",
            density="41k",
            n_perm=n_perm,
            seed=seed,
        )
        for term, values in cognition.items()
    }

    rows = []
    for suffix, physical_lambda in LAMBDAS.items():
        turbulence = project(
            TURBULENCE_DIR
            / f"turbu_diff_{contrast}_lam{suffix}_sch1000_N145_2mm.nii.gz"
        )
        for term, cognitive_values in cognition.items():
            r_value, p_spin = compare_images(
                turbulence,
                cognitive_values,
                nulls=cognitive_nulls[term],
                metric="pearsonr",
            )
            rows.append(
                {
                    "contrast": contrast,
                    "lambda_suffix": suffix,
                    "lambda": physical_lambda,
                    "term": term,
                    "pearson_r": float(r_value),
                    "p_spin": float(p_spin),
                    "n_perm": n_perm,
                    "seed": seed,
                }
            )

    frame = pd.DataFrame(rows)
    frame["p_FDR_within_contrast"] = multipletests(
        frame["p_spin"].to_numpy(), method="fdr_bh"
    )[1]
    frame["significant_FDR_0_05"] = frame["p_FDR_within_contrast"] < 0.05
    return frame


def main() -> None:
    args = parse_args()
    if args.n_perm < 1:
        raise ValueError("--n-perm must be a positive integer.")
    contrasts = ("hc_ad", "mci_ad") if args.contrast == "all" else (args.contrast,)
    result = pd.concat(
        [run_contrast(c, args.n_perm, args.seed) for c in contrasts],
        ignore_index=True,
    )
    output_dir = args.output_dir.expanduser().resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    output = output_dir / "N145_turbulence_neurosynth_spatial_correlations.csv"
    result.to_csv(output, index=False)
    print(result.to_string(index=False))
    print(f"Saved {output}")


if __name__ == "__main__":
    main()
