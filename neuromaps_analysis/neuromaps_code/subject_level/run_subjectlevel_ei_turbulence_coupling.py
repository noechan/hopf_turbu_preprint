#!/usr/bin/env python3
"""Run the confirmatory participant-level E:I--turbulence analysis.

The scientific workflow is intentionally visible in ``main``. Implementation
details live in the neighboring ``ei_subjectlevel`` package:

1. Load the canonical ComBat-harmonized N145 turbulence maps.
2. Project turbulence and E:I maps to the left fsaverage surface.
3. Calculate one spatial E:I--turbulence coupling value per participant.
4. Test the ordered disease-stage trend and planned group contrasts.
5. Save participant values, inference tables, provenance, and the figure.
"""

from __future__ import annotations

import argparse
import os
from pathlib import Path
import tempfile


NEUROMAPS_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_CONFIG = NEUROMAPS_ROOT / "config" / "ei_subjectlevel_coupling_lam001.json"
DEFAULT_CACHE = NEUROMAPS_ROOT / "neuromaps-data" / "cache"
os.environ.setdefault("NEUROMAPS_DATA", str(DEFAULT_CACHE))
os.environ.setdefault(
    "MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "neuromaps_matplotlib")
)

import matplotlib

matplotlib.use("Agg")
import numpy as np
import pandas as pd

from ei_subjectlevel.inputs import (
    load_and_validate_inputs,
    read_config,
    resolve_from_config,
)
from ei_subjectlevel.outputs import (
    build_metadata,
    print_summary,
    save_metadata,
    save_tables,
)
from ei_subjectlevel.plotting import plot_participant_results
from ei_subjectlevel.spatial import compute_coupling, prepare_left_hemisphere_surfaces
from ei_subjectlevel.statistics import group_descriptives, run_all_models


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    parser.add_argument(
        "--skip-spins",
        action="store_true",
        help="Skip spatial normalization; intended only for quick diagnostics.",
    )
    parser.add_argument(
        "--between-subject-permutations",
        type=int,
        help="Override the configured Freedman-Lane permutation count.",
    )
    parser.add_argument(
        "--validate-inputs-only",
        action="store_true",
        help=(
            "Validate the canonical harmonized table, cohort, retained group "
            "map, atlas, and metadata without projecting surfaces or writing outputs."
        ),
    )
    return parser.parse_args()


def add_coupling_to_participants(
    participants: pd.DataFrame,
    participant_surface: np.ndarray,
    ei_surface: np.ndarray,
    model_config: dict,
    surface_config: dict,
    skip_spins: bool,
) -> tuple[pd.DataFrame, int]:
    """Calculate and attach the four participant-level coupling measures."""

    empirical_r, fisher_z, spin_z, spin_p, n_vertices = compute_coupling(
        participant_surface,
        ei_surface,
        spatial_rotations=int(model_config["spatial_rotations"]),
        seed=int(model_config["seed"]),
        atlas_name=surface_config["atlas"],
        density=surface_config["density"],
        hemisphere=surface_config["hemisphere"],
        skip_spins=skip_spins,
    )
    result = participants.copy()
    result["pearson_r"] = empirical_r
    result["fisher_z"] = fisher_z
    result["spin_z"] = spin_z
    result["spin_p_two_sided"] = spin_p
    return result, n_vertices


def main() -> None:
    args = parse_args()
    config_path = args.config.expanduser().resolve()
    config = read_config(config_path)

    # STEP 1 — Input data and provenance ------------------------------------
    print("\n[1/5] Loading and validating the canonical N145 inputs")
    inputs = load_and_validate_inputs(config_path, config)
    if args.validate_inputs_only:
        print(
            "Input validation passed: canonical harmonized N145 table, cohort, "
            "retained HC-AD map, demographics, E:I inputs, and atlas are consistent."
        )
        return

    surface_config = config["surface"]
    if (
        surface_config["hemisphere"] != "left"
        or int(surface_config["source_parcels"]) != 500
    ):
        raise ValueError("Participant-level correlations must use the 500 left parcels.")

    # STEP 2 — Common cortical representation -------------------------------
    print("\n[2/5] Projecting turbulence and E:I maps to the left cortical surface")
    participant_surface, ei_surface = prepare_left_hemisphere_surfaces(
        inputs.harmonized_nodes[inputs.parcel_columns].to_numpy(dtype=float),
        inputs.atlas_img,
        inputs.atlas_volume,
        inputs.paths.ei_nifti,
        density=surface_config["density"],
    )

    # STEP 3 — One coupling value per participant ---------------------------
    print("\n[3/5] Calculating participant E:I--turbulence spatial coupling")
    model_config = config["statistical_model"]
    participants, n_vertices = add_coupling_to_participants(
        inputs.participants,
        participant_surface,
        ei_surface,
        model_config,
        surface_config,
        skip_spins=args.skip_spins,
    )

    # STEP 4 — Between-participant confirmatory inference -------------------
    print("\n[4/5] Testing the ordered stage trend and planned group contrasts")
    n_between = (
        args.between_subject_permutations
        if args.between_subject_permutations is not None
        else int(model_config["between_subject_permutations"])
    )
    if n_between < 1000:
        raise ValueError("Use at least 1,000 between-subject permutations.")
    outcomes = [model_config["primary_outcome"]]
    if not args.skip_spins:
        outcomes.append(model_config["sensitivity_outcome"])
    comparisons, trends = run_all_models(
        participant_table=participants,
        outcomes=outcomes,
        group_order=list(config["group_order"]),
        n_perm=n_between,
        seed=int(model_config["seed"]),
    )
    descriptives = group_descriptives(
        participants,
        list(config["group_order"]),
        config["group_display_names"],
        seed=int(model_config["seed"]),
    )

    # STEP 5 — Tables, figure, and reproducibility metadata -----------------
    print("\n[5/5] Saving tables, figure, and provenance metadata")
    output_dir = resolve_from_config(config_path, config["output_dir"])
    table_paths = save_tables(
        participants, comparisons, trends, descriptives, output_dir
    )
    pdf_path = resolve_from_config(config_path, config["figure_pdf"])
    png_path = resolve_from_config(config_path, config["figure_png"])
    if not args.skip_spins:
        plot_participant_results(
            participants,
            list(config["group_order"]),
            config["group_display_names"],
            trends[0],
            trends[1],
            pdf_path,
            png_path,
            seed=int(model_config["seed"]),
        )

    metadata = build_metadata(
        config,
        inputs,
        participants,
        n_vertices,
        table_paths,
        pdf_path,
        png_path,
        skip_spins=args.skip_spins,
    )
    metadata_path = save_metadata(metadata, output_dir)

    print_summary(
        descriptives,
        trends,
        comparisons,
        table_paths["participant_scores"],
        metadata_path,
    )


if __name__ == "__main__":
    main()
