"""Write participant-level analysis tables and reproducibility metadata."""

from __future__ import annotations

import json
from pathlib import Path

import neuromaps
import numpy as np
import pandas as pd

from .inputs import AnalysisInputs, file_sha256


def save_tables(
    participants: pd.DataFrame,
    comparisons: list[dict[str, object]],
    trends: list[dict[str, object]],
    descriptives: pd.DataFrame,
    output_dir: Path,
) -> dict[str, Path]:
    """Write the four tabular analysis outputs."""

    output_dir.mkdir(parents=True, exist_ok=True)
    paths = {
        "participant_scores": output_dir
        / "participant_EI_turbulence_coupling_lam001_N145.csv",
        "group_descriptives": output_dir
        / "group_descriptives_EI_turbulence_lam001_N145.csv",
        "planned_comparisons": output_dir
        / "confirmatory_tests_EI_turbulence_lam001_N145.csv",
        "ordered_trends": output_dir
        / "ordered_trend_EI_turbulence_lam001_N145.csv",
    }
    participants.to_csv(paths["participant_scores"], index=False)
    descriptives.to_csv(paths["group_descriptives"], index=False)
    pd.DataFrame(comparisons).to_csv(paths["planned_comparisons"], index=False)
    pd.DataFrame(trends).to_csv(paths["ordered_trends"], index=False)
    return paths


def build_metadata(
    config: dict,
    inputs: AnalysisInputs,
    participants: pd.DataFrame,
    n_vertices: int,
    table_paths: dict[str, Path],
    pdf_path: Path,
    png_path: Path,
    skip_spins: bool,
) -> dict:
    """Assemble a reproducibility record without changing the analysis."""

    paths = inputs.paths
    return {
        "analysis_name": config["analysis_name"],
        "analysis_definition": config["analysis_definition"],
        "configuration": config,
        "cohort": {
            "harmonized_n": len(inputs.harmonized_nodes),
            "final_n": len(participants),
            "final_group_counts": participants["Group"].value_counts().to_dict(),
        },
        "validation": {
            "reference_hc_ad_map": str(paths.validation_reference_hc_ad),
            "max_absolute_error": inputs.validation_max_error,
        },
        "inputs": {
            "harmonized_node_table": {
                "path": str(paths.harmonized_node_table),
                "sha256": file_sha256(paths.harmonized_node_table),
            },
            "harmonization_script": {
                "path": str(paths.harmonization_script),
                "sha256": file_sha256(paths.harmonization_script),
            },
            "demographics": {
                "path": str(paths.demographics),
                "sha256": file_sha256(paths.demographics),
            },
            "n145_cohort_table": {
                "path": str(paths.n145_cohort_table),
                "sha256": file_sha256(paths.n145_cohort_table),
            },
            "ei_values": {
                "path": str(paths.ei_values),
                "sha256": file_sha256(paths.ei_values),
            },
            "ei_nifti": {
                "path": str(paths.ei_nifti),
                "sha256": file_sha256(paths.ei_nifti),
            },
        },
        "surface_vertices_used": n_vertices,
        "statistical_hemisphere": config["surface"]["hemisphere"],
        "source_parcels_used": int(config["surface"]["source_parcels"]),
        "software": {
            "neuromaps": getattr(neuromaps, "__version__", "unknown"),
            "numpy": np.__version__,
            "pandas": pd.__version__,
        },
        "outputs": {
            **{name: str(path) for name, path in table_paths.items()},
            "figure_pdf": str(pdf_path) if not skip_spins else None,
            "figure_png": str(png_path) if not skip_spins else None,
        },
    }


def save_metadata(metadata: dict, output_dir: Path) -> Path:
    """Write the JSON provenance record."""

    metadata_path = output_dir / "subjectlevel_EI_turbulence_lam001_metadata.json"
    with metadata_path.open("w", encoding="utf-8") as stream:
        json.dump(metadata, stream, indent=2)
        stream.write("\n")
    return metadata_path


def print_summary(
    descriptives: pd.DataFrame,
    trends: list[dict[str, object]],
    comparisons: list[dict[str, object]],
    participant_path: Path,
    metadata_path: Path,
) -> None:
    """Print a compact terminal summary after a successful run."""

    print("\nGroup descriptives:")
    print(
        descriptives[["group", "sample_size", "mean_r", "sd_r"]].to_string(
            index=False
        )
    )
    print("\nOrdered trends:")
    print(pd.DataFrame(trends).to_string(index=False))
    print("\nPlanned comparisons:")
    print(pd.DataFrame(comparisons).to_string(index=False))
    print(f"Saved participant scores: {participant_path}")
    print(f"Saved metadata: {metadata_path}")
