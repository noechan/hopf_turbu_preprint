#!/usr/bin/env python3
"""Export the spatially normalized E:I-turbulence coupling panel alone."""

from __future__ import annotations

import json
import os
from pathlib import Path
import tempfile

os.environ.setdefault(
    "MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "neuromaps_matplotlib")
)
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from plot_ei_surface_spin_coupling_lam001 import (
    CONFIG_PATH,
    COMPARISON_PATH,
    PARTICIPANT_PATH,
    TREND_PATH,
    load_inputs,
    render_spin_panel,
    sha256,
)


NEUROMAPS_ROOT = Path(__file__).resolve().parents[2]
PDF_PATH = (
    NEUROMAPS_ROOT
    / "output"
    / "pdf"
    / "Figure_spatially_normalized_EI_turbulence_coupling_lam001.pdf"
)
PNG_PATH = (
    NEUROMAPS_ROOT
    / "output"
    / "figures"
    / "Figure_spatially_normalized_EI_turbulence_coupling_lam001.png"
)
METADATA_PATH = (
    NEUROMAPS_ROOT
    / "output"
    / "figures"
    / "Figure_spatially_normalized_EI_turbulence_coupling_lam001_metadata.json"
)


def main() -> None:
    _, participants, trend, ad_comparison = load_inputs()
    plt.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
            "font.size": 9,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "axes.grid": False,
        }
    )
    figure, axis = plt.subplots(figsize=(5.7, 5.5), facecolor="white")
    figure.subplots_adjust(left=0.16, right=0.98, bottom=0.16, top=0.75)
    render_spin_panel(axis, participants, trend, seed=1234)

    PDF_PATH.parent.mkdir(parents=True, exist_ok=True)
    PNG_PATH.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(PDF_PATH, bbox_inches="tight", facecolor="white")
    figure.savefig(PNG_PATH, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(figure)

    metadata = {
        "figure_definition": (
            "Standalone participant-level spatially normalized E:I-turbulence "
            "coupling at physical lambda=0.01"
        ),
        "panel_letter_included": False,
        "inputs": {
            "configuration": {"path": str(CONFIG_PATH), "sha256": sha256(CONFIG_PATH)},
            "participant_scores": {
                "path": str(PARTICIPANT_PATH),
                "sha256": sha256(PARTICIPANT_PATH),
            },
            "ordered_trends": {"path": str(TREND_PATH), "sha256": sha256(TREND_PATH)},
            "planned_comparisons": {
                "path": str(COMPARISON_PATH),
                "sha256": sha256(COMPARISON_PATH),
            },
        },
        "outcome": "spin_z",
        "ordered_trend": trend,
        "hc_negative_vs_ad_positive": ad_comparison,
        "outputs": {"pdf": str(PDF_PATH), "png": str(PNG_PATH)},
    }
    with METADATA_PATH.open("w", encoding="utf-8") as stream:
        json.dump(metadata, stream, indent=2)
        stream.write("\n")

    print(f"Saved PDF: {PDF_PATH}")
    print(f"Saved PNG: {PNG_PATH}")
    print(f"Saved metadata: {METADATA_PATH}")


if __name__ == "__main__":
    main()
