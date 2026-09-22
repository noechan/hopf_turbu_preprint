#!/usr/bin/env python3
"""Render the participant-level coupling equations for presentation slides."""

from __future__ import annotations

import os
from pathlib import Path
import tempfile


os.environ.setdefault(
    "MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "neuromaps_matplotlib")
)

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


ROOT = Path(__file__).resolve().parents[2]
OUTPUT_DIR = ROOT / "output" / "equations"

EQUATIONS = {
    "participant_correlation_ri": (
        r"$r_i = \mathrm{corr}\!\left("
        r"\mathrm{E\!:\!I}_{1,\ldots,V},\,"
        r"\mathrm{turbulence}_{i,1,\ldots,V}"
        r"\right)$"
    ),
    "spatially_normalized_zspin": (
        r"$z_{\mathrm{spin},i} = "
        r"\dfrac{r_{i,\mathrm{observed}} - "
        r"\mathrm{mean}\!\left(r_{i,\mathrm{rotated}}\right)}"
        r"{\mathrm{SD}\!\left(r_{i,\mathrm{rotated}}\right)}$"
    ),
    "ordered_stage_zspin_model": (
        r"$z_{\mathrm{spin}} = \beta_0 + \beta_1\,\mathrm{stage} "
        r"+ \beta_2\,\mathrm{age} + \beta_3\,\mathrm{gender} "
        r"+ \beta_4\,\mathrm{education}$"
    ),
}


def render_equation(stem: str, equation: str) -> None:
    figure = plt.figure(figsize=(12, 1.7), facecolor="none")
    figure.text(
        0.5,
        0.5,
        equation,
        ha="center",
        va="center",
        fontsize=32,
        color="#111111",
    )
    png_path = OUTPUT_DIR / f"{stem}.png"
    svg_path = OUTPUT_DIR / f"{stem}.svg"
    figure.savefig(
        png_path,
        dpi=300,
        transparent=True,
        bbox_inches="tight",
        pad_inches=0.10,
    )
    figure.savefig(
        svg_path,
        transparent=True,
        bbox_inches="tight",
        pad_inches=0.10,
    )
    plt.close(figure)
    print(f"Saved: {png_path}")
    print(f"Saved: {svg_path}")


def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    plt.rcParams.update(
        {
            "font.family": "STIXGeneral",
            "mathtext.fontset": "stix",
            "svg.fonttype": "path",
        }
    )
    for stem, equation in EQUATIONS.items():
        render_equation(stem, equation)


if __name__ == "__main__":
    main()
