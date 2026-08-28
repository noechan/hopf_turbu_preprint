#!/usr/bin/env python3
"""Plot a saved AT(N)-style N145 brain-dynamics--MOCA regression."""

from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
import tempfile
from typing import Final

os.environ.setdefault(
    "MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "turbu_n145_matplotlib")
)
os.environ.setdefault("MPLBACKEND", "Agg")

import matplotlib.pyplot as plt
import pandas as pd


PYTHON_DIR: Final = Path(__file__).resolve().parent.parent
SCH1000_ROOT: Final = PYTHON_DIR.parent.parent
PREDICTOR_SPECS: Final = {
    "information-flow": {
        "source": "InfoFlow_lam_0_01",
        "value": "InfoFlow_lam_0_01",
        "slug": "information_flow",
        "display": "Information flow",
        "file_display": "InformationFlow",
        "transform": "identity",
    },
    "turbulence": {
        "source": "Turbu_lam_0_01",
        "value": "Turbu_lam_0_01",
        "slug": "turbulence",
        "display": "Turbulence",
        "file_display": "Turbulence",
        "transform": "identity",
    },
    "one-minus-information-transfer": {
        "source": "InfoTransfer_lam_0_01",
        "value": "OneMinus_InfoTransfer_lam_0_01",
        "slug": "one_minus_information_transfer",
        "display": "1 - Information transfer",
        "file_display": "OneMinusInformationTransfer",
        "transform": "one_minus",
    },
}
DEFAULT_HARMONIZED_FILE: Final = (
    SCH1000_ROOT
    / "harmonization_allfeat"
    / "recomputed"
    / "Turbu_ComBat_ADNI3_allfeatures_N145.xlsx"
)
DEFAULT_CLINICAL_FILE: Final = (
    PYTHON_DIR / "data" / "Turbu_ComBat_clin_ADNI3_allfeatures_N145.csv"
)
DEFAULT_OUTPUT_DIR: Final = (
    SCH1000_ROOT
    / "figures_N145"
    / "sch1000"
    / "Abeta_Status"
    / "harmonized_allfeat"
    / "python"
    / "cognition_moca"
)
GROUP_ORDER: Final = ["HC_ABneg", "HC_ABpos", "MCI_ABpos", "AD_ABpos"]
GROUP_LABELS: Final = {
    "HC_ABneg": r"HC$^-$",
    "HC_ABpos": r"HC$^+$",
    "MCI_ABpos": r"MCI$^+$",
    "AD_ABpos": r"AD$^+$",
}
GROUP_COLORS: Final = {
    "HC_ABneg": "#0173B2",
    "HC_ABpos": "#DE8F05",
    "MCI_ABpos": "#029E73",
    "AD_ABpos": "#D55E00",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--predictor",
        choices=tuple(PREDICTOR_SPECS),
        default="information-flow",
        help="Saved dynamical predictor to plot (default: information-flow)",
    )
    parser.add_argument("--harmonized-file", type=Path, default=DEFAULT_HARMONIZED_FILE)
    parser.add_argument("--clinical-file", type=Path, default=DEFAULT_CLINICAL_FILE)
    parser.add_argument("--results-dir", type=Path, default=None)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    return parser.parse_args()


def require_file(path: Path, label: str) -> None:
    if not path.is_file():
        raise FileNotFoundError(f"Missing {label}: {path}")


def normalize_ptid(frame: pd.DataFrame, label: str) -> pd.DataFrame:
    if "PTID" not in frame.columns:
        raise ValueError(f"{label} is missing PTID")
    result = frame.copy()
    result["PTID"] = result["PTID"].astype("string").str.strip()
    if result["PTID"].duplicated().any():
        raise ValueError(f"{label} contains duplicate PTIDs")
    return result


def main() -> None:
    args = parse_args()
    spec = PREDICTOR_SPECS[args.predictor]
    predictor_source = spec["source"]
    predictor = spec["value"]
    predictor_slug = spec["slug"]
    predictor_display = spec["display"]
    if args.results_dir is None:
        args.results_dir = (
            SCH1000_ROOT
            / "statistical_analysis"
            / "cognition_moca"
            / "results"
            / f"N145_MOCA_{predictor_slug}_regression"
        )
    file_stem = f"N145_MOCA_{predictor_slug}_HC3_age_sex_education"
    result_path = args.results_dir / f"{file_stem}.csv"
    prediction_path = (
        args.results_dir
        / f"{file_stem}_prediction_line.csv"
    )
    for path, label in (
        (args.harmonized_file, "harmonized input"),
        (args.clinical_file, "clinical input"),
        (result_path, "saved regression result"),
        (prediction_path, "saved prediction line"),
    ):
        require_file(path, label)

    harmonized = normalize_ptid(pd.read_excel(args.harmonized_file), "harmonized input")
    clinical = normalize_ptid(pd.read_csv(args.clinical_file), "clinical input")
    required_harmonized = {"PTID", "Group", predictor_source}
    required_clinical = {"PTID", "Group", "MOCA"}
    if not required_harmonized.issubset(harmonized.columns):
        raise ValueError("Harmonized input lacks the required plotting columns")
    if not required_clinical.issubset(clinical.columns):
        raise ValueError("Clinical input lacks the required plotting columns")
    frame = harmonized[["PTID", "Group", predictor_source]].merge(
        clinical[["PTID", "Group", "MOCA"]].rename(columns={"Group": "Clinical_Group"}),
        on="PTID",
        how="left",
        validate="one_to_one",
    )
    if not frame["Group"].equals(frame["Clinical_Group"]):
        raise ValueError("Group labels disagree between harmonized and clinical inputs")
    frame[predictor_source] = pd.to_numeric(
        frame[predictor_source], errors="raise"
    )
    if spec["transform"] == "one_minus":
        frame[predictor] = 1 - frame[predictor_source]
    else:
        frame[predictor] = frame[predictor_source]
    frame["MOCA"] = pd.to_numeric(frame["MOCA"], errors="coerce")
    frame = frame.dropna(subset=[predictor, "MOCA"])
    if len(frame) != 144:
        raise ValueError(f"Expected 144 participants with MOCA; found {len(frame)}")

    result = pd.read_csv(result_path).iloc[0]
    prediction = pd.read_csv(prediction_path)
    args.output_dir.mkdir(parents=True, exist_ok=True)

    fig, ax = plt.subplots(figsize=(8.0, 6.5))
    for group in GROUP_ORDER:
        group_frame = frame.loc[frame["Group"] == group]
        ax.scatter(
            group_frame[predictor],
            group_frame["MOCA"],
            s=48,
            alpha=0.72,
            color=GROUP_COLORS[group],
            edgecolor="white",
            linewidth=0.5,
            label=f"{GROUP_LABELS[group]} (N={len(group_frame)})",
        )
    ax.plot(
        prediction[predictor],
        prediction["Predicted_MOCA"],
        color="black",
        linewidth=2.5,
        zorder=5,
    )
    annotation = (
        rf"$\beta={result['Beta_MOCA_per_predictor_SD']:.2f}$"
        "\n"
        rf"$R^2={result['R2_Full_Model']:.3f}$"
    )
    ax.text(
        0.03,
        0.97,
        annotation,
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=12,
    )
    ax.set_xlabel(rf"{predictor_display} ($\lambda=0.01$)", fontsize=14)
    ax.set_ylabel("MOCA", fontsize=14)
    ax.set_title(f"{predictor_display} and cognitive performance", fontsize=15)
    ax.grid(axis="y", alpha=0.22)
    ax.spines[["top", "right"]].set_visible(False)
    ax.legend(frameon=False, fontsize=10, loc="lower right")
    fig.tight_layout()

    stem = (
        f"MOCA_vs_{spec['file_display']}_lambda_0_01_HC3_age_sex_education"
    )
    pdf_path = args.output_dir / f"{stem}.pdf"
    png_path = args.output_dir / f"{stem}.png"
    fig.savefig(pdf_path, bbox_inches="tight")
    fig.savefig(png_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    (args.output_dir / f"{stem}_provenance.json").write_text(
        json.dumps(
            {
                "analysis_result": str(result_path.resolve()),
                "prediction_line": str(prediction_path.resolve()),
                "harmonized_input": str(args.harmonized_file.resolve()),
                "clinical_input": str(args.clinical_file.resolve()),
                "group_role": "scatter-point colour only; not a regression covariate",
                "outputs": [str(pdf_path.resolve()), str(png_path.resolve())],
            },
            indent=2,
        )
        + "\n",
        encoding="utf-8",
    )
    print(f"Saved: {pdf_path.resolve()}")
    print(f"Saved: {png_path.resolve()}")


if __name__ == "__main__":
    main()
