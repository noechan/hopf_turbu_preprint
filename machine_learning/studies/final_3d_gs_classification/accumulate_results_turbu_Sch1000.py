"""Accumulate turbulence-only LinearSVM and PolySVM results."""

from __future__ import annotations

import argparse
from pathlib import Path
import sys


PATH_REPO = Path(__file__).resolve().parents[2]
if str(PATH_REPO) not in sys.path:
    sys.path.insert(0, str(PATH_REPO))

from src.utils.result_accumulation import accumulate_results


RESULTS_DIR = PATH_REPO / "Results" / "final_3d_gs_classification_turbu_sch1000"
FEATURE_LABELS = {"all_features_turbu_combat": "Turbu"}
GROUP_LABELS = {
    "HCneg_vs_HCpos": "HC (AB-) vs HC (AB+)",
    "HCneg_vs_MCIpos": "HC (AB-) vs MCI (AB+)",
    "HCneg_vs_ADpos": "HC (AB-) vs AD (AB+)",
    "MCIpos_vs_ADpos": "MCI (AB+) vs AD (AB+)",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--check-only",
        action="store_true",
        help="Validate and summarize results without writing the Excel workbook.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    results, output_path = accumulate_results(
        results_dir=RESULTS_DIR,
        feature_labels=FEATURE_LABELS,
        group_labels=GROUP_LABELS,
        required_classifiers=("LinearSVM", "PolySVM"),
        expected_folds=40,
        check_only=args.check_only,
    )

    preview_columns = [
        "Classifier",
        "Group Comparison ID",
        "Feature Set",
        "Balanced Accuracy (Test) Mean",
        "AUC (Test) Mean",
    ]
    print(results[preview_columns].to_string(index=False))
    print(f"Validated {len(results)} classifier/comparison/feature result rows.")
    if output_path is None:
        print("Check complete; no workbook was written.")
    else:
        print(f"Saved accumulated results to: {output_path}")


if __name__ == "__main__":
    main()
