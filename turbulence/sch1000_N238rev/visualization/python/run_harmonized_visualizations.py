"""Rebuild Python figures from the canonical N145 ComBat outputs.

The global stage reads the harmonized all-feature workbook directly. The
radar stage reads the current age-, sex-, and education-adjusted node-wise
Freedman--Lane results. The AT(N) stage first runs both saved turbulence
regression specifications and then creates separate figures for each model.
The subject-atn stage fits and plots the age-, sex-, and education-adjusted
information-capability and susceptibility associations with AT(N) biomarkers.
The cognition stage fits and plots the corresponding empirical-dynamics--MOCA
regressions. Neurosynth and AHBA analyses are maintained as the separate,
repository-relative workflow under `neuromaps_analysis/`.
"""

from __future__ import annotations

import argparse
import os
from pathlib import Path
import subprocess
import sys
import tempfile


PYTHON_DIR = Path(__file__).resolve().parent
SCRIPTS_DIR = PYTHON_DIR / "scripts"
SCH1000_ROOT = PYTHON_DIR.parent.parent
HARMONIZED_INPUT = (
    SCH1000_ROOT
    / "harmonization_allfeat"
    / "recomputed"
    / "Turbu_ComBat_ADNI3_allfeatures_N145.xlsx"
)
NODEWISE_ANALYSIS_DIR = (
    SCH1000_ROOT / "statistical_analysis" / "nodewise_metastability"
)
NODEWISE_RESULTS_DIR = (
    NODEWISE_ANALYSIS_DIR / "results" / "N145_nodewise_metastability"
)

GLOBAL_SCRIPTS = [
    "InfoFlow_combat_allfeat.py",
    "turbu_lambdas_range_combat_allfeat.py",
    "InfoTransfer_combat_allfeat.py",
]

RADAR_SCRIPT = NODEWISE_ANALYSIS_DIR / "run_yeo7_top30_N145.py"
ATN_ANALYSIS_SCRIPT = (
    SCH1000_ROOT
    / "statistical_analysis"
    / "atn_biomarkers"
    / "run_atn_regressions_N145.py"
)
ATN_PLOT_SCRIPT = SCRIPTS_DIR / "plot_atn_regressions.py"
SUBJECT_ATN_ANALYSIS_SCRIPT = (
    SCH1000_ROOT
    / "statistical_analysis"
    / "perturbation_measures"
    / "run_atn_subjectlevel_perturbation_N145.py"
)
SUBJECT_ATN_PLOT_SCRIPT = (
    SCRIPTS_DIR / "plot_atn_subjectlevel_perturbation.py"
)
SUBJECT_ATN_SHAPE_ANALYSIS_SCRIPT = (
    SCH1000_ROOT
    / "statistical_analysis"
    / "perturbation_measures"
    / "run_atn_subjectlevel_shape_sensitivity_N145.py"
)
SUBJECT_ATN_SHAPE_PLOT_SCRIPT = (
    SCRIPTS_DIR / "plot_atn_subjectlevel_shape_sensitivity.py"
)
COGNITION_ANALYSIS_SCRIPT = (
    SCH1000_ROOT
    / "statistical_analysis"
    / "cognition_moca"
    / "run_moca_information_flow_regression_N145.py"
)
COGNITION_PLOT_SCRIPT = SCRIPTS_DIR / "plot_moca_information_flow_regression.py"
COGNITION_COMPARISON_SCRIPT = (
    SCH1000_ROOT
    / "statistical_analysis"
    / "cognition_moca"
    / "compare_moca_single_predictor_regressions_N145.py"
)
COGNITION_LMG_SCRIPT = (
    SCH1000_ROOT
    / "statistical_analysis"
    / "cognition_moca"
    / "run_moca_lmg_relative_importance_N145.py"
)
RADAR_INPUTS = [
    NODEWISE_RESULTS_DIR
    / "N145_nodewise_metastability_HC_ABneg_vs_AD_ABpos_"
    "FreedmanLane_age_sex_education.csv",
    NODEWISE_RESULTS_DIR
    / "N145_nodewise_metastability_MCI_ABpos_vs_AD_ABpos_"
    "FreedmanLane_age_sex_education.csv",
]


def require_files(paths: list[Path], label: str) -> None:
    missing = [path for path in paths if not path.is_file()]
    if missing:
        formatted = "\n".join(f"  - {path}" for path in missing)
        raise FileNotFoundError(f"Missing {label}:\n{formatted}")


def run_script(script: Path, arguments: list[str] | None = None) -> None:
    env = os.environ.copy()
    env.setdefault("MPLBACKEND", "Agg")
    env.setdefault("PYTHONDONTWRITEBYTECODE", "1")
    matplotlib_cache = Path(tempfile.gettempdir()) / "turbu_n145_matplotlib"
    matplotlib_cache.mkdir(parents=True, exist_ok=True)
    env.setdefault("MPLCONFIGDIR", str(matplotlib_cache))
    require_files([script], "Python pipeline script")
    arguments = arguments or []
    print(f"\nRunning {script.name} {' '.join(arguments)}", flush=True)
    subprocess.run(
        [sys.executable, str(script), *arguments], check=True, env=env
    )


def run_scripts(names: list[str]) -> None:
    for name in names:
        run_script(SCRIPTS_DIR / name)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--stage",
        choices=(
            "global",
            "radar",
            "atn",
            "subject-atn",
            "subject-atn-shape",
            "cognition",
            "all",
        ),
        default="global",
        help=(
            "global reads the N145 all-feature workbook; radar reads fresh "
            "age/sex/education-adjusted node-wise outputs; atn rebuilds "
            "Figure 4 c/d/j/k/l; subject-atn rebuilds the education-adjusted "
            "subject-level Hopf--AT(N) models and plots; subject-atn-shape "
            "runs the Figure 4 within-group and pooled-quadratic sensitivity "
            "analysis; cognition fits and "
            "plots the MOCA models; "
            "all runs every stage (default: global)"
        ),
    )
    args = parser.parse_args()

    if args.stage in {"global", "all"}:
        require_files([HARMONIZED_INPUT], "harmonized N145 input")
        run_scripts(GLOBAL_SCRIPTS)

    if args.stage in {"radar", "all"}:
        require_files(RADAR_INPUTS, "covariate-adjusted node-wise result")
        require_files([RADAR_SCRIPT], "Yeo-7 top-30% script")
        run_script(RADAR_SCRIPT)

    if args.stage in {"atn", "all"}:
        require_files([HARMONIZED_INPUT], "harmonized N145 input")
        run_script(ATN_ANALYSIS_SCRIPT, ["--model", "both"])
        run_script(ATN_PLOT_SCRIPT, ["--model", "age_sex"])
        run_script(ATN_PLOT_SCRIPT, ["--model", "age_sex_education"])

    if args.stage in {"subject-atn", "all"}:
        run_script(SUBJECT_ATN_ANALYSIS_SCRIPT)
        run_script(SUBJECT_ATN_PLOT_SCRIPT)

    if args.stage in {"subject-atn-shape", "all"}:
        run_script(SUBJECT_ATN_SHAPE_ANALYSIS_SCRIPT)
        run_script(SUBJECT_ATN_SHAPE_PLOT_SCRIPT)

    if args.stage in {"cognition", "all"}:
        require_files([HARMONIZED_INPUT], "harmonized N145 input")
        run_script(COGNITION_ANALYSIS_SCRIPT, ["--predictor", "information-flow"])
        run_script(COGNITION_PLOT_SCRIPT, ["--predictor", "information-flow"])
        run_script(COGNITION_ANALYSIS_SCRIPT, ["--predictor", "turbulence"])
        run_script(COGNITION_PLOT_SCRIPT, ["--predictor", "turbulence"])
        run_script(
            COGNITION_ANALYSIS_SCRIPT,
            ["--predictor", "one-minus-information-transfer"],
        )
        run_script(
            COGNITION_PLOT_SCRIPT,
            ["--predictor", "one-minus-information-transfer"],
        )
        run_script(COGNITION_COMPARISON_SCRIPT)
        run_script(COGNITION_LMG_SCRIPT)

    print("\nRequested Python visualization stage completed.")


if __name__ == "__main__":
    main()
