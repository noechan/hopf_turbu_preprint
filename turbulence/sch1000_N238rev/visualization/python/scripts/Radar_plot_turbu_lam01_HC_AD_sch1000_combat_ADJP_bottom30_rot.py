#!/usr/bin/env python3
"""Compatibility entry point for the canonical Figure 3 Yeo-7 workflow.

The historical implementation read a superseded ``permFDR.xlsx`` file. This
wrapper deliberately rebuilds both radar panels from the current N145
age/sex/education-adjusted Freedman--Lane results so the two panels cannot
silently diverge.
"""

from __future__ import annotations

from pathlib import Path
import subprocess
import sys


SCH1000_ROOT = Path(__file__).resolve().parents[3]
CANONICAL_SCRIPT = (
    SCH1000_ROOT
    / "statistical_analysis"
    / "nodewise_metastability"
    / "run_yeo7_top30_N145.py"
)


if __name__ == "__main__":
    print(f"Delegating to canonical radar workflow: {CANONICAL_SCRIPT}")
    subprocess.run([sys.executable, str(CANONICAL_SCRIPT)], check=True)
