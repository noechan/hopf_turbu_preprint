#!/usr/bin/env python3
"""Convenience entry point for age-, sex-, and education-adjusted AT(N) models."""

from __future__ import annotations

import sys

from run_atn_regressions_N145 import main


if __name__ == "__main__":
    main(["--model", "age_sex_education", *sys.argv[1:]])
