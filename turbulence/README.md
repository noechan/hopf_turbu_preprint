# Turbulence full-denoising analyses

This directory contains the publication version of the empirical turbulence
pipeline. Historical, exploratory, generated, and atlas-mismatched branches
from the working repository are intentionally omitted.

## Cohort preservation policy

Scripts must remain separated when they represent different diagnostic or
amyloid-status cohorts. In particular, HC Aβ−, HC Aβ+, MCI Aβ+, and AD Aβ+
analyses are scientifically distinct even when their implementations are very
similar. Do not remove or merge a cohort-specific script based only on filename
or code similarity. Any future shared implementation must retain explicit
cohort wrappers and must verify PTID membership and subject counts.

## Layout

- `sch1000_N238rev/`: current Schaefer-1000 N238-revision pipeline, organized
  into preparation, turbulence calculation, harmonization, statistical
  analysis, data export, and visualization stages.
- `helper_functions/`: shared MATLAB utilities.

Set `ADNI3_ROOT` for restricted inputs. Visible `/path/to/...` strings are
deliberate placeholders for institution-specific Hopf or imaging resources;
no personal filesystem paths are retained.

## Python environment

The Python environment is captured in:

- `requirements.txt`: direct packages used by the current Python scripts.


To create a fresh environment:

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install -r requirements.txt
```

Use Python 3.12 and validate all inputs before accepting regenerated outputs.

## Data and generated results

Participant-level data and generated results are excluded from this release.
See the repository-level `DATA.md` for the expected local inputs.
