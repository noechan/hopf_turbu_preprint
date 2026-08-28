# Data exports

This directory contains final data-export scripts and their retained manuscript
outputs. It is separate from `prepare_data/`, which ends after construction of
the diagnosis- and amyloid-specific time-series inputs.

- The top-level MATLAB scripts export ML-ready and parcel-level turbulence data.
- `data_for_ML_4Staging.m` reads the four standalone amyloid-status turbulence
  outputs directly and exports 32 non-constant features for N152.
- `export_harmonization_inputs.m` writes the N152 all-feature and node-level
  raw tables consumed by `harmonization_allfeat/run_harmonization.py`. It uses
  PTIDs embedded in the standalone result files and therefore does not depend
  on diagnosis-wide column ordering. It exports physical lambda values 0.01,
  0.03 and 0.06 for the Neuromaps analysis (0.01 is also used by manuscript
  Figure 4). Historical `lam1`, `lam3` and `lam6` filenames denoted the
  physical values multiplied by 100, not MATLAB array indices.
- `neuromaps_exports/` creates Schaefer-1000 group-difference maps from the
  harmonized tables and writes them to the configured neuromaps annotation
  directory on the mounted ADNI drive.
