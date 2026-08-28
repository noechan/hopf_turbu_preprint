# Schaefer-1000 N238rev pipeline

This directory contains the publication Schaefer-1000 workflow. Code is
resolved from this repository; restricted participant and time-series data are
resolved from an authorized local ADNI directory.

## Configuration

`pipeline_paths.m` derives all repository paths from its own location. The
external data root defaults to:

```text
/path/to/ADNI3
```

To use another mount without editing scripts, set `ADNI3_ROOT` before starting
MATLAB.

## Environment and read-only tests

Use MATLAB R2026a. The R2025b installation does not contain the Signal
Processing and Statistics toolboxes needed by this pipeline.

```matlab
addpath('/path/to/hopf_turbu_preprint/turbulence/sch1000_N238rev', '-begin')
preflight_sch1000()
validate_sch1000_data(false)
validate_sch1000_data(true)  % reads and checks every 1000 x 197 time series
```

These tests do not modify any data or results.

## Python environment

The repository-level virtual environment contains the pinned visualization and
harmonization dependencies:

```bash
cd /path/to/hopf_turbu_preprint
.venv-turbulence/bin/python -m pip install -r turbulence/requirements-lock.txt
```

For non-interactive plot checks, set `MPLBACKEND=Agg` so that scripts do not
wait for a display window.

## Cohorts preserved by the pipeline

The full selected sample has 238 participants: 109 HC, 90 MCI, and 39 AD.
Amyloid-specific downstream analyses keep HC Aβ−, HC Aβ+, MCI Aβ+, and AD
Aβ+ as separate groups. Do not combine these scripts solely because their
code structure is similar.

## Main execution order

1. `prepare_data/01_cohort_selection/s1_*`: select QC-passing batch PTIDs and
   CONN IDs.
2. `prepare_data/02_ptid_assembly/s2_*`: concatenate PTIDs in the same batch
   order used for the time series and export the combined participant table.
3. `prepare_data/03_timeseries_extraction/s3_*`: extract every supported
   parcellation, including Schaefer-1000, from the extended CONN outputs.
4. `prepare_data/04_combine_modalities/s4_*`: combine MPRAGE and IR-FSPGR time
   series.
5. `prepare_data/05_amyloid_groups/s5_*`: derive amyloid-specific time-series
   and PTID groups across all parcellations.
6. `calculate_turbu/amyloid_status/*`: calculate the manuscript HC Aβ−,
   HC Aβ+, MCI Aβ+, and AD Aβ+ turbulence observables directly from the
   time series created in stage 5.
7. `data_export/export_harmonization_inputs.m`: export the N152 all-feature and
   node-level raw tables for ComBat, preserving the PTIDs embedded in each
   result file.
8. `harmonization_allfeat/run_harmonization.py`: validate the raw tables and
   recompute the N145 ComBat tables.
9. `data_export/data_for_ML_4Staging.m`: validate the recomputed N145
   all-feature ComBat table and export ML-ready MAT, CSV, and XLSX copies with
   the original `PTID`, `Group`, and 32-feature layout.
   `data_export/data_for_ML_4Staging_with_InfoCap_Susceptibility.m` performs
   the corresponding PTID-matched combined export, placing `Info_Cap` and
   `Susceptibility` before the same 32 recomputed dynamics features.
10. `statistical_analysis/turbu_harm_stats/run_group_permutations_N145_age_sex_education.R`:
    run the manuscript Freedman--Lane group comparisons on the N145 ComBat
    turbulence table, with age, sex, and education as nuisance covariates.
11. Statistics and visualization.

The canonical harmonized-result path is `P.harmonization_results`, currently
`harmonization_allfeat/recomputed/`. This generated directory is ignored by
Git and must be recreated locally.

## Recomputing harmonization from the mounted drive

First run `data_export/export_harmonization_inputs.m` in MATLAB. It writes raw
tables derived directly from the standalone amyloid-status results to
`$ADNI3_ROOT/timeseries/harmonization_inputs/sch1000_N238rev`.

`harmonization_allfeat/run_harmonization.py` reads those feature tables and the
separate covariate/site metadata from the historical neuroHarmonize project.
By default it writes to `harmonization_allfeat/recomputed/`, leaving the active
manuscript tables unchanged:

```bash
cd /path/to/hopf_turbu_preprint
.venv-turbulence/bin/python \
  turbulence/sch1000_N238rev/harmonization_allfeat/run_harmonization.py \
  --exclusions /path/to/authorized/exclusions.csv \
  --validate-only

.venv-turbulence/bin/python \
  turbulence/sch1000_N238rev/harmonization_allfeat/run_harmonization.py \
  --exclusions /path/to/authorized/exclusions.csv
```

The script derives the six single-site MRI exclusions from the input data,
fits ComBat on the remaining 146 participants, and then applies the authorized
local manuscript exclusion. This post-ComBat order is required to reproduce
the retained N145 values.

The post-ComBat exclusion was made because the participant lacked the final
VBM/GMV output required for the common multimodal manuscript cohort. To comply
with the ADNI Data Use Agreement, the participant identifier and exclusion
manifest are not distributed here. Authorized users must provide the local
manifest with `--exclusions /path/to/exclusions.csv`.

Generated participant metadata are stored at
`$ADNI3_ROOT/timeseries/pipeline_metadata/sch1000_N238rev`, not in the code
repository.
