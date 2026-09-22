# Neuromaps analyses

This directory contains the publication workflows for spatially relating the
N145 Schaefer-1000 turbulence difference maps to Neurosynth memory maps and to
15 a priori Alzheimer-related AHBA gene-expression maps. All scripts resolve
repository files relative to their own location.

The script suffixes `lam1`, `lam3`, and `lam6` denote the physical
values $\lambda=0.01$, $0.03$, and $0.06$, respectively. The turbulence maps
are absolute group-mean differences, consistent with the manuscript analysis.

## Environment

```bash
python3.12 -m venv .venv
.venv/bin/python -m pip install --upgrade pip
.venv/bin/python -m pip install -r requirements.txt
```

## Rebuild the aggregate turbulence NIfTI maps

The six non-identifiable parcel-level `.mat` inputs and the Schaefer atlas are
included. Rebuild their MNI152 NIfTI representations with:

```bash
.venv/bin/python \
  neuromaps_code/parcellate/build_turbulence_difference_niftis.py
```


## AHBA gene-expression correlations

The processed gene maps are not redistributed. Set `AHBA_GENE_NIFTI_DIR` to a
directory containing files named `<GENE>_expression_sch1000_2mm.nii.gz`, then
run the desired contrast and scale, for example:

```bash
AHBA_GENE_NIFTI_DIR=/path/to/authorized/gene_niftis \
  .venv/bin/python \
  neuromaps_code/correlations_nulls/AHBA/corr_geneexpr_Turbu_GE15_HC_AD_lam1_fsaverage41k_sch1000.py
```

The six scripts cover HC$^-$ vs. AD$^+$ and MCI$^+$ vs. AD$^+$ at the three
scales. Each uses fsaverage 41k projection, Pearson correlation, and 1,000
Alexander--Bloch rotations with seed 1234. Outputs are written to
`results/gene_spatial_correlations/` and ignored by Git.

## Primary transcriptomic E:I correlation with HC turbulence at lambda = 0.01

The primary E:I map is derived from the same processed Schaefer-1000 abagen
matrix used for the 15-gene analysis. Its gene sets and source-processing
provenance are recorded in `../abagen_analysis/config/ei_sch1000.json`.

From the repository root, build the parcel table and raw/normalized NIfTI maps:

```bash
.venv-neuromaps/bin/python \
  abagen_analysis/abagen-code/build_ei_sch1000.py
```

Visualize the raw and min--max-normalized E:I maps as axial Nilearn slices:

```bash
.venv-neuromaps/bin/python \
  abagen_analysis/abagen-code/plot_ei_sch1000.py
```

The combined PNG/PDF and separate map PNGs are written to
`abagen_analysis/abagen-code/ei_maps/schaefer1000_2mm/figures/`. Use `--force`
to replace figures that have already been generated.

For the six-view MATLAB surface rendering used by the turbulence figures,
the included 32k fsLR/Schaefer GIFTI assets are resolved through
`setup_sch1000_paths.m`. Run:

```matlab
addpath('turbulence/sch1000_N238rev/visualization/render')
render_ei_sch1000
```

This reads `EI_values_sch1000.csv` directly and writes raw and min--max PNG,
PDF, and editable FIG outputs under
`abagen_analysis/abagen-code/ei_maps/schaefer1000_2mm/figures/surface_matlab/`.

Build the raw HC amyloid-negative group-mean turbulence map at
$\lambda=0.01$ (54 participants):

```bash
.venv-neuromaps/bin/python \
  neuromaps_analysis/neuromaps_code/parcellate/build_hc_turbulence_map.py
```

The default source is the mounted ADNI MATLAB v7.3 file. On another machine,
set `HC_TURBULENCE_MAT` to its local path. The builder verifies that MATLAB
lambda index 10 is the physical value 0.01 before averaging subjects.

Then correlate the E:I map directly with that HC mean map:

```bash
.venv-neuromaps/bin/python \
  neuromaps_analysis/neuromaps_code/correlations_nulls/AHBA/\
corr_EI_Turbu_HC_lam001_fsaverage41k_sch1000.py
```

The analysis uses fsaverage 41k projection, Pearson correlation, and 1,000
Alexander--Bloch rotations with seed 1234, matching the existing 15-gene
workflow. Configuration is stored in
`config/ei_hc_turbulence_lam001.json`; results are written to
`results/ei_hc_turbulence_lam001/` and ignored by Git. This is an HC-only
spatial association, not an HC--AD or MCI--AD group-difference analysis.

The earlier group-difference E:I workflow remains available in
`config/ei_turbulence_lam001.json` and
`corr_EI_Turbu_lam001_fsaverage41k_sch1000.py` for the separate question of
whether E:I follows a disease-related spatial contrast.

## E:I correlation across amyloid-staging groups

To compare the HC amyloid-negative result with the other raw group means at
the same physical scale, build all four maps and run the shared-spin analysis:

```bash
.venv-neuromaps/bin/python \
  neuromaps_analysis/neuromaps_code/parcellate/build_group_turbulence_maps.py

XDG_CACHE_HOME=/private/tmp/codex-neuromaps-cache \
  .venv-neuromaps/bin/python \
  neuromaps_analysis/neuromaps_code/correlations_nulls/AHBA/\
corr_EI_Turbu_groups_lam001_fsaverage41k_sch1000.py
```

The four tests use one shared set of 1,000 Alexander--Bloch rotations. The
summary reports the unadjusted spin-test p-values and Benjamini--Hochberg FDR
values across HC amyloid-negative, HC amyloid-positive, MCI amyloid-positive,
and AD amyloid-positive groups. Configuration and portable source-path
overrides are in `config/ei_group_turbulence_lam001.json`.

## E:I correlation with HC amyloid-negative baseline contrasts

For the disease-focused analysis, compute signed spatial contrasts as the HC
amyloid-negative mean minus each comparison-group mean, then correlate those
contrast maps with E:I:

```bash
.venv-neuromaps/bin/python \
  neuromaps_analysis/neuromaps_code/parcellate/\
build_hc_baseline_turbulence_contrasts.py

XDG_CACHE_HOME=/private/tmp/codex-neuromaps-cache \
  .venv-neuromaps/bin/python \
  neuromaps_analysis/neuromaps_code/correlations_nulls/AHBA/\
corr_EI_Turbu_HC_baseline_contrasts_lam001_fsaverage41k_sch1000.py
```

The three maps are HC A-beta-negative minus HC A-beta-positive, MCI
A-beta-positive, and AD A-beta-positive. The sign convention is retained in
the configuration, parcel tables, NIfTI filenames, and result metadata.

Create the manuscript-ready PDF and 300-dpi PNG summary figure with:

```bash
XDG_CACHE_HOME=/private/tmp/codex-neuromaps-cache \
  .venv-neuromaps/bin/python \
  neuromaps_analysis/neuromaps_code/\
plot_EI_HC_baseline_contrasts_lam001.py
```

Panels A--C show the surface-vertex density and fitted association for each
contrast; panel D shows the observed correlations against their spatial-null
distributions. The figure script verifies that each recomputed Pearson
correlation matches the saved analysis summary before plotting.

An alternative compact boxplot figure combines the HC amyloid-negative
group-mean result with the three signed HC baseline contrasts:

```bash
.venv-neuromaps/bin/python \
  neuromaps_analysis/neuromaps_code/plot_EI_turbulence_boxplots_lam001.py
```

Each box represents the distribution of Pearson correlations across 1,000
spatial rotations, with the observed correlation overlaid as a point. The
plot deliberately separates the HC group-mean map from the three contrast
maps because they answer related but distinct spatial questions.

## Confirmatory participant-level E:I coupling test

For a step-by-step explanation of the scientific calculations and their code
locations, see
`neuromaps_code/subject_level/PARTICIPANT_LEVEL_ANALYSIS.md`.

The confirmatory workflow reads the canonical, harmonized physical
lambda=0.01 N145 node table produced by
`turbulence/sch1000_N238rev/harmonization_allfeat/run_harmonization.py`, checks
it against the existing N145 machine-learning PTID set and retained HC-minus-AD
map, projects every participant to fsaverage 41k, and calculates
participant-level E:I--turbulence coupling. Run the harmonization stage first,
then run:

```bash
XDG_CACHE_HOME=/private/tmp/codex-neuromaps-cache \
  .venv-neuromaps/bin/python \
  neuromaps_analysis/neuromaps_code/subject_level/\
run_subjectlevel_ei_turbulence_coupling.py
```

To check the upstream workbook and cohort linkage without projecting surfaces
or rewriting analysis outputs, append `--validate-inputs-only`.

The primary outcome is the Fisher-transformed participant correlation. A
spatially normalized z score from 1,000 shared Alexander--Bloch rotations is a
sensitivity outcome. Both are tested for an ordered HC amyloid-negative,
HC amyloid-positive, MCI amyloid-positive, AD amyloid-positive trend using
10,000 Freedman--Lane permutations with age, gender, and education covariates.
The three planned HC amyloid-negative comparisons receive BH-FDR correction.
Participant-level outputs remain under `results/` and are ignored by Git.

Create the corresponding manuscript figure with the normative E:I surface map
and the spatially normalized participant-level result using:

```bash
XDG_CACHE_HOME=/private/tmp/codex-neuromaps-cache \
  .venv-neuromaps/bin/python \
  neuromaps_analysis/neuromaps_code/subject_level/\
plot_ei_surface_spin_coupling_lam001.py
```

Panel A renders the raw E:I expression ratio on lateral and medial fsaverage
41k surfaces. Panel B shows participant spin-z coupling values and the adjusted
ordered trend across the four amyloid-clinical stages. The script reads the
saved confirmatory tables and does not repeat the spatial rotations or
between-subject permutations.

Export the spatially normalized participant panel alone, without a panel
letter, using:

```bash
.venv-neuromaps/bin/python \
  neuromaps_analysis/neuromaps_code/subject_level/\
plot_spin_coupling_lam001.py
```

## Data provenance

- Turbulence difference maps were derived from the current N145 ComBat tables.
- AHBA expression maps must be regenerated or supplied locally.
