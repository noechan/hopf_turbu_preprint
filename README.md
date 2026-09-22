# Whole-brain turbulent dynamics across biological stages of Alzheimer's disease

This repository contains the analysis code accompanying the manuscript
*Turbulent dynamics and in silico perturbation responsiveness are disrupted across
biological stages of Alzheimer's disease*.

The publication workflow combines empirical multiscale turbulent-dynamics
measures, group- and subject-level Hopf models, covariate-adjusted statistical
analyses, machine-learning classification, and spatial correlations with
Neurosynth and Allen Human Brain Atlas maps.

## Repository structure

```text
turbulence/          Data preparation, empirical dynamics, ComBat, statistics,
                     and manuscript visualizations
hopf/group_level/    Group-level Hopf fitting and perturbation simulations
hopf/subject_level/  Subject-level Hopf fitting and perturbation simulations
hopf/gec/            Group generative effective-connectivity estimation
machine_learning/    Logistic-regression, mRMR, ROC, and SHAP workflows
abagen_analysis/     AHBA preprocessing and regional gene-expression maps
neuromaps_analysis/  Neurosynth and gene-expression spatial analyses
```

## Data availability and privacy

Participant-level ADNI data and derived participant-level tables are not
distributed. Researchers must obtain access through ADNI and construct the
expected local inputs described in [DATA.md](DATA.md). The repository contains
only source code and small group-level spatial maps used by the spatial
analyses.


## Software

- MATLAB R2026a with Signal Processing and Statistics and Machine Learning
  Toolboxes for the empirical and Hopf analyses.
- R with `readxl`, `permuco`, and the packages listed in the analysis READMEs.
- Python 3.12 for turbulence, harmonization, and neuromaps.
- Python 3.9 for the machine-learning environment used in the manuscript.

Install the Python environments separately:

```bash
python3.12 -m venv .venv-turbulence
.venv-turbulence/bin/pip install -r turbulence/requirements.txt

python3.9 -m venv .venv-ml
.venv-ml/bin/pip install -r machine_learning/requirements.txt

python3.12 -m venv .venv-neuromaps
.venv-neuromaps/bin/pip install -r neuromaps_analysis/requirements.txt
```

## Local paths

Set the external ADNI root before running analyses that consume restricted
inputs:

```bash
export ADNI3_ROOT=/your/local/path/to/ADNI3
```

The turbulence pipeline resolves repository paths from the script location.
Some original Hopf preprocessing scripts retain visible `/path/to/...`
placeholders because the input layout is institution-specific; replace these
with the corresponding authorized local paths. No personal filesystem paths
are retained in this publication copy.

For the AHBA analysis, set:

```bash
export AHBA_GENE_NIFTI_DIR=/your/local/path/to/schaefer1000_gene_niftis
```

The MATLAB cortical-rendering entry points resolve their bundled Schaefer and
Desikan--Killiany surface assets through `setup_sch1000_paths.m`; no
machine-specific rendering path is required.

## Manuscript workflow

The main execution order is:

1. Run `turbulence/sch1000_N238rev/prepare_data/` stages 1--5.
2. Run the four scripts in
   `turbulence/sch1000_N238rev/calculate_turbu/amyloid_status/`.
3. Export ComBat inputs with
   `turbulence/sch1000_N238rev/data_export/export_harmonization_inputs.m`.
4. Run `harmonization_allfeat/run_harmonization.py`, providing the authorized
   local post-ComBat exclusion manifest.
5. Run the education-adjusted R and Python analyses under
   `turbulence/sch1000_N238rev/statistical_analysis/`.
6. Run the group- and subject-level Hopf workflows under `hopf/`.
7. Run the two classification entry points documented in
   `machine_learning/README.md`.
8. Run the spatial analyses documented in
   `neuromaps_analysis/README.md`.
