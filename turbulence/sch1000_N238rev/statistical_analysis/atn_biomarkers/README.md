# N145 turbulence--AT(N) biomarker regressions

This folder contains the statistical analysis underlying manuscript Figure 4
panels c, d, j, k, and l. Statistics are computed here and are no longer fitted
inside the visualization script.

## Models

The same five outcomes are fitted independently with HC3 robust ordinary least
squares:

```text
Age/sex model:
outcome ~ Turbu_lam_0_01 + AGE + SEX_NUM

Education-adjusted model:
outcome ~ Turbu_lam_0_01 + AGE + SEX_NUM + EDUCATION
```

The outcomes are amyloid Centiloid, bilateral hippocampal GMV, and mesial,
metatemporal, and temporoparietal tau. BH-FDR is applied to the five turbulence
coefficient p-values separately within each model. Sex is coded female=0 and
male=1. Group determines scatter-point colour only and is not a regression
covariate.

## Inputs

- Turbulence and canonical N145 group membership:
  `harmonization_allfeat/recomputed/Turbu_ComBat_ADNI3_allfeatures_N145.xlsx`.
- Amyloid, age, sex, and regional tau: mounted
  `ADNI3_N238rev_with_ABETA_Status_CL24_tau_regional.xlsx`.
- Bilateral hippocampal GMV: mounted post-ComBat `ADNI3_VBM_postCOMBAT.csv`.
- Education: `edu` joined by PTID from the controlled local copy at
  `data/covariates/covariates_ADNI3_ABeta_N152.csv`.

Amyloid and hippocampal models use N=145; regional tau models use N=134. The
script validates these samples and the four canonical group counts. It does
not save participant-level merged data.

## Run commands

From the repository root, run both models:

```bash
.venv/bin/python \
  sch1000_N238rev/statistical_analysis/atn_biomarkers/run_atn_regressions_N145.py
```

The explicit single-model entry points are:

```bash
.venv/bin/python \
  sch1000_N238rev/statistical_analysis/atn_biomarkers/run_atn_regressions_N145_age_sex.py

.venv/bin/python \
  sch1000_N238rev/statistical_analysis/atn_biomarkers/run_atn_regressions_N145_age_sex_education.py
```

Add `--validate-only` to validate without fitting or writing outputs.

Aggregate results, saved prediction lines, a direct model-comparison table,
and JSON provenance are written under `results/N145_atn_regressions/`.

## Current results

| Outcome | Age/sex beta | Age/sex pFDR | Age/sex/education beta | Age/sex/education pFDR |
|---|---:|---:|---:|---:|
| Amyloid Centiloid | -690.2256 | 0.02681 | -649.7013 | 0.03917 |
| Bilateral hippocampal GMV | 0.6640 | 0.01685 | 0.6118 | 0.02799 |
| Mesial tau | -6.6544 | 0.01420 | -6.5582 | 0.01601 |
| Metatemporal tau | -8.1308 | 0.01685 | -7.9704 | 0.02212 |
| Temporoparietal tau | -8.0828 | 0.02681 | -7.9844 | 0.03103 |

All five associations survive FDR correction in both specifications. Adding
education attenuates the absolute turbulence coefficient by approximately
1.2--7.9% without changing the direction or significance conclusion.

## Stage-specific sensitivity analysis

`run_atn_turbulence_stage_sensitivity_N145.py` tests whether the pooled
associations are also present within disease stages. It fits the
age/sex/education-adjusted model separately in HC Aβ-, HC Aβ+, MCI Aβ+, and
AD Aβ+, a group-adjusted common-slope model, and a formal
turbulence-by-group interaction model. Turbulence is standardized over each
outcome's complete analysis sample so that slopes are expressed per SD.

Run from the repository root (with `ADNI3_ROOT` pointing to the mounted ADNI3
directory):

```bash
.venv-neuromaps/bin/python \
  turbulence/sch1000_N238rev/statistical_analysis/atn_biomarkers/run_atn_turbulence_stage_sensitivity_N145.py
```

None of the 20 within-stage slopes survives BH-FDR, and none of the five
group-adjusted common slopes is significant after BH-FDR. The five joint
interaction tests are also non-significant. Aggregate CSV files and a forest
plot are written under
`results/N145_ATN_turbulence_stage_sensitivity/`.

## Figures

The unified visualization entry point reruns both models and creates separate
figure directories for the two specifications:

```bash
.venv/bin/python \
  sch1000_N238rev/visualization/python/run_harmonized_visualizations.py \
  --stage atn
```

The visualizer reads the saved aggregate results and prediction lines; it does
not refit the regressions.
