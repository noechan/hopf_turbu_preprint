# Figure 4 AT(N) regression provenance

Figure 4 panels c, d, j, k, and l are now generated in two explicit stages:

1. `statistical_analysis/atn_biomarkers/run_atn_regressions_N145.py` fits and
   saves the statistical models.
2. `visualization/python/scripts/plot_atn_regressions.py` reads the aggregate
   results and saved prediction lines and creates the figures without refitting.

## Inputs

- Turbulence and group membership come from the current N145 all-feature
  ComBat workbook.
- Age, sex, amyloid Centiloid (`CL_pvc`), and three regional PVC tau measures
  come from the mounted ADNI biomarker workbook.
- Left and right hippocampal GMV come from the mounted post-ComBat VBM table;
  bilateral GMV is their row-wise mean.
- Education comes from the `edu` column of the controlled local copy at
  `data/covariates/covariates_ADNI3_ABeta_N152.csv` and is joined by unique
  PTID.

The workflow validates HC Aβ−=51, HC Aβ+=37, MCI Aβ+=31, and AD Aβ+=26.
Amyloid and hippocampal analyses use N=145; tau analyses use the 134
participants with regional tau data. Participant-level merged data are not
written into the repository.

## Statistical models

Both specifications are retained:

```text
outcome ~ Turbu_lam_0_01 + AGE + SEX_NUM
outcome ~ Turbu_lam_0_01 + AGE + SEX_NUM + EDUCATION
```

Sex is coded female=0 and male=1. Ordinary least squares uses HC3 robust
standard errors. BH-FDR is applied jointly to the five turbulence-coefficient
p-values separately within each covariate specification. Group determines
point colour only and is not an additional model covariate.

The displayed beta is the unstandardized turbulence coefficient, and the
displayed R² is the full-model R². Prediction lines fix age at the model-sample
mean and sex at its mode; the education-adjusted line additionally fixes
education at its model-sample mean.

All five AT(N) associations survive FDR correction in both specifications.
Adding education slightly attenuates the absolute turbulence coefficients but
does not change their direction or inferential conclusion.

## Outputs and run command

Aggregate statistics are stored under:

```text
statistical_analysis/atn_biomarkers/results/N145_atn_regressions/
```

Figures are separated by model under:

```text
figures_N145/sch1000/Abeta_Status/harmonized_allfeat/python/atn_regressions/
├── age_sex/
└── age_sex_education/
```

From the repository root, rerun statistics and both figure sets with:

```bash
.venv/bin/python \
  sch1000_N238rev/visualization/python/run_harmonized_visualizations.py \
  --stage atn
```
