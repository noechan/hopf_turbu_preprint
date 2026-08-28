# N145 cross-sectional MOCA--brain-dynamics associations

This folder contains the active cross-sectional analysis of whether harmonized
brain-dynamical measures are associated with MOCA in the canonical N145
amyloid-status cohort.

## Estimands

The script fits two separate models for every selected dynamical predictor:

```text
Total association:
MOCA ~ z(dynamics) + Age + Sex + Education

Association beyond disease stage:
MOCA ~ z(dynamics) + Age + Sex + Education + Group
```

The first model estimates the demographic-adjusted association across the
biological-stage continuum. The second asks whether a dynamical measure is
associated with MOCA after accounting for the four-level amyloid-status group.
It should not be interpreted as the same estimand with a merely optional
covariate: group may lie on the disease pathway.

The dynamical predictor is standardized within the analysis sample, so the
reported coefficient is the expected difference in raw MOCA points per
one-standard-deviation increase in the predictor. Age and education remain in
their original units; sex and group are categorical.

## Inference and multiplicity

The predictor is tested using a two-sided Freedman--Lane residual-permutation
test. The default is 100,000 permutations. Parametric 95% confidence intervals,
partial R-squared, full-model R-squared, and the incremental R-squared over the
covariate-only model are saved as descriptive effect-size information.

The manuscript predictor set contains three prespecified measures:

- turbulence at `lambda=0.01`;
- information flow at `lambda=0.01`;
- `1 - information transfer` at `lambda=0.01`.

BH-FDR is applied across these three predictors separately within each model.
The output also contains a within-measure-family FDR column. With the primary
set, the latter is descriptive because each family contains one predictor.

The optional `all-scales` set includes all 32 harmonized dynamical features.
For that exploratory analysis, the global FDR column covers all 32 predictors,
and the within-family column corrects across scales separately for turbulence,
information flow, and `1 - information transfer`.

## Input provenance

- Current dynamical values and canonical group membership are read from
  `harmonization_allfeat/recomputed/Turbu_ComBat_ADNI3_allfeatures_N145.xlsx`.
- MOCA is joined by PTID from
  `visualization/python/data/Turbu_ComBat_clin_ADNI3_allfeatures_N145.csv`.
  Only `PTID`, `Group`, and `MOCA` are used from this file. Its stale dynamical
  columns are explicitly ignored.
- Age, sex, and education are joined by PTID from the controlled local copy at
  `data/covariates/covariates_ADNI3_ABeta_N152.csv`.

The script validates PTID uniqueness, group agreement across all three inputs,
the canonical group counts (51/37/31/26), numeric covariates, and the selected
dynamical columns before fitting any model.

## Run

From the repository root, validate the inputs without writing files:

```bash
Rscript sch1000_N238rev/statistical_analysis/cognition_moca/run_moca_dynamics_N145.R \
  --validate-only
```

Run the confirmatory analysis:

```bash
Rscript sch1000_N238rev/statistical_analysis/cognition_moca/run_moca_dynamics_N145.R
```

For a short technical smoke test, add `--n-perm 100`. To run the exploratory
all-scale analysis:

```bash
Rscript sch1000_N238rev/statistical_analysis/cognition_moca/run_moca_dynamics_N145.R \
  --predictor-set all-scales
```

The script also accepts `--model-set total` or
`--model-set group-adjusted`. Results are written under
`results/N145_MOCA_associations/`; primary and exploratory filenames remain
distinct.

## Current confirmatory results

The 100,000-permutation analysis uses 144 participants. MOCA is missing for
one HC Aβ+ participant, leaving group counts of 51/36/31/26.

| Predictor | Total-association beta | Total pFDR | Group-adjusted beta | Group-adjusted pFDR |
|---|---:|---:|---:|---:|
| Turbulence, `lambda=0.01` | 1.146 | 0.00497 | -0.031 | 0.91090 |
| Information flow, `lambda=0.01` | 1.522 | 0.00051 | 0.155 | 0.91090 |
| `1 -` information transfer, `lambda=0.01` | 1.283 | 0.00308 | 0.200 | 0.91090 |
All three prespecified measures are positively associated with MOCA in the
model adjusted for age, sex, and education, and all survive FDR correction.
None remains
significant after additionally adjusting for amyloid-status group. Thus, the
cross-sectional association is consistent with variation across disease stages
and does not currently provide evidence for an association within/beyond stage.

## Focused hierarchical information-flow model

`run_moca_hierarchical_information_flow_N145.R` implements the prespecified
blockwise model that asks whether information flow at `lambda=0.01` explains
MOCA beyond demographics, amyloid-status stage, and turbulence:

```text
M0: MOCA ~ Age + Sex + Education + Group
M1: M0 + z(Turbulence at lambda=0.01)
M2: M1 + z(Information flow at lambda=0.01)
```

The primary comparison is `M1 -> M2`. The script reports change in R-squared,
the conventional nested-model F test, a 100,000-draw two-sided Freedman--Lane
test for the added predictor, and HC3 heteroscedasticity-robust coefficient
inference. It also saves model coefficients, assumption diagnostics, and
flagged influence observations.

Run it from the repository root:

```bash
Rscript sch1000_N238rev/statistical_analysis/cognition_moca/run_moca_hierarchical_information_flow_N145.R
```

Results are written to
`results/N145_MOCA_hierarchical_information_flow/`.

The current 100,000-permutation run used 144 participants (51 HC Aβ-, 36 HC
Aβ+, 31 MCI Aβ+, and 26 AD Aβ+). The demographic-and-stage model explained
65.52% of MOCA variance. Adding turbulence did not improve the model
(`delta R2 = 0.00003`, Freedman--Lane `p = 0.91081`). Adding information flow
after turbulence increased explained variance by 0.28 percentage points
(`delta R2 = 0.00275`). Its coefficient was 0.522 MOCA points per predictor SD,
but the increment was not significant (Freedman--Lane `p = 0.30124`; HC3
`p = 0.41025`, HC3 95% CI `[-0.728, 1.772]`). Because the final information-flow
block contains one prespecified test, no FDR correction is required for that
primary comparison.

## AT(N)-style information-flow regression

`run_moca_information_flow_regression_N145.py` mirrors the active AT(N)
regression design and fits:

```text
MOCA ~ z(InformationFlow at lambda=0.01) + Age + Sex + Education
```

It uses HC3 heteroscedasticity-consistent inference. Diagnostic group is not a
regression covariate; it is retained only to colour the corresponding scatter
plot. The script reports VIF values for information flow and every included
covariate, along with residual linearity, normality, and heteroscedasticity
diagnostics. Each model contains one outcome and one dynamical predictor. The
raw HC3 result is therefore the within-model test; BH-FDR is additionally
applied across the three parallel dynamical models in the combined comparison.

From the repository root, run the analysis and visualization together with:

```bash
.venv/bin/python \
  sch1000_N238rev/visualization/python/run_harmonized_visualizations.py \
  --stage cognition
```

Aggregate results are saved under
`results/N145_MOCA_information_flow_regression/`; the PDF and PNG are saved
under the active `figures_N145/.../python/cognition_moca/` directory.

The current HC3 analysis uses 144 participants and finds a positive association
between information flow and MOCA (`beta = 1.522` MOCA points per predictor SD,
HC3 95% CI `[0.661, 2.382]`, `p = 0.000529`). Information flow adds 7.97
percentage points to model R-squared beyond age, sex, and education; the full
model has `R2 = 0.2271`. Its VIF is 1.19, and all non-intercept VIFs are at or
below 1.19, so there is no evidence of problematic multicollinearity in this
specification. Residual heteroscedasticity remains detectable, motivating HC3
inference; the corresponding 100,000-draw Freedman--Lane analysis independently
gives `p = 0.00017`.

The same entry points also fit turbulence and `1 - information transfer` at
`lambda=0.01`. Use `--predictor turbulence` or
`--predictor one-minus-information-transfer` with the statistical or plotting
script. The unified `--stage cognition` command runs and plots all three
prespecified single-predictor models and keeps their outputs separate.
In the current turbulence model, `beta = 1.146` MOCA points per predictor SD
(HC3 95% CI `[0.416, 1.876]`, `p = 0.00209`). Turbulence adds 4.90 percentage
points to R-squared beyond age, sex, and education, and the full model has
`R2 = 0.1963`. The turbulence VIF is 1.09, so multicollinearity is negligible.
The `1 - information transfer` model gives `beta = 1.283` MOCA points per
predictor SD (HC3 95% CI `[0.327, 2.240]`, `p = 0.00854`), adds 5.98 percentage
points to R-squared, and has `R2 = 0.2071`. Its VIF is 1.12.

Across the three-model family, BH-FDR gives `pFDR = 0.00159` for information
flow, `pFDR = 0.00313` for turbulence, and `pFDR = 0.00854` for
`1 - information transfer`; all three remain significant. The combined
comparison is written to
`results/N145_MOCA_single_predictor_regression_comparison.csv`.

## LMG relative importance

`run_moca_lmg_relative_importance_N145.py` formally compares how much MOCA
variance is attributable to turbulence, information flow, and
`1 - information transfer` at `lambda=0.01`. The default model adjusts for age,
sex, and education, matching the AT(N)-style analyses. LMG importance averages
each measure's incremental R-squared over all six possible entry orders. A
10,000-draw paired participant bootstrap provides percentile confidence
intervals and tests whether information flow contributes more than each of the
other two measures; BH-FDR covers those two comparisons.

Run from the repository root:

```bash
.venv/bin/python \
  sch1000_N238rev/statistical_analysis/cognition_moca/run_moca_lmg_relative_importance_N145.py
```

Add `--include-group` for the separate beyond-stage sensitivity estimand.
Aggregate importance, pairwise tests, all subset-model R-squared values, the
bootstrap distribution, and provenance are saved under
`results/N145_MOCA_LMG_relative_importance/`.

In the current N144 analysis, information flow has the largest LMG point
contribution (`LMG R2 = 0.0410`, 50.7% of the joint dynamics contribution),
followed by `1 - information transfer` (`0.0227`, 28.1%) and turbulence
(`0.0171`, 21.2%). The information-flow-minus-turbulence difference is `0.0239`
(bootstrap 95% CI `[-0.0141, 0.0667]`), and the information-flow-minus-
`1 - information transfer` difference is `0.0183` (95% CI
`[-0.0156, 0.0537]`). Neither comparison is significant after BH-FDR
(`pFDR = 0.2656` for both). Therefore, information flow ranks highest
descriptively, but the present sample does not establish that its relative
contribution is statistically greater than either alternative.
