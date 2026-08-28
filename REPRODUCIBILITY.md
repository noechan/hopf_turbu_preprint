# Manuscript analysis specifications

## Cohort and harmonization

- Raw harmonization sample: N=152.
- Six recursively identified single-site MRI participants are removed before
  fitting ComBat, leaving N=146.
- ComBat preserves age, sex, education, and categorical group while treating
  MRI site as the batch variable.
- One authorized post-ComBat multimodal-cohort exclusion produces N=145.

## Group-level empirical dynamics

- Five planned stage contrasts are tested at each spatial scale.
- Freedman--Lane residual permutations: 100,000.
- Covariates: age, sex, and years of education.
- Pairwise BH-FDR family: the five planned contrasts within each scale and
  measure.
- Information transfer is reported as `1 - InformationTransfer`.

## Node-level metastability

- Scale: `lambda=0.01`.
- Schaefer atlas: 1,000 parcels.
- Freedman--Lane permutations: 10,000.
- Covariates: age, sex, and education.
- BH-FDR across 1,000 nodes separately within each planned contrast.

## MOCA associations

- Predictors: turbulence, information flow, and
  `1 - InformationTransfer`, each at `lambda=0.01`.
- Covariates: age, sex, and education.
- The predictor is standardized; MOCA remains in raw points.
- BH-FDR across the three prespecified predictors.
- LMG relative importance with 10,000 bootstrap resamples is descriptive.

## AT(N) associations

- HC3 robust ordinary least squares.
- Covariates: age, sex, and education.
- Outcomes: amyloid burden, bilateral hippocampal volume, and three regional
  tau summaries.
- BH-FDR across the five outcomes separately for each dynamical predictor.

## Machine learning

- Four pairwise classifications: HC Aβ− vs HC Aβ+, HC Aβ− vs MCI Aβ+,
  HC Aβ− vs AD Aβ+, and MCI Aβ+ vs AD Aβ+.
- Two feature sets: empirical dynamics alone and empirical dynamics plus
  subject-level information capability and susceptibility.
- Forty stratified random 75/25 train/test partitions.
- mRMR is fitted on training data only.
- Inner leave-one-out cross-validation selects 5--15 features, regularization
  strength, and L1/L2 penalty.
- Class-balanced `liblinear` logistic regression.
- The Youden threshold is estimated from training predictions.
- SHAP is evaluated on held-out observations.

## Spatial analyses

- Turbulence maps are absolute N145 group-mean differences.
- Neurosynth: three memory terms by three scales, Pearson correlation,
  fsaverage-41k projection, 1,000 Alexander--Bloch rotations, seed 1234, and
  BH-FDR across nine tests within contrast.
- Gene analysis: 15 prespecified genes, `lambda=0.01` manuscript contrasts,
  fsaverage-41k projection, Pearson correlation, and 1,000 Alexander--Bloch
  rotations with seed 1234.
