# Manuscript node-level specification

Figure 4 uses node-level metastability, defined as the temporal standard
deviation of the local Kuramoto order parameter, at physical
`lambda=0.01`. The analysis contains 1,000 Schaefer parcels and two planned
contrasts: HC A$\beta^-$ vs. AD A$\beta^+$ and MCI A$\beta^+$ vs.
AD A$\beta^+$.

The publication inference is implemented in
`statistical_analysis/nodewise_metastability/run_nodewise_metastability_N145_age_sex_education.R`.
It uses two-sided Freedman--Lane residual permutations (10,000), adjusts for
age, sex, and years of education, and applies Benjamini--Hochberg FDR across
the 1,000 parcels separately within each contrast. The age/sex-only script is
retained as a sensitivity analysis; the historical unadjusted MATLAB analysis
is not part of this publication copy.

Surface maps display the complete $-\log_{10}(p_{\mathrm{FDR}})$ vector. The
Yeo-7 summary retains nodes at or below the 30th percentile of the adjusted
permutation p-value distribution, including ties, and reports both raw network
counts and atlas-size-normalised representation ratios.

Historical suffixes must not be interpreted as array indices: `lam1`, `lam3`,
and `lam6` denote physical lambda values 0.01, 0.03, and 0.06. Active node-wise
filenames use the explicit physical value `lambda_0_01`.
