# N145 node-level metastability statistics

This folder contains the covariate-adjusted inferential analysis for the two
node-level contrasts shown in manuscript Figure 4:

1. HC Aβ− versus AD Aβ+;
2. MCI Aβ+ versus AD Aβ+.

The measure is the temporal standard deviation of the local Kuramoto order
parameter at `lambda=0.01`, stored for 1,000 Schaefer parcels in the current
N145 node-level ComBat table.

## Statistical model

The manuscript analysis fits, independently at each parcel:

```text
NodeMetastability ~ Group + Age + Sex + Education
```

Inference uses a two-sided Freedman--Lane residual-permutation test for the
group coefficient. A common participant-row permutation is applied to all
1,000 parcels at each iteration, preserving cross-parcel covariance. Raw
permutation p-values use the plus-one correction. Benjamini--Hochberg FDR is
then applied across the 1,000 parcels separately within each planned contrast.

The default is 10,000 permutations. A reduced-covariate sensitivity analysis
omits education:

```text
NodeMetastability ~ Group + Age + Sex
```

## Run commands

From the publication repository's `turbulence/` directory, validate the
reduced-covariate sensitivity inputs:

```bash
Rscript sch1000_N238rev/statistical_analysis/nodewise_metastability/run_nodewise_metastability_N145.R \
  --validate-only
```

Run that sensitivity analysis:

```bash
Rscript sch1000_N238rev/statistical_analysis/nodewise_metastability/run_nodewise_metastability_N145.R
```

Run the age-, sex-, and education-adjusted manuscript analysis:

```bash
Rscript sch1000_N238rev/statistical_analysis/nodewise_metastability/run_nodewise_metastability_N145_age_sex_education.R
```

Regenerate the manuscript Yeo-7 summary from the age-, sex-, and
education-adjusted results:

```bash
.venv-turbulence/bin/python \
  sch1000_N238rev/statistical_analysis/nodewise_metastability/run_yeo7_top30_N145.py
```

The Yeo-7 script applies the manuscript selection rule separately within each
contrast: parcels at or below the 30th percentile of `P_FDR_BH` are retained,
including ties at the threshold. It saves parcel assignments, raw network
counts, atlas-size-normalised representation ratios, and PDF/PNG figures under
`figures_N145/sch1000/Abeta_Status/harmonized_allfeat/RSN_top30_FDR_age_sex_education/`.
The bundled `RSN7vector.mat` is checked against the expected Schaefer-1000
Yeo-7 network sizes. By default, the script also verifies the frozen Figure 3
selected totals and network counts before writing plots. Use
`--skip-manuscript-check` only for an explicitly exploratory result set.

For a short technical test, add `--n-perm 100`. For higher p-value precision,
the scripts also accept `--n-perm 100000`, although this is substantially more
computationally expensive.

Outputs are written under `results/N145_nodewise_metastability/`. Each result
contains raw group means, the adjusted group coefficient, t statistic, partial
eta squared, permutation p-value, BH-adjusted p-value, and significance flag.

## Current 10,000-permutation results

| Contrast | Age + sex | Age + sex + education | Direction |
|---|---:|---:|---|
| HC Aβ− vs AD Aβ+ | 1,000 / 1,000 FDR-significant parcels | 1,000 / 1,000 | HC Aβ− > AD Aβ+ at every parcel |
| MCI Aβ+ vs AD Aβ+ | 974 / 1,000 | 917 / 1,000 | MCI Aβ+ > AD Aβ+ at every parcel |

The education-adjusted MCI-versus-AD result is a strict subset of the age/sex
result: 57 parcels lose FDR significance and no new parcel becomes significant.
The manuscript Figure 4 workflow uses the age-, sex-, and education-adjusted
CSVs.

## Current Yeo-7 top-30% results

Because BH-adjusted p-values are tied, the threshold retains 320 parcels for
HC Aβ− versus AD Aβ+ and 325 parcels for MCI Aβ+ versus AD Aβ+. DMN has the
largest raw parcel count in both selected sets (80 and 92 parcels), while LIM
has the largest atlas-size-normalised representation ratio (2.76 and 1.85).
VIS is under-represented in HC Aβ− versus AD Aβ+ (38 parcels; ratio 0.73) and
absent from the MCI Aβ+ versus AD Aβ+ selected set.

## Suggested Methods wording

> Node-level metastability at lambda = 0.01 was analysed separately at each of
> the 1,000 Schaefer parcels using two-sided Freedman--Lane residual
> permutation tests (10,000 permutations), with group as the predictor of
> interest and age, sex, and years of education as nuisance covariates. Benjamini--Hochberg FDR
> correction was applied across the 1,000 parcels separately within each
> planned contrast. A reduced-covariate sensitivity analysis omitted years of
> education.
