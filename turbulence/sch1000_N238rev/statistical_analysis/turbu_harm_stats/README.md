# N145 harmonized turbulence statistics

The manuscript group analysis is
`run_group_permutations_N145_age_sex_education.R`. It uses the N145 cohort,
five planned contrasts, Freedman--Lane permutations, and the model:

```text
Outcome ~ Group + Age + Sex + Education
```

## RStudio

Open `turbu_harm_stats.Rproj`, open the R script, and click **Source**. Do not
select the entire script and click **Run**, because sending a large selection
line-by-line can leave the Console at the continuation (`+`) prompt.

Correct execution begins with messages like:

```text
Validated canonical N145 cohort: HC_ABneg=51 | HC_ABpos=37 | MCI_ABpos=31 | AD_ABpos=26
Validated 10 turbulence outcomes with complete age, sex, and education covariates.
Omnibus: Turbu_lam_0_27
```

## Terminal

From this project directory, the equivalent command is:

```bash
Rscript run_group_permutations_N145_age_sex_education.R
```

Use `--validate-only` to check the inputs without fitting models, or
`--n-perm 100` for a short technical smoke test. The canonical analysis uses
the default 100,000 permutations. Aggregate outputs are written to
`results/N145_group_permutations/` inside this project.

Manuscript outputs use the suffix `age_sex_education`. The age/sex-only script
is retained as a sensitivity analysis. The manuscript run summary is written as
`N145_group_permutation_run_summary_age_sex_education.csv`.

## Alternative FDR family: contrasts within each lambda

`recalculate_pairwise_FDR_by_lambda_N145.R` preserves the fitted
Freedman-Lane models and recalculates only the pairwise BH-FDR family. It
corrects across the five planned AD-stage contrasts separately at each of the
ten lambda values. The original comparison-wise result remains unchanged.

Run it after the canonical permutation analysis:

```bash
Rscript recalculate_pairwise_FDR_by_lambda_N145.R
```

The alternative table is written as
`results/N145_group_permutations/N145_turbulence_pairwise_FreedmanLane_age_sex_FDR_by_lambda.csv`.
The output contains only the lambda-wise FDR columns, preventing downstream
scripts from accidentally mixing the two multiplicity definitions. This
alternative family should be used only when the scientific question and
reporting plan define the five stage contrasts within a lambda as the family.

To apply the same lambda-wise FDR family to the education-adjusted secondary
analysis, run:

```bash
Rscript recalculate_pairwise_FDR_by_lambda_N145.R \
  --input-file results/N145_group_permutations/N145_turbulence_pairwise_FreedmanLane_age_sex_education.csv \
  --output-file results/N145_group_permutations/N145_turbulence_pairwise_FreedmanLane_age_sex_education_FDR_by_lambda.csv
```

## Information flow, cascade, and one-minus information transfer

`run_information_measures_permutations_N145.R` applies the same canonical N145
cohort, five planned stage contrasts, Freedman--Lane method, permutation count,
random seed, and age/sex adjustment to the other manuscript measures. It tests:

- nine information-flow scales (lambda 0.24 through 0.01);
- ten `1 - InformationTransfer` scales (lambda 0.27 through 0.01); and
- the scalar information-cascade measure.

The `1 - InformationTransfer` transformation is applied to the harmonized N145
values before model fitting. Consequently, its two-sided permutation p-values
are equivalent to the raw-value analysis, while the signs of the t statistic
and adjusted Hedges' g are reversed and directly describe `1 - transfer`.

For pairwise tests, BH-FDR is calculated across the five planned stage
contrasts separately within each scale. Information flow, one-minus information
transfer, and information cascade remain separate multiplicity families.

Run the age/sex sensitivity analysis with:

```bash
Rscript run_information_measures_permutations_N145.R
```

Run the manuscript education-adjusted analysis with:

```bash
Rscript run_information_measures_permutations_N145_age_sex_education.R
```

Use `--validate-only` to validate the 9/10/1 outcome structure without fitting,
or `--n-perm 100` for a short technical smoke test. The canonical analyses use
the default 100,000 permutations. Measure-specific omnibus and pairwise CSVs,
plus a run summary documenting the transformation and FDR families, are written
to `results/N145_group_permutations/`.
