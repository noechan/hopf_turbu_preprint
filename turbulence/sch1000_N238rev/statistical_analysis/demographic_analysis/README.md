# N145 demographic and clinical analysis

`run_demographic_tests_N145.R` reproduces the demographic and clinical table
for the canonical HC Aβ−, HC Aβ+, MCI Aβ+, and AD Aβ+ N145 cohort.

The active clinical export supplies age, sex, cognition, and Centiloid values.
Years of education are joined by PTID from the same metadata CSV used for
ComBat and the education-adjusted permutation analyses. Where the clinical
export and harmonisation metadata disagree, the canonical harmonisation value
is used and the discrepancy is recorded in the private run summary.

The active script implements the tests described in the manuscript table:

- Kruskal--Wallis tests for continuous variables;
- Pearson's chi-square test for sex;
- Dunn pairwise tests with Holm correction after a significant continuous
  omnibus test;
- epsilon squared for significant Kruskal--Wallis tests.

From the publication repository's `turbulence/` directory, validate inputs with:

```bash
Rscript sch1000_N238rev/statistical_analysis/demographic_analysis/run_demographic_tests_N145.R \
  --validate-only
```

Run the analysis with:

```bash
Rscript sch1000_N238rev/statistical_analysis/demographic_analysis/run_demographic_tests_N145.R
```

Outputs are written to `results/N145_demographics/`, including descriptives,
omnibus tests, Dunn--Holm post-hoc tests, sex counts, provenance, and a complete
LaTeX table fragment that can be included in the manuscript.
