# Statistical analysis

Statistics must consume either the recomputed ComBat tables in
`harmonization_allfeat/recomputed/` or the self-describing standalone
amyloid-status result files. Diagnosis-wide results must not be split with a
separate current PTID list because historical MCI column order is different.

## Demographic and clinical table

`demographic_analysis/run_demographic_tests_N145.R` recomputes the complete
N145 demographic and clinical table. It uses Kruskal--Wallis tests for
continuous variables, Pearson's chi-square test for sex, and Holm-corrected
Dunn post-hoc tests after significant continuous omnibus tests. Education is
joined by PTID from the local copy of the same harmonisation metadata used by
ComBat and the education-adjusted permutation analyses. See
`demographic_analysis/README.md` for commands, provenance, and output details.

## Active harmonized analysis

`turbu_harm_stats/run_group_permutations_N145_age_sex_education.R` is the
manuscript turbulence group analysis. It reads
`harmonization_allfeat/recomputed/Turbu_ComBat_ADNI3_allfeatures_N145.xlsx`
directly and joins age, sex, and education from the authorized metadata at
`data/covariates/`. It validates the four N145 group counts before running any
model.

For each of the ten turbulence scales, it fits:

```text
Turbu_lam_* ~ Group + Age + Sex + Education
```

The metadata coding is converted to a factor as `0 = Female` and `1 = Male`.

The omnibus group test uses the marginal factor test from `permuco::aovperm`
with sum coding; the five manuscript pairwise contrasts use
`permuco::lmperm`. Both use Freedman-Lane permutations. The
pairwise p-values are BH-FDR corrected across the five planned stage contrasts
separately within each spatial scale, matching the manuscript analysis.
Scanner is not included after ComBat harmonization. Only aggregate result
tables may be published.

## Node-level metastability

The manuscript node-wise analysis is
`nodewise_metastability/run_nodewise_metastability_N145_age_sex_education.R`.
It uses a two-sided vectorized Freedman--Lane test of `Group` while adjusting
for age, sex, and education, followed by BH-FDR across the 1,000 parcels
separately within each Figure 4 contrast at `lambda=0.01`. The age/sex-only
entry point is retained as a reduced-covariate sensitivity analysis.

`nodewise_metastability/run_yeo7_top30_N145.py` regenerates the Figure 4 Yeo-7
network summaries directly from the age-, sex-, and education-adjusted CSVs.
It preserves the manuscript's lower-30th-percentile rule with threshold ties
and reports both raw parcel counts and atlas-size-normalised representation
ratios.

## AT(N) biomarker regressions

`atn_biomarkers/run_atn_regressions_N145.py` is the statistical source for
Figure 4 panels c, d, j, k, and l. It fits both the age/sex and the
age/sex/education specifications with HC3 robust ordinary least squares and
applies BH-FDR across the five turbulence-coefficient tests separately within
each model. Education is joined by PTID from the canonical ComBat metadata
column `edu`. Only aggregate results and prediction lines are saved; the
participant-level merged table remains in memory. The corresponding Python
visualizer reads these saved outputs and no longer refits the models.

## Subject-level perturbation measures

`perturbation_measures/run_subjectlevel_perturbation_group_N145.R` is the
active group analysis for participant-level information capacity and
susceptibility. It fits education-adjusted Freedman--Lane models and applies
BH-FDR across the five planned stage comparisons separately within each
measure. `perturbation_measures/run_moca_subjectlevel_perturbation_N145.py`
tests the cross-sectional association of each measure with MOCA using HC3
multiple linear regression. It reports both an age-, sex-, and
education-adjusted total association and a separate group-adjusted sensitivity
estimand. `perturbation_measures/run_atn_subjectlevel_perturbation_N145.py`
tests each measure separately against the five prespecified AT(N) outcomes,
with age, sex, and education as nuisance covariates and BH-FDR across outcomes
within measure. See `perturbation_measures/README.md` for inputs, correction
families, and run commands.

Install the two R dependencies once:

```r
install.packages(c("readxl", "permuco"))
```

From the publication repository's `turbulence/` directory, validate the inputs:

```bash
Rscript sch1000_N238rev/statistical_analysis/turbu_harm_stats/run_group_permutations_N145_age_sex_education.R \
  --validate-only
```

Then run the canonical 100,000-permutation analysis:

```bash
Rscript sch1000_N238rev/statistical_analysis/turbu_harm_stats/run_group_permutations_N145_age_sex_education.R
```

For a quick technical smoke test before the full run, use `--n-perm 100`.
Results are written under
`statistical_analysis/turbu_harm_stats/results/N145_group_permutations/`. Set
`ADNI3_ROOT` or pass `--metadata-file` if the external drive is mounted
elsewhere.

An alternative post-processing script,
`turbu_harm_stats/recalculate_pairwise_FDR_by_lambda_N145.R`, applies BH-FDR
across the five planned AD-stage contrasts separately within each lambda. It
uses the canonical pairwise permutation table without refitting models and
writes a distinct `_FDR_by_lambda.csv` result so the two correction families
remain auditable.

The active cross-sectional MOCA analysis is implemented in
`cognition_moca/run_moca_dynamics_N145.R`. It joins MOCA and the canonical
harmonization demographics by PTID while taking every dynamical value directly
from the recomputed N145 ComBat table. Its manuscript default tests three
prespecified global measures with two-sided Freedman--Lane inference in both a
demographic-adjusted total-association model and a group-adjusted model. See
`cognition_moca/README.md` for the estimands, correction families, and run
commands.

The post-ComBat exclusion is retained to match the common multimodal N145
manuscript cohort; its participant identifier remains in the authorized local
manifest and is not distributed.
