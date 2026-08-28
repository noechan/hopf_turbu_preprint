# Harmonization

This stage applies ComBat to raw turbulence features after the four standalone
amyloid-status calculations have completed.

## Inputs

Run `data_export/export_harmonization_inputs.m` in MATLAB first. It reads the
self-describing HC Aβ−, HC Aβ+, MCI Aβ+, and AD Aβ+ result files and writes
the all-feature and lambda=0.01, 0.03 and 0.06 node N152 workbooks to:

```text
$ADNI3_ROOT/timeseries/harmonization_inputs/sch1000_N238rev
```

The tables contain PTID and Group explicitly. They do not reconstruct groups
from diagnosis-wide output ordering. The node tables are named by physical
scale. Lambda 0.01, 0.03 and 0.06 are extracted from MATLAB array indices 10,
9 and 8, respectively.
Explicit physical-lambda names are used by the active pipeline. Historical
`lam1`, `lam3` and `lam6` names mean physical lambda 0.01, 0.03 and 0.06, but
fallback to them is opt-in because the saved workbooks may be from an older
calculation run.

The covariate and MRI-site CSV files remain in the separate historical
neuroHarmonize data project. They are metadata inputs, not turbulence results.

## Validation and execution

From the repository root:

```bash
.venv/bin/python \
  sch1000_N238rev/harmonization_allfeat/run_harmonization.py \
  --exclusions /path/to/authorized_exclusions.csv \
  --validate-only

.venv/bin/python \
  sch1000_N238rev/harmonization_allfeat/run_harmonization.py \
  --exclusions /path/to/authorized_exclusions.csv
```

After changing only the manuscript node scale, preserve the existing global
all-feature table and fit just the lambda=0.01 node table with:

```bash
.venv/bin/python \
  sch1000_N238rev/harmonization_allfeat/run_harmonization.py \
  --exclusions /path/to/authorized_exclusions.csv \
  --node-only
```

Validation checks the four raw group counts (54, 39, 33, 26), PTID uniqueness,
identical participants across both tables, feature counts, finite values, metadata
coverage, single-site filtering, and the final N145 group counts. It does not
fit ComBat or write output.

The full run writes two complete tables to `harmonization_allfeat/recomputed/`
by default: one all-feature table and one Schaefer-1000 node table for
lambda=0.01. Pairwise comparisons must subset these complete tables by `Group`;
they are not separately written or harmonized.

The runner recursively removes MRI sites represented by only one participant,
fits ComBat on N146, and then applies the documented post-ComBat exclusion in
`harmonization_exclusions.csv` to obtain N145. Age, gender, and education enter
as numeric covariates. Group enters categorically using three dummy variables,
with HC Aβ− as the reference, so no linear spacing between disease stages is
assumed.

The post-ComBat exclusion reflects a missing final VBM/GMV output required for
the common multimodal manuscript sample, rather than failed fMRI preprocessing.
The authorized exclusion manifest is deliberately external because it contains
an ADNI participant identifier.

The August 2026 audit found material historical-input differences and an
ambiguous node filename: historical `lam1` meant lambda=0.01, whereas MATLAB
array index 1 is lambda=0.27. The active exporter now uses the explicit
`lambda_0_01` name and regenerates both raw feature tables together.
