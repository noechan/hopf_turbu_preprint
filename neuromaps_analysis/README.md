# Neuromaps analyses

This directory contains the publication workflows for spatially relating the
N145 Schaefer-1000 turbulence difference maps to Neurosynth memory maps and to
15 a priori Alzheimer-related AHBA gene-expression maps. All scripts resolve
repository files relative to their own location.

The historical script suffixes `lam1`, `lam3`, and `lam6` denote the physical
values $\lambda=0.01$, $0.03$, and $0.06$, respectively. The turbulence maps
are absolute group-mean differences, consistent with the manuscript analysis.

## Environment

```bash
python3.12 -m venv .venv
.venv/bin/python -m pip install --upgrade pip
.venv/bin/python -m pip install -r requirements.txt
```

## Rebuild the aggregate turbulence NIfTI maps

The six non-identifiable parcel-level `.mat` inputs and the Schaefer atlas are
included. Rebuild their MNI152 NIfTI representations with:

```bash
.venv/bin/python \
  neuromaps_code/parcellate/build_turbulence_difference_niftis.py
```

## Neurosynth spatial correlations

The three thresholded, parcellated Neurosynth maps used in the planned
analysis are included. The following command evaluates three terms at three
scales for each of the two contrasts. Alexander--Bloch spin tests use 1,000
rotations by default, and Benjamini--Hochberg FDR is applied across the nine
tests within each contrast.

```bash
.venv/bin/python \
  neuromaps_code/correlations_nulls/run_neurosynth_correlations_N145.py
```

Use `--contrast hc_ad` or `--contrast mci_ad` to run one contrast. Results are
written to `results/neurosynth/`, which is ignored by Git.

## AHBA gene-expression correlations

The processed gene maps are not redistributed. Set `AHBA_GENE_NIFTI_DIR` to a
directory containing files named `<GENE>_expression_sch1000_2mm.nii.gz`, then
run the desired contrast and scale, for example:

```bash
AHBA_GENE_NIFTI_DIR=/path/to/authorized/gene_niftis \
  .venv/bin/python \
  neuromaps_code/correlations_nulls/AHBA/corr_geneexpr_Turbu_GE15_HC_AD_lam1_fsaverage41k_sch1000.py
```

The six scripts cover HC$^-$ vs. AD$^+$ and MCI$^+$ vs. AD$^+$ at the three
scales. Each uses fsaverage 41k projection, Pearson correlation, and 1,000
Alexander--Bloch rotations with seed 1234. Outputs are written to
`results/gene_spatial_correlations/` and ignored by Git.

## Data provenance

- Turbulence difference maps were derived from the current N145 ComBat tables.
- Neurosynth maps are redistributed only as the derived thresholded and
  parcellated inputs needed to reproduce this analysis.
- AHBA expression maps must be regenerated or supplied locally; see the main
  repository `DATA.md` and `THIRD_PARTY_NOTICES.md`.
