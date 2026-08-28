# Neuromaps analyses

This directory contains the publication workflows for spatially relating the
N145 Schaefer-1000 turbulence difference maps to Neurosynth memory maps and to
15 a priori Alzheimer-related AHBA gene-expression maps. All scripts resolve
repository files relative to their own location.

The script suffixes `lam1`, `lam3`, and `lam6` denote the physical
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
- AHBA expression maps must be regenerated or supplied locally.
