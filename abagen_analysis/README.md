# Allen Human Brain Atlas analyses

This directory contains the maintained builders for the transcriptomic
excitation--inhibition (E:I) maps and the Alzheimer-risk gene-expression maps
used by the spatial analyses. Downloaded Allen Human Brain Atlas data and
generated expression matrices are deliberately excluded from Git.

## Environment

Use the neuromaps environment from the repository root:

```bash
python3.12 -m venv .venv-neuromaps
.venv-neuromaps/bin/python -m pip install --upgrade pip
.venv-neuromaps/bin/python -m pip install -r neuromaps_analysis/requirements.txt
```

`abagen` downloads its public source data into
`abagen_analysis/abagen-data/` when the builders are first run. This cache can
be several gigabytes and is ignored by Git.

## Primary Schaefer-1000 E:I map

The analysis settings and gene sets are recorded in
`config/ei_sch1000.json`. From the repository root, run:

```bash
.venv-neuromaps/bin/python \
  abagen_analysis/abagen-code/build_ei_sch1000.py
```

The builder estimates expression independently in all 1,000 parcels without
left--right mirroring. It uses RNA-seq-guided probe selection, no sample-wise
normalization, and scaled robust sigmoid gene normalization. The bilateral map
is retained for visualization, while inferential analyses use the 500
left-hemisphere parcels.

## Desikan--Killiany validation map

The paper-style validation settings are recorded in `config/ei_dk68.json`.
Run:

```bash
.venv-neuromaps/bin/python \
  abagen_analysis/abagen-code/build_ei_dk68.py
```

Expression is estimated in the 34 left-hemisphere cortical parcels. Those
values are reflected to the right hemisphere only to create the bilateral
DK68 rendering vector.

## Alzheimer-risk gene maps

`abagen-code/GE_ADrisk_sch1000_2mm.py` builds the 15 Schaefer-1000 maps used
by the gene-expression spin tests. By default, NIfTI files are written below
`abagen-code/gene_niftis/schaefer1000_2mm/`. To use another location, set:

```bash
export AHBA_GENE_NIFTI_DIR=/path/to/authorized/gene_niftis
```

Generated CSV, NIfTI, report, and figure products are excluded from the
publication repository and can be regenerated from the scripts and JSON
configuration files.
