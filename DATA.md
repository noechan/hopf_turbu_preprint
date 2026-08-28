# Restricted inputs and expected data layout

Data used in this study were obtained from the Alzheimer's Disease
Neuroimaging Initiative (ADNI). Participant-level raw and derived data are not
redistributed in this repository. Users must obtain ADNI access and comply with
the applicable Data Use Agreement.

Set `ADNI3_ROOT` to the authorized local ADNI3 project directory. The empirical
turbulence pipeline expects the following broad structure:

```text
$ADNI3_ROOT/
├── participants/cross_sectional/
├── timeseries/
│   ├── prepro_fulldenoising/
│   ├── inputs_fulldenoising/
│   ├── outputs_fulldenoising/
│   └── harmonization_inputs/sch1000_N238rev/
└── code/ADNI3_neuroHarmonize_site/data/raw/turbu/
```

The ComBat metadata must contain `PTID`, age, sex/gender, years of education,
MRI site, and amyloid-status group. The manuscript groups contain 51 HC Aβ−,
37 HC Aβ+, 31 MCI Aβ+, and 26 AD Aβ+ participants after harmonization and
the common multimodal-cohort exclusion.
