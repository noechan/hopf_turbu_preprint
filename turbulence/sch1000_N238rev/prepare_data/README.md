# Prepare-data workflow

This directory contains active MATLAB code only. Run stages in numeric order
and stop if any reported subject count differs from the expected N238rev
cohort.

1. `01_cohort_selection/`: select QC-passing PTIDs and CONN IDs separately
   for HC, MCI, and AD.
2. `02_ptid_assembly/`: concatenate batch PTIDs in the same order used by the
   time series and export the combined cohort table.
3. `03_timeseries_extraction/`: extract Schaefer-400, Schaefer-1000,
   Schaefer-100, DBS80, and Glasser-360 data from the extended CONN outputs.
4. `04_combine_modalities/`: combine MPRAGE and IR-FSPGR batches separately
   for HC, MCI, and AD.
5. `05_amyloid_groups/`: retain separate HC Aβ−, HC Aβ+, MCI Aβ+, and AD Aβ+
   groups across all parcellations.

Final-result exports are kept outside this preparation workflow under
`../data_export/`. This includes the neuromaps group-difference exports.

Generated PTID, CONN-ID, CSV, and Excel metadata are stored outside the code
repository at:

```text
$ADNI3_ROOT/timeseries/pipeline_metadata/sch1000_N238rev
```

New extraction writes to `timeseries/inputs_fulldenoising`. The directory
`timeseries/inputs_fulldenoising_clean` is a frozen reference and must not be
overwritten by preparation scripts.

Superseded Schaefer-1000-only extraction and amyloid-split scripts are omitted
from the publication repository.
