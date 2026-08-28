# Group-level Hopf workflow

Run the cohort-specific scripts in this order after configuring the visible
local path placeholders:

1. `Compute_Hopf_Freq_*` and `Empirical_corrfcn_*` create the empirical
   frequency and correlation-function inputs.
2. `hopf_DTI_Grange_*_GEC.m` evaluates the global-coupling grid for each of the
   four amyloid-status groups.
3. `get_working_point_G_ADNI3_HC_MCI_AD.m` selects the group working points.
4. `pert_infocapacity_susc_*_GEC.m` runs the perturbation simulations.
5. `get_InfoCap_Susc_ADNI3_HC_MCI_AD.m` aggregates the trial-level model
   outputs; `export_infocap_susc_dataframe.m` prepares aggregate statistics.

The simulation scripts use 100 trials in the manuscript workflow. Generated
G-range, GEC, trial, plot, and export directories are intentionally absent and
ignored by Git.
