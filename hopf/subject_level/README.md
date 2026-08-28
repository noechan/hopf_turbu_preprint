# Subject-level Hopf workflow

This workflow estimates information capability and susceptibility separately
for every participant. Configure `helper/get_paths.m`, then run:

1. the cohort-specific `Compute_Hopf_Freq_*_SUB_split.m` and
   `Empirical_corrfcn_*_SUB_split.m` scripts;
2. `hopf_DTI_Grange_*_GEC_SUB.m` to fit participant-specific models;
3. the matching `get_optG_*_SUB.m` scripts;
4. `pert_infocapacity_susc_*_GEC_SUB.m` for the perturbation simulations;
5. `get_InfoCap_Susc_ADNI3_HC_MCI_AD_SUB.m`, followed by
   `build_cap_susc_tables_ML.m`, to assemble downstream inputs.

The active N145 group, MOCA, and AT(N) inference is implemented under
`turbulence/sch1000_N238rev/statistical_analysis/perturbation_measures/`, where
age, sex, and years of education are included as nuisance covariates.
Participant-level exports must remain outside Git.
