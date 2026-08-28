# Generative effective connectivity

The four `global_GEC_*_sch1000_mod_eps.m` scripts estimate group-specific
effective-connectivity matrices for HC A$\beta^-$, HC A$\beta^+$,
MCI A$\beta^+$, and AD A$\beta^+$. `hopf_int.m` provides the model integration,
and the two `QC_GEC_*` scripts document the empirical/simulated covariance,
asymmetry, stability, and structural-connectivity quality-control checks.

Configure each script's restricted input paths before running. Fitted matrices
and participant-level source data are generated products and are not included
in the publication repository.
