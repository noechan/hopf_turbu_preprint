# Hopf whole-brain modelling

The manuscript uses two related workflows:

- `group_level/` fits one model per amyloid-status group and evaluates
  perturbation responses across simulation trials;
- `subject_level/` fits participant-specific models to estimate individual
  information capability and susceptibility;
- `gec/` contains the generative effective-connectivity estimation used by the
  Hopf workflows.

The four groups are kept as explicit scripts because HC A$\beta^-$,
HC A$\beta^+$, MCI A$\beta^+$, and AD A$\beta^+$ are distinct analyses.
Participant time series, structural connectivity, fitted matrices, and
simulation outputs are restricted inputs and are not distributed.

Several preprocessing scripts contain `/path/to/...` placeholders. Replace
only those paths with authorized local locations; do not commit the resulting
participant-level files. The downstream covariate-adjusted tests are under
`turbulence/sch1000_N238rev/statistical_analysis/perturbation_measures/`.
