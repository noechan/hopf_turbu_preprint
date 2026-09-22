# Machine-learning classification

This module reproduces the four manuscript pairwise classifications using
either empirical turbulent-dynamics features alone or the same features plus
subject-level information capability and susceptibility. The classifier is
selected in the corresponding YAML configuration file.

Participant-level Excel inputs are deliberately not distributed. After
obtaining authorized ADNI data, place the two input workbooks under
`machine_learning/Data/turbu_hopf/` using the filenames declared near the top
of the entry-point scripts.

## Environment

The manuscript environment used Python 3.9:

```bash
python3.9 -m venv .venv
.venv/bin/pip install -r requirements.txt
```

## Run classifications

From `machine_learning/`, run:

```bash
.venv/bin/python \
  studies/final_3d_gs_classification/all_classifiers_run_pipeline_turbu_sch1000.py

.venv/bin/python \
  studies/final_3d_gs_classification/all_classifiers_run_pipeline_turbu_hopf_sch1000_constrained.py
```

Both scripts run 40 repeated stratified 75/25 train/test partitions. Feature
selection and model tuning are performed within training data. In the combined
model, information capability and susceptibility are retained independently of
mRMR selection of the empirical dynamics features.

Generated results are written below `Results/`.

## Accumulate results

First validate that both `LinearSVM` and `PolySVM` have all 40 folds for every
comparison, without creating or replacing a workbook:

```bash
.venv/bin/python \
  studies/final_3d_gs_classification/accumulate_results_turbu_Sch1000.py \
  --check-only

.venv/bin/python \
  studies/final_3d_gs_classification/accumulate_results_turbu_Hopf_Sch1000_constrained.py \
  --check-only
```

If both checks pass, omit `--check-only` to create
`group_comparison_results.xlsx` in the corresponding results directory:

```bash
.venv/bin/python \
  studies/final_3d_gs_classification/accumulate_results_turbu_Sch1000.py

.venv/bin/python \
  studies/final_3d_gs_classification/accumulate_results_turbu_Hopf_Sch1000_constrained.py
```

The workbooks contain one sheet per group comparison. Metric means and standard
deviations remain numeric for downstream analysis, with an additional formatted
`Mean ± SD` column for manuscript-ready reporting.

## Manuscript plots

After the classifications finish, regenerate the overlaid ROC panels with:

```bash
.venv/bin/python \
  studies/final_3d_gs_classification/plots_compare_roc_turbu_vs_hopf.py
```

Regenerate SHAP summaries with:

```bash
.venv/bin/python \
  studies/final_3d_gs_classification/plots_paper_turbu_combat.py
.venv/bin/python \
  studies/final_3d_gs_classification/plots_paper_Hopf_IC_S_constrained.py
```
