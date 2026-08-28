import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path

# ============================================================
# Directories
# ============================================================
python_dir = Path(__file__).resolve().parent.parent
sch1000_root = python_dir.parent.parent
data_dir = python_dir / "data"
plots_dir = (
    sch1000_root / "figures_N145" / "sch1000" / "Abeta_Status"
    / "harmonized_allfeat" / "python"
)

os.makedirs(plots_dir, exist_ok=True)

# ============================================================
# Load ComBat-corrected turbulence data
# ============================================================
combat_file = (
    sch1000_root
    / "harmonization_allfeat"
    / "recomputed"
    / "Turbu_ComBat_ADNI3_allfeatures_N145.xlsx"
)
df = pd.read_excel(combat_file)

# ============================================================
# Rename groups and keep the disease-staging order used throughout N145.
# ============================================================
rename_map = {
    "HC_ABneg": "HC-",
    "HC_ABpos": "HC+",
    "MCI_ABpos": "MCI+",
    "AD_ABpos": "AD+"
}

df["Group"] = df["Group"].map(rename_map)

group_order = ["HC-", "HC+", "MCI+", "AD+"]

# ============================================================
# Purple gradient palette
# ============================================================
purple_gradient = sns.color_palette("Purples", n_colors=len(group_order))
palette = dict(zip(group_order, purple_gradient))

# ============================================================
# Clean scientific plotting (no grid)
# ============================================================
sns.set(style="ticks", context="talk")
sns.set_style({"axes.grid": False})

# ============================================================
# Lambda columns (NO averaging) – original InfoTransfer columns
# ============================================================
lambda_cols = {
    "InfoTransfer_lam_0_01": "0.01",
    "InfoTransfer_lam_0_03": "0.03",
    "InfoTransfer_lam_0_06": "0.06",
    "InfoTransfer_lam_0_09": "0.09",
    "InfoTransfer_lam_0_12": "0.12",
    "InfoTransfer_lam_0_15": "0.15",
    "InfoTransfer_lam_0_18": "0.18",
    "InfoTransfer_lam_0_21": "0.21",
    "InfoTransfer_lam_0_24": "0.24",
    "InfoTransfer_lam_0_27": "0.27",
}

# ============================================================
# Precompute shared y-ticks for 1 - InfoTransfer at λ = 0.01 and 0.03
# ============================================================
low_lam_cols = ["InfoTransfer_lam_0_01", "InfoTransfer_lam_0_03"]

low_vals = df[low_lam_cols].values.astype(float).flatten()
low_vals = low_vals[~np.isnan(low_vals)]

# Transform to 1 - InfoTransfer
low_vals_1m = 1.0 - low_vals

low_min = np.min(low_vals_1m)
low_max = np.max(low_vals_1m)


# ============================================================
# One separate boxplot per λ for 1 - InfoTransfer
# ============================================================
for col_name, lam_label in lambda_cols.items():

    # Extract data for this λ
    df_lam = df[["Group", col_name]].dropna().copy()
    df_lam = df_lam[df_lam["Group"].isin(group_order)]

    # Compute 1 - InfoTransfer
    df_lam["OneMinusInfoTransfer"] = 1.0 - df_lam[col_name].astype(float)

    # Create plot
    plt.figure(figsize=(6, 6))
    sns.boxplot(
        data=df_lam,
        x="Group",
        y="OneMinusInfoTransfer",
        order=group_order,
        hue="Group",
        palette=palette,
        dodge=False,
        legend=False
    )

    ax = plt.gca()

    # Use shared ticks for low λ, otherwise per-plot ticks
    vals = df_lam["OneMinusInfoTransfer"].values
    vmin = np.min(vals)
    vmax = np.max(vals)
    yticks = np.linspace(vmin, vmax, 4)

    ax.set_yticks(yticks)
    ax.set_yticklabels([f"{t:.3f}" for t in yticks])

    # ------------------------------------------------------------

    plt.title(f"1 - Information Transfer by group\nλ = {lam_label}", pad=12)
    plt.xlabel("Group")
    plt.ylabel("1 - Information Transfer")
    plt.xticks(rotation=0)

    sns.despine()
    plt.tight_layout()

    # Output filename
    lam_tag = lam_label.replace(".", "_")
    out_file = os.path.join(
        plots_dir,
        f"OneMinus_Info_Transfer_Boxplot_lam{lam_tag}_ComBat_PurpleGradient.pdf"
    )

    plt.savefig(out_file, dpi=300)
    plt.close()
