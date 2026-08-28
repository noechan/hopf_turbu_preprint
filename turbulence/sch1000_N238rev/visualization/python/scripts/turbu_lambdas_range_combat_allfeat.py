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
# Blue gradient palette
# ============================================================
blue_gradient = sns.color_palette("Blues", n_colors=len(group_order))
palette = dict(zip(group_order, blue_gradient))

# ============================================================
# Clean scientific plotting (no grid)
# ============================================================
sns.set(style="ticks", context="talk")
sns.set_style({"axes.grid": False})

# ============================================================
# Lambda columns (NO averaging)
# ============================================================
lambda_cols = {
    "Turbu_lam_0_01": "0.01",
    "Turbu_lam_0_03": "0.03",
    "Turbu_lam_0_06": "0.06",
    "Turbu_lam_0_09": "0.09",
    "Turbu_lam_0_12": "0.12",
    "Turbu_lam_0_15": "0.15",
    "Turbu_lam_0_18": "0.18",
    "Turbu_lam_0_21": "0.21",
    "Turbu_lam_0_24": "0.24",
    "Turbu_lam_0_27": "0.27",
}

# ============================================================
# Precompute shared y-ticks for λ = 0.01 and 0.03
# ============================================================
low_lam_cols = ["Turbu_lam_0_01", "Turbu_lam_0_03"]
low_vals = df[low_lam_cols].values.flatten()
low_vals = low_vals[~np.isnan(low_vals)]

low_min = np.min(low_vals)
low_max = np.max(low_vals)

# 4 tick positions shared by λ=0.01 and λ=0.03
low_ticks = np.linspace(low_min, low_max, 4)

# ============================================================
# One separate boxplot per λ
# ============================================================
for col_name, lam_label in lambda_cols.items():

    # Extract data for this λ
    df_lam = df[["Group", col_name]].dropna().copy()
    df_lam = df_lam[df_lam["Group"].isin(group_order)]
    df_lam = df_lam.rename(columns={col_name: "Turbulence"})

    # Create plot
    plt.figure(figsize=(6, 6))
    sns.boxplot(
        data=df_lam,
        x="Group",
        y="Turbulence",
        order=group_order,
        hue="Group",
        palette=palette,
        dodge=False,
        legend=False
    )

    ax = plt.gca()

    # ------------------------------------------------------------
    # Set y-ticks:
    #  - For λ = 0.01 and 0.03: shared ticks (low_ticks)
    #  - For others: 4 ticks based on min/max of that λ
    # ------------------------------------------------------------
    if lam_label in ["0.01", "0.03"]:
        ax.set_yticks(low_ticks)
        ax.set_yticklabels([f"{t:.3f}" for t in ax.get_yticks()])

    else:
        vals = df_lam["Turbulence"].values
        vmin = np.min(vals)
        vmax = np.max(vals)
        yticks = np.linspace(vmin, vmax, 4)
        ax.set_yticks(yticks)
        ax.set_yticklabels([f"{t:.3f}" for t in ax.get_yticks()])

    # ------------------------------------------------------------

    plt.title(f"Brain turbulence by group\nλ = {lam_label}", pad=12)
    plt.xlabel("Group")
    plt.ylabel("Turbulence")
    plt.xticks(rotation=0)

    sns.despine()
    plt.tight_layout()

    # Output filename
    lam_tag = lam_label.replace(".", "_")
    out_file = os.path.join(
        plots_dir,
        f"Brain_Turbulence_Boxplot_lam{lam_tag}_ComBat_BlueGradient.pdf"
    )

    plt.savefig(out_file, dpi=300)
    plt.close()
