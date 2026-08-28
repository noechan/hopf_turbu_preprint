import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os
from pathlib import Path

# === File paths ===
python_dir = Path(__file__).resolve().parent.parent
sch1000_root = python_dir.parent.parent
data_path = (
    sch1000_root
    / "harmonization_allfeat"
    / "recomputed"
    / "Turbu_ComBat_ADNI3_allfeatures_N145.xlsx"
)
save_dir = (
    sch1000_root / "figures_N145" / "sch1000" / "Abeta_Status"
    / "harmonized_allfeat" / "python"
)

# Create save directory if it does not exist
os.makedirs(save_dir, exist_ok=True)

# === Load data ===
df = pd.read_excel(data_path)

# === Lambda values and corresponding InfoFlow column names ===
lambda_values = [0.24, 0.21, 0.18, 0.15, 0.12, 0.09, 0.06, 0.03, 0.01]

infoflow_columns = [
    f"InfoFlow_lam_0_{str(lam).split('.')[1].zfill(2)}"
    for lam in lambda_values
]

# === NEW GROUP LABELS ===
rename_map = {
    "HC_ABneg": "HC-",
    "HC_ABpos": "HC+",
    "MCI_ABpos": "MCI+",
    "AD_ABpos": "AD+",
}

df["Group"] = df["Group"].map(rename_map)

group_order = ["HC-", "HC+", "MCI+", "AD+"]

# === Palettes ===
infoflow_palette = sns.color_palette("Oranges", n_colors=len(group_order))

infoflow_colors = dict(zip(group_order, infoflow_palette))

# === Extract InfoFlow data by group ===
grouped_data = {
    group: df[df["Group"] == group][infoflow_columns].to_numpy().T
    for group in group_order
}

# === InfoFlow curve plotting function ===
def plot_infoflow_curve(grouped_data, lambda_values, colors):
    fig, ax = plt.subplots(figsize=(8, 6))

    for group, data in grouped_data.items():
        means = np.mean(data, axis=1)
        stds = np.std(data, axis=1)

        ax.plot(
            lambda_values,
            means,
            "-o",
            color=colors[group],
            label=group,
            markerfacecolor=colors[group],
        )
        ax.fill_between(
            lambda_values,
            means - stds,
            means + stds,
            color=colors[group],
            alpha=0.3,
        )

    ax.set_xlabel("Lambda values")
    ax.set_ylabel("Information flow")
    ax.set_title("Information Cascade Flow Curve (ComBat)")
    ax.legend(title="Group")

    ax.set_xticks(lambda_values)
    ax.set_xticklabels([f"{lam:.2f}" for lam in lambda_values])
    ax.set_xlim(min(lambda_values) - 0.01, max(lambda_values) + 0.01)

    plt.tight_layout()
    return fig, ax


# === Plot InfoFlow curve ===
fig, ax = plot_infoflow_curve(grouped_data, lambda_values, infoflow_colors)
fig.savefig(os.path.join(save_dir, "InfoFlow_Curve_ComBat_4groups.pdf"), dpi=300)
plt.close(fig)

# === Information Cascade Boxplot ===
info_cascade_df = df[["Group", "InfoCascade"]].copy()
info_cascade_df = info_cascade_df[info_cascade_df["Group"].isin(group_order)]

plt.figure(figsize=(6, 6))
sns.boxplot(
    data=info_cascade_df,
    x="Group",
    y="InfoCascade",
    order=group_order,
    hue="Group",
    palette=infoflow_palette,
    dodge=False,
    legend=False,
)

plt.title("Information Cascade by Group (ComBat)")
plt.xlabel("Group")
plt.ylabel("Information Cascade")
plt.tight_layout()

plt.savefig(os.path.join(save_dir, "InfoCascade_Boxplot_ComBat_4groups_sq.pdf"), dpi=300)
plt.close()
