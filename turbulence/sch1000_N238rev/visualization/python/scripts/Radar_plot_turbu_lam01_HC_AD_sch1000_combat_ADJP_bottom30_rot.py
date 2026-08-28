import numpy as np
import matplotlib.pyplot as plt
from scipy.io import loadmat
from collections import Counter
import pandas as pd
import os
from pathlib import Path

# === Define paths ===
python_dir = Path(__file__).resolve().parent.parent
sch1000_root = python_dir.parent.parent
harmonized_figures = (
    sch1000_root / "figures_N145" / "sch1000" / "Abeta_Status"
    / "harmonized_allfeat"
)
nodewise_dir = harmonized_figures / "nodewise_stats"
plot_path = harmonized_figures / "python"

os.makedirs(plot_path, exist_ok=True)

# === Node-wise permutation/FDR file (HC_ABneg vs AD_ABpos) ===
nodewise_file = (
    nodewise_dir
    / "NodeWise_Turbu_lambda_0_01_HC_ABneg_vs_AD_ABpos_permFDR.xlsx"
)

node_df = pd.read_excel(nodewise_file)

required_cols = {"Node", "MeanDiff", "Pvalue", "AdjP", "Significant"}
if not required_cols.issubset(node_df.columns):
    raise ValueError(
        f"Expected columns {required_cols} not all found in {nodewise_file}.\n"
        f"Found columns: {list(node_df.columns)}"
    )

print(f"Node-wise file loaded: {node_df.shape[0]} rows")

# === Load RSN labels (Yeo7vector for Schaefer-1000) ===
rsn_mat = loadmat(sch1000_root / "RSN7vector.mat")
RSN_labels = rsn_mat["Yeo7vector"].squeeze()   # shape: (1000,)

n_nodes = RSN_labels.shape[0]
print(f"RSN_labels length: {n_nodes} nodes")

# === Map 'Node' (e.g., 'Schaefer_1') to RSN labels ===
node_df["Schaefer_idx"] = (
    node_df["Node"]
    .astype(str)
    .str.extract(r"Schaefer_(\d+)", expand=False)
    .astype(int) - 1
)

if node_df["Schaefer_idx"].isna().any():
    raise ValueError("Some 'Node' entries could not be parsed as 'Schaefer_#'.")

if node_df["Schaefer_idx"].max() >= n_nodes or node_df["Schaefer_idx"].min() < 0:
    raise ValueError(
        "Parsed Schaefer indices are out of bounds for RSN_labels. "
        f"Index range: {node_df['Schaefer_idx'].min()}–{node_df['Schaefer_idx'].max()}, "
        f"RSN_labels length: {n_nodes}"
    )

# Attach RSN label (1–7) to each node
node_df["RSN"] = RSN_labels[node_df["Schaefer_idx"].to_numpy()]

# === Select nodes from bottom 30% of AdjP ===
q_thresh = 0.30
adj_p = node_df["AdjP"].to_numpy()

thresh_val = np.quantile(adj_p, q_thresh)
mask_bottom30 = node_df["AdjP"] <= thresh_val
bottom_df = node_df[mask_bottom30].copy()

print(f"Bottom 30% of AdjP: {bottom_df.shape[0]} nodes (q={q_thresh}, thresh={thresh_val:.4g})")

if bottom_df.empty:
    raise ValueError("No nodes found in the bottom 30% of AdjP. Cannot build radar plot.")

# === Count bottom-30% nodes per RSN ===
# Original RSN coding: 1–7 = [VIS, SM, DAT, VAT, LIM, CNT, DMN]
rsn_counts = Counter(bottom_df["RSN"])

values_original = [rsn_counts.get(i, 0) for i in range(1, 8)]
# values_original = [VIS, SM, DAT, VAT, LIM, CNT, DMN]
print("Counts per RSN (VIS, SM, DAT, VAT, LIM, CNT, DMN):", values_original)

RSN_names = ["VIS", "SM", "DAT", "VAT", "LIM", "CNT", "DMN"]
values_plot = values_original

# === Radar plot angles: VIS at top, manuscript order counterclockwise ===
angles = np.linspace(0, 2 * np.pi, len(RSN_names), endpoint=False)
angles = angles.tolist()

# Close the loop for the polygon
values = values_plot + values_plot[:1]
angles_closed = angles + angles[:1]

# === Create radar plot ===
fig, ax = plt.subplots(figsize=(6, 6), subplot_kw=dict(polar=True))
ax.set_theta_offset(np.pi / 2)
ax.set_theta_direction(1)

ax.plot(angles_closed, values, linewidth=2, color='purple')
ax.fill(angles_closed, values, alpha=0.25, color='purple')

ax.set_xticks(angles)
ax.set_xticklabels(RSN_names)

# === Radial ticks: 25, 75, 125, 175 ===
radial_ticks = [25, 50, 75,100, 125]
max_data = max(values_plot) if len(values_plot) > 0 else 0
max_radius = max(max_data, max(radial_ticks))

ax.set_ylim(0, max_radius)

# Only keep ticks up to the max radius
radial_ticks_used = [v for v in radial_ticks if v <= max_radius]
ax.set_yticks(radial_ticks_used)
ax.set_yticklabels([str(v) for v in radial_ticks_used])

# Optional: move radial labels so they do not overlap with DMN axis
ax.set_rlabel_position(225)  # degrees

ax.set_title(
    'Radar λ = 0.01 (node-wise perm/FDR)\n'
    'Bottom 30% AdjP Turbulence Nodes: HC_ABneg vs AD_ABpos',
    pad=20
)
ax.grid(True)

# === Save and show ===
out_file = os.path.join(
    plot_path,
    "Radar_HCneg_vs_ADpos_lambda001_nodewise_permFDR_bottom30AdjP_DMN_topright_purple.pdf"
)

plt.tight_layout()
plt.savefig(out_file)
plt.close(fig)

print(f"Saved radar plot to: {out_file}")
