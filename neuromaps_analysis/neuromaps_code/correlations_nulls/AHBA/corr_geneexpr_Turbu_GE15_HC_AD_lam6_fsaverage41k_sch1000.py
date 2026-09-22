
import os
from pathlib import Path
import tempfile

# Resolve all project files from neuromaps_analysis, independently of the
# directory from which this script is launched.
NEUROMAPS_ROOT = Path(__file__).resolve().parents[3]
NEUROMAPS_DATA = NEUROMAPS_ROOT / "neuromaps-data"
DEFAULT_CACHE = NEUROMAPS_DATA / "cache"
os.environ.setdefault("NEUROMAPS_DATA", str(DEFAULT_CACHE))
os.environ.setdefault(
    "MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "neuromaps_matplotlib")
)

# Import dependencies after NEUROMAPS_DATA is configured.
from neuromaps import transforms, nulls
from neuromaps.stats import compare_images
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# Set up base path and lam values
base_path = Path(os.environ.get(
    "AHBA_GENE_NIFTI_DIR",
    NEUROMAPS_ROOT.parent
    / "abagen_analysis/abagen-code/gene_niftis/schaefer1000_2mm",
)).expanduser()
output_dir = NEUROMAPS_ROOT / "results" / "gene_spatial_correlations"
output_dir.mkdir(parents=True, exist_ok=True)
gene_names = [ "APP", "ADAM10","BACE1","PSEN1", "PSEN2",
               "APOE", "CLU", "SORL1",
               "ABCA7","BIN1","CD2AP", "CD33", "PTK2B","PICALM", "RIN3"]

# Load turbu image and project to fsLR 32k
turbu_MNI152_2mm_path = (
    NEUROMAPS_DATA / "annotations" / "turbu" / "MNI152" / "N145_sch1000"
    / "turbu_diff_hc_ad_lam6_sch1000_N145_2mm.nii.gz"
)


def surface_to_array(surface):
    """Concatenate in-memory GIFTI hemispheres (neuromaps 0.0.5 compatible)."""
    return np.hstack([
        np.asarray(hemisphere.agg_data(), dtype=float).squeeze()
        for hemisphere in surface
    ])


if not turbu_MNI152_2mm_path.is_file():
    raise FileNotFoundError(f"Missing turbulence map: {turbu_MNI152_2mm_path}")

turbu_fsaverage41k = transforms.mni152_to_fsaverage(
    str(turbu_MNI152_2mm_path), "41k"
)
turbu_fsaverage41k_data = surface_to_array(turbu_fsaverage41k)

print("turbu_fsaverage41k_data shape:", turbu_fsaverage41k_data.shape)

turbu_rotated = nulls.alexander_bloch(turbu_fsaverage41k_data, atlas='fsaverage', density='41k',
                               n_perm=1000, seed=1234)
print(turbu_rotated.shape)

# Initialize dictionary to store correlation results
results = {}

# Loop through each gene NIfTI image
for gene in gene_names:
    gene_path = base_path / f'{gene}_expression_sch1000_2mm.nii.gz'
    if not gene_path.is_file():
        raise FileNotFoundError(f"Missing gene-expression map: {gene_path}")

    # Project gene image to fsaverage 41k
    gene_fsaverage41k = transforms.mni152_to_fsaverage(str(gene_path), '41k')
    gene_fsaverage41k_data = surface_to_array(gene_fsaverage41k)

    # Compute Pearson's r using precomputed surrogates
    r, p, nulldist = compare_images(turbu_fsaverage41k_data, gene_fsaverage41k_data,
                                    nulls=turbu_rotated, metric='pearsonr',
                                    return_nulls=True)
    print(f'{gene}: rho = {r}, pspin = {p}')

    # Store results
    results[gene] = {'r': r, 'p': p, 'nulldist': nulldist}


# Prepare a list to hold all rows
rows = []

for gene_key, stats in results.items():
    gene = gene_key
    r_val = stats['r']
    p_val = stats['p']
    nulls = stats['nulldist']

    # Add one row for each null value (for full export)
    for i, null_r in enumerate(nulls):
        rows.append({
            'Gene': gene,
            'Empirical_r': r_val if i == 0 else '',
            'p_spin': p_val if i == 0 else '',
            'Null_r': null_r
        })

# Convert to DataFrame
df = pd.DataFrame(rows)

# Save to CSV
csv_path = output_dir / (
    "correlation_results_with_nulls_genesAD_turbulam6_HC_AD_"
    "fsaverage41k_GE15_sch1000_v2.csv"
)
df.to_csv(csv_path, index=False)
print("Saved CSV with r, p, and null distributions.")


# Set up plot
fig, ax = plt.subplots(figsize=(10, 6))

# Extract data
positions = list(range(1, len(gene_names) + 1))
null_dists = [results[gene]['nulldist'] for gene in gene_names]
empirical_rs = [results[gene]['r'] for gene in gene_names]

# Create boxplots
ax.boxplot(null_dists, positions=positions, patch_artist=True, boxprops=dict(facecolor='lightgray'))

# Overlay scatter plot for empirical r values
ax.scatter(positions, empirical_rs, color='orange', marker='o', s=60, label='Empirical r', zorder=3)

# Format plot
ax.set_ylabel("Pearson's r")
ax.set_xticks(positions)
ax.set_xticklabels([f"{gene}" for gene in gene_names], rotation=45, ha='right')
ax.set_title("Correlation between Turbulence Lambda 6 and Gene expression (fsaverage 41k)")
ax.legend()
plt.tight_layout()
figure_path = output_dir / (
    "correlation_plot_genesAD_turbulam6_HC_AD_"
    "fsaverage41k_GE15_sch1000_v2.pdf"
)
plt.savefig(figure_path, dpi=300)
plt.close(fig)
