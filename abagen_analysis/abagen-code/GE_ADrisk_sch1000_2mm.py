import nibabel as nib
import numpy as np
import pandas as pd
from nilearn import datasets
import abagen
import os
from pathlib import Path

script_dir = Path(__file__).resolve().parent
output_dir = Path(
    os.environ.get(
        "AHBA_GENE_NIFTI_DIR",
        script_dir / "gene_niftis" / "schaefer1000_2mm",
    )
).expanduser().resolve()
output_dir.mkdir(parents=True, exist_ok=True)
# Fetch the Allen Human Brain Atlas (AHBA) microarray data
# files = abagen.fetch_microarray(donors='all', verbose=0)

# Load the Schaefer100 (7 networks, 1mm resolution) atlas
schaefer = datasets.fetch_atlas_schaefer_2018(n_rois=1000, yeo_networks=7, resolution_mm=2)
atlas_img = nib.load(schaefer.maps)  # Load atlas as NIfTI
atlas_data = atlas_img.get_fdata()  # Get atlas 3D array

# Extract gene expression mapped to the Schaefer1000 atlas
expression_data, report = abagen.get_expression_data(
    atlas=atlas_img,  # Pass the NIfTI file
    ibf_threshold=0.5,  # Standard thresholding
    missing="centroids",  # Missing data parameter
    norm_matched=False,  # Recommended when using missing data
    return_report=True,
)



expression_data.to_csv(script_dir / "full_GE_schaefer1000_2mm.csv")
gene_names = expression_data.columns.tolist()
pd.Series(gene_names).to_csv(
    script_dir / "full_gene_symbols_schaefer1000_2mm.csv",
    index=False,
    header=["GeneSymbol"],
)
print(type(report))
with (script_dir / "abagen_report_schaefer1000_2mm.txt").open(
    "w", encoding="utf-8"
) as f:
    f.write(report)


# Define genes of interest

genes_of_interest = [
    "APP", "ADAM10", "BACE1", "PSEN1", "PSEN2",
    "APOE","CLU","SORL1",
    "ABCA7","BIN1","CD2AP","CD33","PTK2B","PICALM","RIN3",
]

# Filter the gene expression data for selected genes
gene_expression_data = expression_data[genes_of_interest]
gene_expression_data.to_csv(script_dir / "gene_expression_AD_schaefer1000_2mm.csv")



# Function to generate NIfTI for a given gene
def create_gene_nifti(gene_name, atlas_data, atlas_img):
    """Generates a NIfTI image for a given gene expression profile."""
    nifti_data = np.zeros_like(atlas_data)  # Empty brain template

    # Map gene expression to brain regions
    for roi_idx in range(1, 1001):  # Schaefer1000 has 1000 ROIs, indexed 1-1001
        mask = atlas_data == roi_idx  # Find voxels belonging to this ROI
        if roi_idx in gene_expression_data.index:
            nifti_data[mask] = gene_expression_data.loc[roi_idx, gene_name]  # Assign expression value

    # Save NIfTI
    nifti_img = nib.Nifti1Image(nifti_data, affine=atlas_img.affine)
    output_path = output_dir / f"{gene_name}_expression_sch1000_2mm.nii.gz"
    nib.save(nifti_img, output_path)
    return nifti_img


# Generate NIfTI images for all genes of interest
for gene in genes_of_interest:
    if gene in gene_expression_data.columns:
        create_gene_nifti(gene, atlas_data, atlas_img)

print("NIfTI images for selected genes have been generated.")
