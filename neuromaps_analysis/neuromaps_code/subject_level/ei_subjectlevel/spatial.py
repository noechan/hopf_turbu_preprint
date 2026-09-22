"""Project parcel maps to the surface and quantify E:I--turbulence coupling."""

from __future__ import annotations

from pathlib import Path
import tempfile

import nibabel as nib
from neuromaps import nulls, transforms
import numpy as np


def surface_to_array(surface: tuple) -> np.ndarray:
    """Join left and right GIFTI data into one finite one-dimensional array."""

    data = np.hstack(
        [
            np.asarray(hemisphere.agg_data(), dtype=float).squeeze()
            for hemisphere in surface
        ]
    )
    if data.ndim != 1 or not np.isfinite(data).all():
        raise ValueError("Surface map must be a finite one-dimensional array.")
    return data


def project_subjects_to_surface(
    matrix: np.ndarray,
    atlas_img: nib.Nifti1Image,
    atlas: np.ndarray,
    density: str,
) -> np.ndarray:
    """Project every participant's 1,000 parcel values to fsaverage."""

    projected = np.empty((matrix.shape[0], 81924), dtype=np.float32)
    label_masks = [atlas == label for label in range(1, 1001)]
    header = atlas_img.header.copy()
    header.set_data_dtype(np.float32)
    with tempfile.TemporaryDirectory(prefix="ei_subject_surface_") as temp_dir:
        map_path = Path(temp_dir) / "subject.nii.gz"
        for subject_index, values in enumerate(matrix):
            mapped = np.zeros(atlas.shape, dtype=np.float32)
            for mask, value in zip(label_masks, values, strict=True):
                mapped[mask] = value
            nib.save(nib.Nifti1Image(mapped, atlas_img.affine, header), map_path)
            projected[subject_index] = surface_to_array(
                transforms.mni152_to_fsaverage(str(map_path), density)
            )
            if (subject_index + 1) % 25 == 0 or subject_index + 1 == matrix.shape[0]:
                print(f"Projected {subject_index + 1}/{matrix.shape[0]} participants")
    if not np.isfinite(projected).all():
        raise ValueError("Non-finite participant surface values were produced.")
    return projected


def prepare_left_hemisphere_surfaces(
    participant_parcels: np.ndarray,
    atlas_img: nib.Nifti1Image,
    atlas_volume: np.ndarray,
    ei_nifti_path: Path,
    density: str,
) -> tuple[np.ndarray, np.ndarray]:
    """Project both maps and retain the left hemisphere used for inference."""

    participant_surface = project_subjects_to_surface(
        participant_parcels, atlas_img, atlas_volume, density
    )
    ei_surface = surface_to_array(
        transforms.mni152_to_fsaverage(str(ei_nifti_path), density)
    )
    if ei_surface.size % 2:
        raise ValueError("Expected equal left and right surface vertex counts.")
    left_vertices = ei_surface.size // 2
    return participant_surface[:, :left_vertices], ei_surface[:left_vertices]


def _row_unit_center(matrix: np.ndarray) -> np.ndarray:
    """Mean-center and unit-normalize rows so dot products equal Pearson r."""

    centered = matrix - matrix.mean(axis=1, keepdims=True)
    norms = np.linalg.norm(centered, axis=1, keepdims=True)
    if (norms == 0).any() or not np.isfinite(norms).all():
        raise ValueError("Cannot correlate a constant participant map.")
    return centered / norms


def compute_coupling(
    participant_surface: np.ndarray,
    ei_surface: np.ndarray,
    spatial_rotations: int,
    seed: int,
    atlas_name: str,
    density: str,
    hemisphere: str,
    skip_spins: bool,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, int]:
    """Calculate participant Pearson r, Fisher z, and spin-normalized scores."""

    common_mask = ~np.isclose(ei_surface, 0)
    common_mask &= np.all(~np.isclose(participant_surface, 0), axis=0)
    n_vertices = int(common_mask.sum())
    if n_vertices < 1000:
        raise ValueError(f"Too few jointly nonzero surface vertices: {n_vertices}")

    participant_unit = _row_unit_center(
        participant_surface[:, common_mask].astype(float)
    )
    ei_valid = ei_surface[common_mask].astype(float)
    ei_centered = ei_valid - ei_valid.mean()
    ei_unit = ei_centered / np.linalg.norm(ei_centered)
    empirical_r = np.einsum("ij,j->i", participant_unit, ei_unit, optimize=True)
    fisher_z = np.arctanh(np.clip(empirical_r, -1 + 1e-12, 1 - 1e-12))

    if skip_spins:
        nan = np.full(empirical_r.shape, np.nan)
        return empirical_r, fisher_z, nan, nan, n_vertices

    print(
        f"Generating {spatial_rotations} shared Alexander-Bloch rotations "
        f"(seed={seed})"
    )
    full_spin_indices = np.asarray(
        nulls.alexander_bloch(
            None,
            atlas=atlas_name,
            density=density,
            n_perm=spatial_rotations,
            seed=seed,
        )
    )
    if hemisphere != "left":
        raise ValueError("This paper-style analysis requires the left hemisphere.")
    if full_spin_indices.shape != (2 * ei_surface.size, spatial_rotations):
        raise ValueError(
            f"Unexpected bilateral spin-index shape: {full_spin_indices.shape}"
        )
    spin_indices = full_spin_indices[: ei_surface.size]
    if spin_indices.min() < 0 or spin_indices.max() >= ei_surface.size:
        raise ValueError("Left-hemisphere spins reference right-hemisphere vertices.")

    null_r = np.empty(
        (participant_surface.shape[0], spatial_rotations), dtype=np.float32
    )
    chunk_size = 40
    for start in range(0, spatial_rotations, chunk_size):
        stop = min(start + chunk_size, spatial_rotations)
        rotated = ei_surface[spin_indices[:, start:stop]][common_mask].astype(float)
        rotated -= rotated.mean(axis=0, keepdims=True)
        norms = np.linalg.norm(rotated, axis=0, keepdims=True)
        if (norms == 0).any() or not np.isfinite(norms).all():
            raise ValueError("A rotated E:I map is constant or non-finite.")
        null_r[:, start:stop] = np.einsum(
            "ij,jk->ik", participant_unit, rotated / norms, optimize=True
        )
    if not np.isfinite(null_r).all():
        raise ValueError("Participant spatial-null correlations contain non-finite values.")

    null_mean = null_r.mean(axis=1)
    null_sd = null_r.std(axis=1, ddof=1)
    spin_z = (empirical_r - null_mean) / null_sd
    spin_p = (
        np.sum(np.abs(null_r) >= np.abs(empirical_r[:, None]), axis=1) + 1
    ) / (spatial_rotations + 1)
    return empirical_r, fisher_z, spin_z, spin_p, n_vertices
