#!/usr/bin/env python3
"""Recompute the manuscript ComBat tables from mounted ADNI inputs.

The script never overwrites active tables by default. It writes to a
``recomputed`` subdirectory unless ``--output-dir`` is supplied explicitly.
"""

from __future__ import annotations

import argparse
import os
from pathlib import Path

import numpy as np
import pandas as pd
from neuroHarmonize import harmonizationLearn


HERE = Path(__file__).resolve().parent
DEFAULT_ADNI3_ROOT = Path(
    os.environ.get(
        "ADNI3_ROOT",
        "/path/to/ADNI3",
    )
)

EXPECTED_GROUP_COUNTS = {
    "HC_ABneg": 51,
    "HC_ABpos": 37,
    "MCI_ABpos": 31,
    "AD_ABpos": 26,
}

EXPECTED_RAW_GROUP_COUNTS = {
    "HC_ABneg": 54,
    "HC_ABpos": 39,
    "MCI_ABpos": 33,
    "AD_ABpos": 26,
}

GROUP_LEVELS = tuple(EXPECTED_RAW_GROUP_COUNTS)
GROUP_REFERENCE = "HC_ABneg"
REQUIRED_DEMOGRAPHIC_COLUMNS = {"PTID", "age", "gender", "edu"}
REQUIRED_SITE_COLUMNS = {"PTID", "Modality", "ScanSite"}

NODE_SCALE_INPUTS = {
    "0.01": (
        "turbu_by_node_raw_lambda_0_01_sch1000_HC_MCI_AD_ABeta_precombat.xlsx",
        "turbu_by_node_raw_lam1_sch1000_HC_MCI_AD_ABeta_precombat.xlsx",
    ),
    "0.03": (
        "turbu_by_node_raw_lambda_0_03_sch1000_HC_MCI_AD_ABeta_precombat.xlsx",
        "turbu_by_node_raw_lam3_sch1000_HC_MCI_AD_ABeta_precombat.xlsx",
    ),
    "0.06": (
        "turbu_by_node_raw_lambda_0_06_sch1000_HC_MCI_AD_ABeta_precombat.xlsx",
        "turbu_by_node_raw_lam6_sch1000_HC_MCI_AD_ABeta_precombat.xlsx",
    ),
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--adni3-root", type=Path, default=DEFAULT_ADNI3_ROOT)
    parser.add_argument(
        "--feature-dir",
        type=Path,
        help=(
            "Directory containing the four raw feature workbooks. Defaults to "
            "$ADNI3_ROOT/timeseries/harmonization_inputs/sch1000_N238rev."
        ),
    )
    parser.add_argument(
        "--metadata-dir",
        type=Path,
        help=(
            "Directory containing the covariate and site CSV files. Defaults to "
            "$ADNI3_ROOT/code/ADNI3_neuroHarmonize_site/data/raw/turbu."
        ),
    )
    parser.add_argument("--output-dir", type=Path, default=HERE / "recomputed")
    parser.add_argument(
        "--exclusions",
        type=Path,
        required=True,
        help=(
            "Authorized local CSV containing a PTID column with the "
            "post-ComBat manuscript exclusion. Participant identifiers are "
            "not distributed with this repository."
        ),
    )
    parser.add_argument(
        "--validate-only",
        action="store_true",
        help="Validate inputs and cohort filtering without fitting ComBat or writing files.",
    )
    parser.add_argument(
        "--node-only",
        action="store_true",
        help=(
            "Fit only the node tables selected by --node-scales; keep the "
            "existing all-feature harmonization unchanged."
        ),
    )
    parser.add_argument(
        "--node-scales",
        nargs="+",
        choices=tuple(NODE_SCALE_INPUTS),
        default=["0.01"],
        help=(
            "Physical lambda values to harmonize (default: 0.01). The "
            "Neuromaps replication uses: 0.01 0.03 0.06."
        ),
    )
    parser.add_argument(
        "--allow-legacy-input-names",
        action="store_true",
        help=(
            "Allow fallback to historical lam1/lam3/lam6 filenames when the "
            "explicit physical-lambda file is absent. Disabled by default "
            "because legacy files may come from an older calculation run."
        ),
    )
    return parser.parse_args()


def resolve_node_input(
    feature_dir: Path,
    scale: str,
    allow_legacy_input_names: bool,
) -> Path:
    physical_name, legacy_name = NODE_SCALE_INPUTS[scale]
    physical_path = feature_dir / physical_name
    if physical_path.is_file():
        return physical_path
    legacy_path = feature_dir / legacy_name
    if allow_legacy_input_names and legacy_path.is_file():
        print(
            "WARNING: using historical input name "
            f"{legacy_path.name} for physical lambda={scale}."
        )
        return legacy_path
    legacy_note = ""
    if legacy_path.is_file():
        legacy_note = (
            f" A historical file exists at {legacy_path}, but fallback is "
            "disabled; rerun data_export/export_harmonization_inputs.m to "
            "create an explicit physical-lambda input, or opt in with "
            "--allow-legacy-input-names after verifying its provenance."
        )
    raise FileNotFoundError(
        f"No node input found for physical lambda={scale}: expected "
        f"{physical_path}.{legacy_note}"
    )


def read_exclusions(path: Path) -> set[str]:
    table = pd.read_csv(path)
    if "PTID" not in table.columns:
        raise ValueError(f"Exclusion manifest has no PTID column: {path}")
    return set(table["PTID"].dropna().astype(str))


def validate_metadata(
    demographics: pd.DataFrame,
    sites: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    missing_demographics = REQUIRED_DEMOGRAPHIC_COLUMNS - set(demographics.columns)
    if missing_demographics:
        raise ValueError(
            "Missing demographic columns: "
            f"{sorted(missing_demographics)}"
        )
    missing_sites = REQUIRED_SITE_COLUMNS - set(sites.columns)
    if missing_sites:
        raise ValueError(f"Missing site columns: {sorted(missing_sites)}")

    demographics = demographics.assign(PTID=demographics["PTID"].astype(str))
    sites = sites.assign(PTID=sites["PTID"].astype(str))
    if demographics["PTID"].duplicated().any():
        duplicates = demographics.loc[
            demographics["PTID"].duplicated(), "PTID"
        ].tolist()
        raise ValueError(f"Duplicate demographic PTIDs: {duplicates}")

    mri_sites = sites[sites["Modality"] == "MRI"]
    if mri_sites["PTID"].duplicated().any():
        duplicates = mri_sites.loc[mri_sites["PTID"].duplicated(), "PTID"].tolist()
        raise ValueError(f"Duplicate MRI-site PTIDs: {duplicates}")
    return demographics, sites


def filter_common_and_single_site(
    sites: pd.DataFrame,
    demographics: pd.DataFrame,
    features: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    common_ids = (
        set(sites["PTID"].astype(str))
        & set(demographics["PTID"].astype(str))
        & set(features["PTID"].astype(str))
    )

    sites = sites.assign(PTID=sites["PTID"].astype(str))
    demographics = demographics.assign(PTID=demographics["PTID"].astype(str))
    features = features.assign(PTID=features["PTID"].astype(str))

    while True:
        current = sites[sites["PTID"].isin(common_ids)]
        single_site_ids = (
            current[current["Modality"] == "MRI"]
            .groupby("ScanSite")
            .filter(lambda rows: len(rows) == 1)["PTID"]
            .drop_duplicates()
            .tolist()
        )
        if not single_site_ids:
            break
        print(f"Removing single-site MRI subjects: {single_site_ids}")
        common_ids.difference_update(single_site_ids)

    restrict = lambda frame: frame[frame["PTID"].isin(common_ids)].reset_index(drop=True)
    return restrict(sites), restrict(demographics), restrict(features)


def build_covariates(
    features: pd.DataFrame,
    demographics: pd.DataFrame,
    sites: pd.DataFrame,
) -> pd.DataFrame:
    covariates = features[["PTID", "Group"]].copy()
    covariates = covariates.merge(
        demographics[["PTID", "age", "gender", "edu"]], on="PTID", how="left"
    )
    mri_sites = sites[sites["Modality"] == "MRI"]
    covariates = covariates.merge(
        mri_sites[["PTID", "ScanSite"]], on="PTID", how="left"
    )
    if len(covariates) != len(features):
        raise ValueError(
            "Covariate merges changed the participant count: "
            f"{len(features)} feature rows became {len(covariates)} covariate rows"
        )

    group = pd.Categorical(
        covariates.pop("Group"), categories=GROUP_LEVELS, ordered=False
    )
    if group.isna().any():
        raise ValueError("Unknown or missing group labels in the feature table")
    group_dummies = pd.get_dummies(group, prefix="group", dtype=float).drop(
        columns=f"group_{GROUP_REFERENCE}"
    )

    covariates = pd.concat([covariates, group_dummies], axis=1)
    covariates = covariates.rename(columns={"ScanSite": "SITE"}).drop(columns="PTID")
    for column in ("age", "gender", "edu"):
        covariates[column] = pd.to_numeric(covariates[column], errors="raise")

    if covariates.isna().any().any():
        raise ValueError(
            "Missing ComBat covariates:\n" + covariates.isna().sum().to_string()
        )
    print(
        "ComBat covariates: age, gender, education, "
        f"categorical group (reference={GROUP_REFERENCE}), and SITE"
    )
    return covariates


def harmonize_one(
    input_file: Path,
    output_file: Path,
    demographics: pd.DataFrame,
    sites: pd.DataFrame,
    excluded_ptids: set[str],
) -> pd.DataFrame:
    features = pd.read_excel(input_file, engine="openpyxl")
    required = {"PTID", "Group"}
    if not required.issubset(features.columns):
        raise ValueError(f"{input_file} must contain PTID and Group columns")

    sites_filtered, demographics_filtered, features_filtered = (
        filter_common_and_single_site(sites, demographics, features)
    )
    covariates = build_covariates(
        features_filtered, demographics_filtered, sites_filtered
    )
    feature_columns = [c for c in features_filtered.columns if c not in required]
    matrix = features_filtered[feature_columns].to_numpy(dtype=float)

    print(
        f"Fitting ComBat: {matrix.shape[0]} participants x "
        f"{matrix.shape[1]} features; {covariates['SITE'].nunique()} MRI sites"
    )
    _, harmonized = harmonizationLearn(
        matrix,
        covariates,
        eb=True,
        smooth_terms=[],
        ref_batch=None,
    )
    if not np.isfinite(harmonized).all():
        raise ValueError(f"Non-finite harmonized values produced for {input_file}")

    result = features_filtered.copy()
    result[feature_columns] = harmonized

    # Preserve the manuscript workflow exactly: fit ComBat on the 146
    # multi-site subjects, then apply the separately documented exclusion.
    # Excluding before ComBat changes every retained harmonized value.
    if excluded_ptids:
        found = excluded_ptids & set(result["PTID"])
        missing_exclusions = excluded_ptids - found
        if missing_exclusions:
            raise ValueError(
                "Configured post-ComBat exclusions were not present after site "
                f"filtering: {sorted(missing_exclusions)}"
            )
        print(f"Applying post-ComBat manuscript exclusions: {sorted(found)}")
        result = result[~result["PTID"].isin(excluded_ptids)].reset_index(drop=True)

    validate_groups(result, output_file.name)

    output_file.parent.mkdir(parents=True, exist_ok=True)
    result.to_excel(output_file, index=False)
    print(f"Saved {result.shape}: {output_file}")
    return result


def validate_raw_table(
    path: Path,
    expected_feature_count: int,
) -> pd.DataFrame:
    table = pd.read_excel(path, engine="openpyxl")
    required = {"PTID", "Group"}
    if not required.issubset(table.columns):
        raise ValueError(f"{path} must contain PTID and Group columns")

    table = table.assign(PTID=table["PTID"].astype(str))
    if table["PTID"].duplicated().any():
        duplicates = table.loc[table["PTID"].duplicated(), "PTID"].tolist()
        raise ValueError(f"Duplicate PTIDs in {path}: {duplicates}")

    counts = table["Group"].value_counts().to_dict()
    if counts != EXPECTED_RAW_GROUP_COUNTS:
        raise ValueError(
            f"Unexpected raw group counts for {path.name}: {counts}; "
            f"expected {EXPECTED_RAW_GROUP_COUNTS}"
        )

    feature_columns = [column for column in table.columns if column not in required]
    if len(feature_columns) != expected_feature_count:
        raise ValueError(
            f"{path.name} has {len(feature_columns)} features; "
            f"expected {expected_feature_count}"
        )
    matrix = table[feature_columns].to_numpy(dtype=float)
    if not np.isfinite(matrix).all():
        raise ValueError(f"Non-finite raw values in {path}")

    print(f"Validated raw input {table.shape}: {path}")
    return table


def validate_groups(table: pd.DataFrame, label: str) -> None:
    counts = table["Group"].value_counts().to_dict()
    if counts != EXPECTED_GROUP_COUNTS:
        raise ValueError(
            f"Unexpected group counts for {label}: {counts}; "
            f"expected {EXPECTED_GROUP_COUNTS}"
        )


def main() -> None:
    args = parse_args()
    feature_dir = args.feature_dir or (
        args.adni3_root
        / "timeseries"
        / "harmonization_inputs"
        / "sch1000_N238rev"
    )
    metadata_dir = args.metadata_dir or (
        args.adni3_root
        / "code"
        / "ADNI3_neuroHarmonize_site"
        / "data"
        / "raw"
        / "turbu"
    )
    demographics = pd.read_csv(metadata_dir / "covariates_ADNI3_ABeta_N152.csv")
    sites = pd.read_csv(metadata_dir / "sites_ADNI3_ABeta_N304.csv")
    demographics, sites = validate_metadata(demographics, sites)
    excluded_ptids = read_exclusions(args.exclusions)

    allfeatures_path = feature_dir / "turbu_raw_ADNI3_ABeta_N152_sch1000.xlsx"
    node_paths = {
        scale: resolve_node_input(
            feature_dir, scale, args.allow_legacy_input_names
        )
        for scale in dict.fromkeys(args.node_scales)
    }
    raw_tables = [validate_raw_table(allfeatures_path, 32)]
    raw_tables.extend(
        validate_raw_table(path, 1000) for path in node_paths.values()
    )
    reference_ptids = set(raw_tables[0]["PTID"])
    for table in raw_tables[1:]:
        if set(table["PTID"]) != reference_ptids:
            raise ValueError("Raw feature tables do not contain identical PTID sets")

    if args.validate_only:
        for table in raw_tables:
            sites_filtered, demographics_filtered, features_filtered = (
                filter_common_and_single_site(sites, demographics, table)
            )
            build_covariates(
                features_filtered, demographics_filtered, sites_filtered
            )
            retained = features_filtered[
                ~features_filtered["PTID"].isin(excluded_ptids)
            ].reset_index(drop=True)
            validate_groups(retained, "validate-only")
        print("Validation passed; no ComBat model was fitted and no files were written.")
        return

    if not args.node_only:
        harmonize_one(
            allfeatures_path,
            args.output_dir / "Turbu_ComBat_ADNI3_allfeatures_N145.xlsx",
            demographics,
            sites,
            excluded_ptids,
        )

    for scale, input_path in node_paths.items():
        scale_label = scale.replace(".", "_")
        harmonize_one(
            input_path,
            args.output_dir
            / "Turbu_ComBat_ADNI3_HC_MCI_AD_ABeta_"
            f"lambda_{scale_label}_sch1000_N145.xlsx",
            demographics,
            sites,
            excluded_ptids,
        )


if __name__ == "__main__":
    main()
