#!/usr/bin/env python3
"""Summarise the strongest covariate-adjusted node effects by Yeo-7 network.

The manuscript rule is preserved exactly: within each planned contrast,
select every parcel whose BH-adjusted permutation p-value is at or below the
30th percentile. Ties at the percentile threshold are retained, so the
selected set can contain more than 300 of the 1,000 Schaefer parcels.
"""

from __future__ import annotations

import argparse
import math
import os
from pathlib import Path
import tempfile

os.environ.setdefault(
    "MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "turbu_n145_matplotlib")
)
os.environ.setdefault("MPLBACKEND", "Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.io import loadmat


SCRIPT_DIR = Path(__file__).resolve().parent
SCH1000_ROOT = SCRIPT_DIR.parent.parent
RESULTS_DIR = SCRIPT_DIR / "results" / "N145_nodewise_metastability"
DEFAULT_OUTPUT_DIR = (
    SCH1000_ROOT
    / "figures_N145"
    / "sch1000"
    / "Abeta_Status"
    / "harmonized_allfeat"
    / "RSN_top30_FDR_age_sex_education"
)

CONTRASTS = (
    {
        "id": "HC_ABneg_vs_AD_ABpos",
        "label": r"HC$^-$ vs. AD$^+$",
        "input": RESULTS_DIR
        / "N145_nodewise_metastability_HC_ABneg_vs_AD_ABpos_"
        "FreedmanLane_age_sex_education.csv",
    },
    {
        "id": "MCI_ABpos_vs_AD_ABpos",
        "label": r"MCI$^+$ vs. AD$^+$",
        "input": RESULTS_DIR
        / "N145_nodewise_metastability_MCI_ABpos_vs_AD_ABpos_"
        "FreedmanLane_age_sex_education.csv",
    },
)

RSN_NAMES = {
    1: "VIS",
    2: "SM",
    3: "DAT",
    4: "VAT",
    5: "LIM",
    6: "CNT",
    7: "DMN",
}
REQUIRED_COLUMNS = {
    "Node",
    "Adjusted_Group_Beta",
    "P_Permutation",
    "P_FDR_BH",
    "Significant_FDR_0_05",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=DEFAULT_OUTPUT_DIR,
        help=f"Destination for tables and figures (default: {DEFAULT_OUTPUT_DIR})",
    )
    parser.add_argument(
        "--selection-fraction",
        type=float,
        default=0.30,
        help="Lower fraction of BH-adjusted p-values to retain (default: 0.30)",
    )
    parser.add_argument(
        "--mapping-file",
        type=Path,
        default=SCH1000_ROOT / "RSN7vector.mat",
        help="MAT file containing Yeo7vector for Schaefer-1000",
    )
    parser.add_argument(
        "--validate-only",
        action="store_true",
        help="Validate the inputs and print selection counts without saving outputs",
    )
    return parser.parse_args()


def load_mapping(path: Path) -> np.ndarray:
    if not path.is_file():
        raise FileNotFoundError(f"Yeo-7 mapping not found: {path}")
    mat = loadmat(path)
    if "Yeo7vector" not in mat:
        raise KeyError(f"Variable 'Yeo7vector' not found in {path}")
    labels = np.asarray(mat["Yeo7vector"]).squeeze()
    if labels.shape != (1000,):
        raise ValueError(f"Expected 1,000 Yeo-7 labels; found shape {labels.shape}")
    if not np.all(np.isfinite(labels)) or not np.all(labels == labels.astype(int)):
        raise ValueError("Yeo-7 labels must be finite integers")
    labels = labels.astype(int)
    if set(np.unique(labels)) != set(RSN_NAMES):
        raise ValueError(
            f"Expected Yeo-7 labels 1..7; found {sorted(np.unique(labels).tolist())}"
        )
    return labels


def load_node_results(path: Path, labels: np.ndarray) -> pd.DataFrame:
    if not path.is_file():
        raise FileNotFoundError(f"Node-wise result not found: {path}")
    data = pd.read_csv(path)
    missing = REQUIRED_COLUMNS.difference(data.columns)
    if missing:
        raise ValueError(f"Missing columns in {path}: {sorted(missing)}")
    if len(data) != 1000:
        raise ValueError(f"Expected 1,000 node rows in {path}; found {len(data)}")

    extracted = data["Node"].astype(str).str.extract(r"^Schaefer_(\d+)$", expand=False)
    if extracted.isna().any():
        bad = data.loc[extracted.isna(), "Node"].head().tolist()
        raise ValueError(f"Unparseable node identifiers in {path}: {bad}")
    data["Schaefer_Index"] = extracted.astype(int)
    expected = set(range(1, 1001))
    observed = set(data["Schaefer_Index"].tolist())
    if observed != expected or data["Schaefer_Index"].duplicated().any():
        raise ValueError("Node identifiers must contain Schaefer_1 through Schaefer_1000 once")

    for column in ("Adjusted_Group_Beta", "P_Permutation", "P_FDR_BH"):
        data[column] = pd.to_numeric(data[column], errors="raise")
        if not np.all(np.isfinite(data[column])):
            raise ValueError(f"Column {column} contains non-finite values in {path}")
    for column in ("P_Permutation", "P_FDR_BH"):
        if not data[column].between(0, 1, inclusive="both").all():
            raise ValueError(f"Column {column} must lie in [0, 1] in {path}")

    data = data.sort_values("Schaefer_Index").reset_index(drop=True)
    data["Yeo7_ID"] = labels
    data["Yeo7_Network"] = data["Yeo7_ID"].map(RSN_NAMES)
    return data


def analyse_contrast(
    config: dict[str, object],
    labels: np.ndarray,
    selection_fraction: float,
) -> tuple[pd.DataFrame, pd.DataFrame, dict[str, object]]:
    source = Path(config["input"])
    data = load_node_results(source, labels)
    threshold = float(data["P_FDR_BH"].quantile(selection_fraction))
    data["Selected_Top30"] = data["P_FDR_BH"] <= threshold
    data["FDR_Rank_Min"] = data["P_FDR_BH"].rank(method="min", ascending=True).astype(int)
    data["Contrast"] = str(config["id"])
    selected_n = int(data["Selected_Top30"].sum())
    if selected_n == 0:
        raise ValueError(f"The selection rule retained no nodes for {config['id']}")

    atlas_counts = pd.Series(labels).value_counts().sort_index()
    selected_counts = (
        data.loc[data["Selected_Top30"], "Yeo7_ID"].value_counts().sort_index()
    )
    rows: list[dict[str, object]] = []
    for network_id, network_name in RSN_NAMES.items():
        atlas_n = int(atlas_counts.get(network_id, 0))
        selected = int(selected_counts.get(network_id, 0))
        atlas_fraction = atlas_n / len(labels)
        selected_fraction_network = selected / selected_n
        rows.append(
            {
                "Contrast": str(config["id"]),
                "Yeo7_ID": network_id,
                "Yeo7_Network": network_name,
                "Atlas_Node_Count": atlas_n,
                "Atlas_Fraction": atlas_fraction,
                "Selected_Node_Count": selected,
                "Selected_Fraction": selected_fraction_network,
                "Within_Network_Selected_Fraction": selected / atlas_n,
                "Representation_Ratio": selected_fraction_network / atlas_fraction,
            }
        )
    network = pd.DataFrame(rows)

    summary = {
        "Contrast": str(config["id"]),
        "Source_File": str(source.resolve()),
        "Selection_Fraction": selection_fraction,
        "FDR_Quantile_Threshold": threshold,
        "Selected_Node_Count": selected_n,
        "Selected_FDR_Significant_Count": int(
            data.loc[data["Selected_Top30"], "Significant_FDR_0_05"].sum()
        ),
        "All_Node_FDR_Significant_Count": int(data["Significant_FDR_0_05"].sum()),
        "Min_Selected_P_FDR_BH": float(
            data.loc[data["Selected_Top30"], "P_FDR_BH"].min()
        ),
        "Max_Selected_P_FDR_BH": float(
            data.loc[data["Selected_Top30"], "P_FDR_BH"].max()
        ),
        "Largest_Raw_Count_Network": str(
            network.loc[network["Selected_Node_Count"].idxmax(), "Yeo7_Network"]
        ),
        "Largest_Representation_Ratio_Network": str(
            network.loc[network["Representation_Ratio"].idxmax(), "Yeo7_Network"]
        ),
    }
    return data, network, summary


def radar_figure(
    network: pd.DataFrame,
    contrast_label: str,
    selected_n: int,
    output_stem: Path,
) -> None:
    names = network["Yeo7_Network"].tolist()
    values = network["Selected_Node_Count"].to_numpy(dtype=float)
    angles = np.linspace(0, 2 * np.pi, len(names), endpoint=False)
    closed_angles = np.r_[angles, angles[0]]
    closed_values = np.r_[values, values[0]]

    fig, ax = plt.subplots(figsize=(6.2, 6.2), subplot_kw={"polar": True})
    ax.set_theta_offset(np.pi / 2)
    ax.set_theta_direction(1)
    ax.plot(closed_angles, closed_values, color="#6A3D9A", linewidth=2.2)
    ax.fill(closed_angles, closed_values, color="#6A3D9A", alpha=0.23)
    ax.set_xticks(angles)
    ax.set_xticklabels(names, fontsize=11)
    radial_max = max(100, int(math.ceil(values.max() / 20.0) * 20))
    ax.set_ylim(0, radial_max)
    ax.set_yticks(np.arange(20, radial_max + 1, 20))
    ax.set_yticklabels([str(x) for x in range(20, radial_max + 1, 20)], fontsize=8)
    ax.set_rlabel_position(225)
    fig.suptitle(
        f"{contrast_label}, $\\lambda=0.01$\n"
        f"lowest 30% $p_{{\\mathrm{{FDR}}}}$ ({selected_n} nodes; ties retained)",
        y=0.985,
        fontsize=12,
    )
    ax.grid(color="#BDBDBD", linewidth=0.7)
    fig.tight_layout(rect=(0.02, 0.02, 0.98, 0.90))
    fig.savefig(output_stem.with_suffix(".pdf"), bbox_inches="tight")
    fig.savefig(output_stem.with_suffix(".png"), dpi=300, bbox_inches="tight")
    plt.close(fig)


def representation_figure(networks: pd.DataFrame, output_stem: Path) -> None:
    contrasts = list(networks["Contrast"].drop_duplicates())
    names = [RSN_NAMES[i] for i in RSN_NAMES]
    x = np.arange(len(names))
    width = 0.36
    colors = ("#4C78A8", "#E45756")

    fig, ax = plt.subplots(figsize=(8.4, 4.8))
    for index, contrast in enumerate(contrasts):
        subset = networks[networks["Contrast"] == contrast].sort_values("Yeo7_ID")
        label = contrast.replace("_ABneg", "−").replace("_ABpos", "+").replace("_vs_", " vs. ")
        ax.bar(
            x + (index - 0.5) * width,
            subset["Representation_Ratio"],
            width,
            color=colors[index],
            label=label,
        )
    ax.axhline(1.0, color="#333333", linestyle="--", linewidth=1)
    ax.set_xticks(x)
    ax.set_xticklabels(names)
    ax.set_ylabel("Representation ratio")
    ax.set_title("Yeo-7 representation among the lowest 30% FDR-adjusted p-values")
    ax.legend(frameon=False)
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    fig.savefig(output_stem.with_suffix(".pdf"), bbox_inches="tight")
    fig.savefig(output_stem.with_suffix(".png"), dpi=300, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    args = parse_args()
    if not 0 < args.selection_fraction < 1:
        raise ValueError("--selection-fraction must be strictly between 0 and 1")
    labels = load_mapping(args.mapping_file)

    node_tables: list[pd.DataFrame] = []
    network_tables: list[pd.DataFrame] = []
    summaries: list[dict[str, object]] = []
    for config in CONTRASTS:
        nodes, network, summary = analyse_contrast(
            config, labels, args.selection_fraction
        )
        node_tables.append(nodes)
        network_tables.append(network)
        summaries.append(summary)
        print(
            f"{config['id']}: selected {summary['Selected_Node_Count']}/1000 nodes "
            f"at p_FDR <= {summary['FDR_Quantile_Threshold']:.12g}; "
            f"largest raw count={summary['Largest_Raw_Count_Network']}; "
            "largest representation ratio="
            f"{summary['Largest_Representation_Ratio_Network']}"
        )

    if args.validate_only:
        print("Validation completed; no files were written.")
        return

    output_dir = args.output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    all_nodes = pd.concat(node_tables, ignore_index=True)
    all_networks = pd.concat(network_tables, ignore_index=True)
    summary_table = pd.DataFrame(summaries)

    nodes_path = output_dir / "N145_Yeo7_top30_node_assignments_age_sex_education.csv"
    network_path = output_dir / "N145_Yeo7_top30_network_summary_age_sex_education.csv"
    contrast_path = output_dir / "N145_Yeo7_top30_contrast_summary_age_sex_education.csv"
    all_nodes.to_csv(nodes_path, index=False)
    all_networks.to_csv(network_path, index=False)
    summary_table.to_csv(contrast_path, index=False)

    for config, network, summary in zip(CONTRASTS, network_tables, summaries):
        radar_figure(
            network,
            str(config["label"]),
            int(summary["Selected_Node_Count"]),
            output_dir
            / f"N145_Yeo7_top30_{config['id']}_age_sex_education",
        )
    representation_figure(
        all_networks,
        output_dir / "N145_Yeo7_top30_representation_ratio_age_sex_education",
    )

    print(f"Saved node assignments: {nodes_path}")
    print(f"Saved network summary: {network_path}")
    print(f"Saved contrast summary: {contrast_path}")
    print(f"Saved PDF and PNG figures under: {output_dir}")


if __name__ == "__main__":
    main()
