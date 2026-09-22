"""Validated accumulation of repeated-classification results into Excel."""

from __future__ import annotations

from pathlib import Path
from typing import Mapping, Sequence

import pandas as pd

from .summary_functions import collect_selected_metrics


DEFAULT_METRICS = (
    "Balanced Accuracy (Test)",
    "AUC (Test)",
    "Sensitivity (Test)",
    "Specificity (Test)",
    "F1 (Test)",
    "Best NF",
)


def _directory_names(path: Path) -> set[str]:
    return {item.name for item in path.iterdir() if item.is_dir()}


def validate_result_tree(
    results_dir: Path,
    feature_labels: Mapping[str, str],
    group_labels: Mapping[str, str],
    required_classifiers: Sequence[str],
    expected_folds: int,
    results_file_name: str,
) -> list[str]:
    """Validate the complete classifier/comparison/feature/fold hierarchy."""
    if not results_dir.is_dir():
        raise FileNotFoundError(f"Results directory does not exist: {results_dir}")
    if expected_folds < 1:
        raise ValueError("expected_folds must be a positive integer.")

    classifiers = sorted(_directory_names(results_dir))
    missing_classifiers = sorted(set(required_classifiers) - set(classifiers))
    if missing_classifiers:
        raise ValueError(
            "Missing required classifier result folders: "
            + ", ".join(missing_classifiers)
        )

    expected_comparisons = set(group_labels)
    expected_features = set(feature_labels)
    expected_fold_names = {f"fold_{index}" for index in range(expected_folds)}
    errors: list[str] = []

    for classifier in classifiers:
        classifier_dir = results_dir / classifier
        comparisons = _directory_names(classifier_dir)
        if comparisons != expected_comparisons:
            missing = sorted(expected_comparisons - comparisons)
            unexpected = sorted(comparisons - expected_comparisons)
            if missing:
                errors.append(f"{classifier}: missing comparisons {missing}")
            if unexpected:
                errors.append(f"{classifier}: unexpected comparisons {unexpected}")
            continue

        for comparison in group_labels:
            comparison_dir = classifier_dir / comparison
            feature_sets = _directory_names(comparison_dir)
            if feature_sets != expected_features:
                missing = sorted(expected_features - feature_sets)
                unexpected = sorted(feature_sets - expected_features)
                if missing:
                    errors.append(
                        f"{classifier}/{comparison}: missing feature sets {missing}"
                    )
                if unexpected:
                    errors.append(
                        f"{classifier}/{comparison}: unexpected feature sets "
                        f"{unexpected}"
                    )
                continue

            for feature_set in feature_labels:
                leaf = comparison_dir / feature_set
                folds = {
                    item.name
                    for item in leaf.iterdir()
                    if item.is_dir() and item.name.startswith("fold_")
                }
                missing_folds = sorted(
                    expected_fold_names - folds,
                    key=lambda value: int(value.removeprefix("fold_")),
                )
                unexpected_folds = sorted(folds - expected_fold_names)
                if missing_folds:
                    errors.append(
                        f"{classifier}/{comparison}/{feature_set}: "
                        f"missing folds {missing_folds}"
                    )
                if unexpected_folds:
                    errors.append(
                        f"{classifier}/{comparison}/{feature_set}: "
                        f"unexpected folds {unexpected_folds}"
                    )
                missing_fold_summaries = [
                    fold
                    for fold in sorted(
                        folds & expected_fold_names,
                        key=lambda value: int(value.removeprefix("fold_")),
                    )
                    if not (leaf / fold / "summary_results.json").is_file()
                ]
                if missing_fold_summaries:
                    errors.append(
                        f"{classifier}/{comparison}/{feature_set}: folds without "
                        f"summary_results.json {missing_fold_summaries}"
                    )
                if not (leaf / results_file_name).is_file():
                    errors.append(
                        f"{classifier}/{comparison}/{feature_set}: missing "
                        f"{results_file_name}"
                    )

    if errors:
        raise ValueError("Result-tree validation failed:\n- " + "\n- ".join(errors))
    return classifiers


def prepare_accumulated_results(
    results_dir: Path,
    feature_labels: Mapping[str, str],
    group_labels: Mapping[str, str],
    required_classifiers: Sequence[str] = ("LinearSVM", "PolySVM"),
    expected_folds: int = 40,
    metrics: Sequence[str] = DEFAULT_METRICS,
    results_file_name: str = "stats_results_folds.json",
) -> pd.DataFrame:
    """Validate and collect results while retaining numeric and display columns."""
    results_dir = Path(results_dir).resolve()
    classifiers = validate_result_tree(
        results_dir,
        feature_labels,
        group_labels,
        required_classifiers,
        expected_folds,
        results_file_name,
    )
    results = collect_selected_metrics(results_dir, list(metrics), results_file_name)
    if results.empty:
        raise ValueError(f"No accumulated results were found below {results_dir}.")

    expected_rows = len(classifiers) * len(group_labels) * len(feature_labels)
    if len(results) != expected_rows:
        raise ValueError(
            f"Collected {len(results)} result rows; expected {expected_rows} from "
            "the validated result hierarchy."
        )

    key_columns = ["Classifier", "Group Comparison", "Feature Set"]
    if results.duplicated(key_columns).any():
        duplicate_rows = results.duplicated(key_columns, keep=False)
        duplicates = results.loc[duplicate_rows, key_columns]
        raise ValueError(
            "Duplicate classifier/comparison/feature rows were collected:\n"
            + duplicates.to_string(index=False)
        )

    numeric_columns = [
        f"{metric} {suffix}" for metric in metrics for suffix in ("Mean", "Std")
    ]
    missing_numeric = [column for column in numeric_columns if column not in results]
    if missing_numeric:
        raise ValueError(f"Missing accumulated metric columns: {missing_numeric}")
    if results[numeric_columns].isna().any().any():
        bad_columns = results[numeric_columns].columns[
            results[numeric_columns].isna().any()
        ].tolist()
        raise ValueError(
            "One or more stored metric strings could not be parsed: "
            + ", ".join(bad_columns)
        )

    results.insert(
        results.columns.get_loc("Group Comparison"),
        "Group Comparison ID",
        results["Group Comparison"],
    )
    results.insert(
        results.columns.get_loc("Feature Set"),
        "Feature Set ID",
        results["Feature Set"],
    )
    results["Group Comparison"] = results["Group Comparison"].map(group_labels)
    results["Feature Set"] = results["Feature Set"].map(feature_labels)
    if results[["Group Comparison", "Feature Set"]].isna().any().any():
        raise ValueError("A collected group comparison or feature set lacks a label.")

    for metric in metrics:
        mean_column = f"{metric} Mean"
        std_column = f"{metric} Std"
        display_column = f"{metric} Mean ± SD"
        results[display_column] = [
            f"{mean:.3f} ± {standard_deviation:.3f}"
            for mean, standard_deviation in zip(
                results[mean_column], results[std_column]
            )
        ]

    identity_columns = [
        "Classifier",
        "Group Comparison ID",
        "Group Comparison",
        "Feature Set ID",
        "Feature Set",
    ]
    metric_columns = [
        column
        for metric in metrics
        for column in (
            f"{metric} Mean",
            f"{metric} Std",
            f"{metric} Mean ± SD",
        )
    ]
    remaining_columns = [
        column
        for column in results.columns
        if column not in identity_columns and column not in metric_columns
    ]
    results = results[identity_columns + metric_columns + remaining_columns]
    return results.sort_values(
        ["Group Comparison ID", "Feature Set ID", "Classifier"]
    ).reset_index(drop=True)


def export_group_comparison_workbook(
    results: pd.DataFrame,
    output_path: Path,
    group_labels: Mapping[str, str],
) -> Path:
    """Write one non-empty worksheet per validated group comparison."""
    output_path = Path(output_path).resolve()
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with pd.ExcelWriter(output_path, engine="xlsxwriter") as writer:
        workbook = writer.book
        header_format = workbook.add_format(
            {
                "bold": True,
                "font_color": "#FFFFFF",
                "bg_color": "#1F4E78",
                "border": 0,
                "align": "center",
                "valign": "vcenter",
                "text_wrap": True,
            }
        )
        numeric_format = workbook.add_format({"num_format": "0.000"})

        for group_id, sheet_name in group_labels.items():
            group_results = results.loc[
                results["Group Comparison ID"] == group_id
            ].copy()
            if group_results.empty:
                raise ValueError(
                    f"Refusing to write an empty worksheet for {group_id}."
                )
            group_results.to_excel(writer, sheet_name=sheet_name, index=False)
            worksheet = writer.sheets[sheet_name]
            worksheet.hide_gridlines(2)
            worksheet.freeze_panes(1, 0)
            worksheet.autofilter(
                0, 0, len(group_results), len(group_results.columns) - 1
            )
            worksheet.set_row(0, 36)
            for column_index, column_name in enumerate(group_results.columns):
                worksheet.write(0, column_index, column_name, header_format)

            for column_index, column_name in enumerate(group_results.columns):
                values = group_results[column_name].astype(str)
                width = min(
                    max(len(column_name), int(values.str.len().max())) + 2,
                    34,
                )
                if pd.api.types.is_numeric_dtype(group_results[column_name]):
                    worksheet.set_column(
                        column_index, column_index, max(width, 12), numeric_format
                    )
                else:
                    worksheet.set_column(column_index, column_index, max(width, 12))
    return output_path


def accumulate_results(
    results_dir: Path,
    feature_labels: Mapping[str, str],
    group_labels: Mapping[str, str],
    output_file_name: str = "group_comparison_results.xlsx",
    required_classifiers: Sequence[str] = ("LinearSVM", "PolySVM"),
    expected_folds: int = 40,
    check_only: bool = False,
) -> tuple[pd.DataFrame, Path | None]:
    """Prepare results and optionally write the group-comparison workbook."""
    results = prepare_accumulated_results(
        results_dir=results_dir,
        feature_labels=feature_labels,
        group_labels=group_labels,
        required_classifiers=required_classifiers,
        expected_folds=expected_folds,
    )
    if check_only:
        return results, None
    output_path = export_group_comparison_workbook(
        results,
        Path(results_dir) / output_file_name,
        group_labels,
    )
    return results, output_path
