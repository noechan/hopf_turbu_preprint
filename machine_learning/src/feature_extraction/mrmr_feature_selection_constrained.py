import numpy as np
from sklearn.feature_selection import mutual_info_regression, mutual_info_classif


def mrmr_feature_selection_deterministic(X, y, k=15, random_state=0):
    """
    Minimum redundancy maximum relevance (mRMR) feature selection (greedy).
    """

    X = np.asarray(X)
    y = np.asarray(y).ravel()

    n_samples, n_features = X.shape

    if n_features <= k:
        return list(range(n_features))

    # Relevance: MI(feature, y)
    relevance = mutual_info_classif(X, y, random_state=random_state)

    selected = []
    not_selected = list(range(n_features))

    first = int(np.argmax(relevance))
    selected.append(first)
    not_selected.remove(first)

    for _ in range(k - 1):
        if not not_selected:
            break

        scores = []
        for j in not_selected:
            redundancy = 0.0
            if selected:
                X_candidate = X[:, j].reshape(-1, 1)
                redundancy = float(np.mean([
                    mutual_info_regression(
                        X_candidate, X[:, sel].ravel(), random_state=random_state
                    )[0]
                    for sel in selected
                ]))

            score = relevance[j] if redundancy < 1e-10 else relevance[j] / redundancy
            scores.append(score)

        scores = np.asarray(scores, dtype=float)

        # Minimal deterministic selection
        next_feature = int(not_selected[int(np.argmax(scores))])

        selected.append(next_feature)
        not_selected.remove(next_feature)

    return selected


def constrained_mrmr_feature_selection(
    X,
    y,
    feature_columns,
    k,
    always_include=None,
    exclude=None,
    random_state=0,
):
    """
    Select exactly k final features with the following constraints:
    - features in always_include are always present
    - features in exclude are never used
    - mRMR is applied only to the remaining eligible features

    Returns
    -------
    selected_indices : list[int]
        Indices referred to the ORIGINAL feature_columns / X columns.
    """

    X = np.asarray(X)
    y = np.asarray(y).ravel()

    always_include = always_include or []
    exclude = exclude or []

    feature_to_idx = {f: i for i, f in enumerate(feature_columns)}

    missing_always = [f for f in always_include if f not in feature_to_idx]
    missing_exclude = [f for f in exclude if f not in feature_to_idx]

    if missing_always:
        raise ValueError(f"Features in ALWAYS_INCLUDE_FEATURES not found: {missing_always}")
    if missing_exclude:
        raise ValueError(f"Features in EXCLUDED_FEATURES not found: {missing_exclude}")

    excluded_idx = {feature_to_idx[f] for f in exclude}
    fixed_idx = [feature_to_idx[f] for f in always_include if feature_to_idx[f] not in excluded_idx]

    if k < len(fixed_idx):
        raise ValueError(
            f"NF={k} is smaller than the number of always-included features "
            f"({len(fixed_idx)})."
        )

    eligible_idx = [
        i for i, f in enumerate(feature_columns)
        if i not in excluded_idx and i not in fixed_idx
    ]

    n_to_select_with_mrmr = k - len(fixed_idx)

    if n_to_select_with_mrmr <= 0:
        return fixed_idx

    if n_to_select_with_mrmr >= len(eligible_idx):
        return fixed_idx + eligible_idx

    X_eligible = X[:, eligible_idx]
    selected_rel_idx = mrmr_feature_selection_deterministic(
        X_eligible,
        y,
        k=n_to_select_with_mrmr,
        random_state=random_state,
    )

    selected_abs_idx = [eligible_idx[i] for i in selected_rel_idx]
    return fixed_idx + selected_abs_idx