import json
import joblib

import numpy as np

from pathlib import Path
from sklearn.model_selection import StratifiedKFold, train_test_split

from .model_evaluation import evaluate_model
from .model_training_constrained import (
    pipeline_gridsearch_3d_with_loocv,
    classifiers
)
from ..pipelines.shap_feature_importance import pipeline_feature_importance
from ..feature_extraction.mrmr_feature_selection_constrained import (
    constrained_mrmr_feature_selection,
)
from ..utils.summary_pipelines import (
    save_results_summary,
    save_statistics_gridsearch_loocv_pipeline_with_outer_kfold_loop,
)


def train_eval_gridsearch_3d_with_loocv(
    x_train, y_train, x_test, y_test, feature_columns, config_params, outputs_folder
):
    # Perform the gridsearch to obtain the best parameters
    print("Starting 3D GridSearch")

    grid_search_params = {
        "USE_MRMR_FEATURE_SELECTION": config_params["USE_MRMR_FEATURE_SELECTION"],
        "MODEL_TYPE": config_params["MODEL_TYPE"],
        "MODEL_KWARGS": config_params["MODEL_KWARGS"],
        "STANDARDIZE_DATA": config_params["STANDARDIZE_DATA"],
        "RANDOM_STATE": config_params["RANDOM_STATE"],
        "ALWAYS_INCLUDE_FEATURES": config_params.get("ALWAYS_INCLUDE_FEATURES", []),
        "EXCLUDED_FEATURES": config_params.get("EXCLUDED_FEATURES", []),
        **config_params["GRID_SEARCH"]
    }

    best_hyperparams, best_score, _ = pipeline_gridsearch_3d_with_loocv(
        x_train, y_train, feature_columns, config_params=grid_search_params
    )
    grid_search_results = {**best_hyperparams, "Score": best_score}
    print("GridSearch Finished")

    # Before model evaluation, perform constrained feature selection (if applies)
    if config_params["USE_MRMR_FEATURE_SELECTION"]:
        selected_indices = constrained_mrmr_feature_selection(
            x_train,
            y_train,
            feature_columns,
            k=best_hyperparams["NF"],
            always_include=config_params.get("ALWAYS_INCLUDE_FEATURES", []),
            exclude=config_params.get("EXCLUDED_FEATURES", []),
            random_state=config_params.get("RANDOM_STATE", 0),
        )

        best_selected_features = [feature_columns[i] for i in selected_indices]

        # Use selected indices for both training and testing
        x_train = x_train[:, selected_indices]
        x_test = x_test[:, selected_indices]
    else:
        # If not using mRMR, unless "NF" in best_hyperparams, take all:
        if "NF" not in best_hyperparams:
            x_train = x_train[:, :]
            x_test = x_test[:, :]
            best_selected_features = feature_columns[:]
        else:  # If NF has been explored, take first NF
            x_train = x_train[:, :best_hyperparams["NF"]]
            x_test = x_test[:, :best_hyperparams["NF"]]
            best_selected_features = feature_columns[:best_hyperparams["NF"]]

    # Train model on entire data and evaluate performance (generates ROCs)
    print("Evaluating Optimized Model")

    # Generate kwargs for the final model
    args_model = {key: value for key, value in best_hyperparams.items() if key != "NF"}

    # Add additional kwargs
    if config_params["MODEL_KWARGS"] is not None:
        args_model = {**args_model, **config_params["MODEL_KWARGS"]}

    # Initialize model
    final_model = classifiers[grid_search_params["MODEL_TYPE"]](**args_model)

    final_model, *results = evaluate_model(
        final_model,
        x_train,
        y_train,
        x_test,
        y_test,
        outputs_folder=outputs_folder,
        random_state=config_params["RANDOM_STATE"],
        group_labels=list(config_params["GROUPS"].keys()),
        standardize_data=config_params["STANDARDIZE_DATA"]
    )

    # From the reported data, obtain relevant information except ROC data
    scores_train, scores_test, _, _, test_acc_pvalue = results

    # Perform feature importance analyses and store in folder
    print("Computing SHAP Feature Importance")
    importance_dict = pipeline_feature_importance(
        final_model,
        x_test,
        y_test,
        best_selected_features,
        outputs_folder,
        random_state=config_params["RANDOM_STATE"],
        n_jobs=config_params["N_JOBS"],
    )

    # Store summary information in a .json file with all the results for later analyses
    save_results_summary(
        scores_train,
        scores_test,
        test_acc_pvalue,
        grid_search_results,
        importance_dict,
        best_selected_features,
        outputs_folder,
    )

    # Also store the model with joblib
    joblib.dump(final_model, outputs_folder / "fit_model.joblib")

    # Finally, for enhanced reproducibility, we also store the parameters in a .json
    with open(outputs_folder / "parameters.json", "w") as fp:
        json.dump(config_params, fp, indent=4)

    print("Pipeline Finished\n")


def train_eval_gridsearch_loocv_with_outer_n_loop(
    x, y, feature_columns, config_params, outputs_folder, n=20
):
    """
    This function runs N iterations of the 3D Gridsearch with LOOCV approach,
    generalized to several classifiers and allowing for flexibility in the
    hyperparameter choice.

    The main goal is to obtain robust estimations of the machine learning algorithm
    methodology against single train/test splits.

    Technical consideration: We will make N calls to the
    train_eval_gridsearch_3d_with_loocv function, each time storing results in a
    different folder.

    Args:
        x: input-array of shape (n_samples, n_features)
        y: labels-array of shape (n_samples, )
        feature_columns: list of str with names of the features (columns in x)
        config_params: dictionary with parameters
        outputs_folder: Path to the folder where results will be output
        n: (int) number of iterations in the loop calling the pipelines
    """
    # In this case, we will not perform a KFold, but a single train/test split for
    # each of the iterations with a different random state.

    np.random.seed(config_params["RANDOM_STATE"])

    ################## HERE'S THE TRAINING FUNCTION #############################

    train_eval_function = train_eval_gridsearch_3d_with_loocv

    ##############################################################################
    for i in range(n):
        print(f"\n -------- Running GridSearch Pipeline for Iteration {i} -------- \n")

        # For each attempt we fix a different random state
        config_params["RANDOM_STATE"] = np.random.randint(100)

        x_train, x_test, y_train, y_test = train_test_split(
            x,
            y,
            test_size=config_params["TEST_SIZE"],
            random_state=config_params["RANDOM_STATE"]
        )

        # For each fold, we generate a new outputs_folder path to store fold-results
        outputs_folder_fold = Path(outputs_folder) / f"fold_{i}"
        outputs_folder_fold.mkdir(exist_ok=True)

        # Call the function to perform the gridsearch and optimization for the fold
        train_eval_function(
            x_train,
            y_train,
            x_test,
            y_test,
            feature_columns,
            config_params,
            outputs_folder_fold,
        )

    print("Finished Running GridSearch pipeline for all folds")

    # Finally, accumulate the results obtained, computing statistics over results
    save_statistics_gridsearch_loocv_pipeline_with_outer_kfold_loop(
        outputs_folder, n, config_params
    )