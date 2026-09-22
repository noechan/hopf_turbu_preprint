import json
import multiprocessing
import os
import sys
import tempfile
from pathlib import Path

# Allow this file to be launched directly by VS Code or with its absolute path.
path_repo = Path(__file__).resolve().parents[2]
if str(path_repo) not in sys.path:
    sys.path.insert(0, str(path_repo))

# Use writable, non-interactive plotting settings when launched from VS Code.
matplotlib_cache = Path(tempfile.gettempdir()) / "adni3_ml_matplotlib_cache"
matplotlib_cache.mkdir(parents=True, exist_ok=True)
os.environ.setdefault("MPLCONFIGDIR", str(matplotlib_cache))
os.environ.setdefault("MPLBACKEND", "Agg")

from src.data_loading.load_from_excel import get_x_arr_for_scikit_from_excel
from src.pipelines.train_eval_pipelines import (
    train_eval_gridsearch_loocv_with_outer_n_loop,
)

# Determine the number of CPU cores to use
n_jobs = multiprocessing.cpu_count() - 1  # Leave one core free

data_files = {
    "all_features_turbu_combat": "ML_Input_ADNI3_4STAGINGBYABETA_ComBat_N145_withPTID_sch1000.xlsx",
}

classifier = "PolySVM"
classifications = ["HCneg_vs_HCpos", "HCneg_vs_MCIpos","HCneg_vs_ADpos","MCIpos_vs_ADpos"]
group_labels = {"HCneg": "HC_ABneg", "HCpos": "HC_ABpos","MCIpos": "MCI_ABpos", "ADpos": "AD_ABpos"}

if __name__ == "__main__":
    # Main execution
    excel_folder = path_repo / "Data" / "turbu_hopf"
    param_folder = path_repo / "Parameters"
    for classification in classifications:
        for data_type in data_files.keys():
            print(f"Training {classifier} with parameters_{classifier}.json")

            results_folder = (
                    path_repo / "Results" / "final_3d_gs_classification_turbu_sch1000" /
                    classifier / classification / data_type
            )
            results_folder.mkdir(exist_ok=True, parents=True)

            data_file = excel_folder / data_files[data_type]

            # Load the parameters file
            param_file = path_repo / "Parameters" / f"parameters_{classifier}.json"
            with open(param_file, "r") as file:
                parameters = json.load(file)
            parameters["N_JOBS"] = n_jobs

            if classification == "HCneg_vs_HCpos":
                parameters["GROUPS"] = {"HC_ABneg": 0, "HC_ABpos": 1}
            elif classification == "HCneg_vs_MCIpos":
                parameters["GROUPS"] = {"HC_ABneg": 0, "MCI_ABpos": 1}
            elif classification == "HCneg_vs_ADpos":
                parameters["GROUPS"] = {"HC_ABneg": 0, "AD_ABpos": 1}
            elif classification == "MCIpos_vs_ADpos":
                parameters["GROUPS"] = {"MCI_ABpos": 0, "AD_ABpos": 1}
            else:
                pass

            print(f"Training K-fold CV with {classifier} Grid Search and LOOCV {data_file} Data")

            # We load the data from the file
            x, y, feature_columns = get_x_arr_for_scikit_from_excel(
                data_file,
                parameters["GROUPS"],
                parameters["ID_KEY"],
                parameters["GROUP_KEY"],
            )

            # And we run the train_eval_knn gridsearch which performs gridsearch,
            # trains on best parameters, and evaluates the model and reports the
            # performances.
            train_eval_gridsearch_loocv_with_outer_n_loop(
                x, y, feature_columns, parameters, results_folder,
                n=40
            )
