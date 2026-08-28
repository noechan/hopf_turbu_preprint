%% Prepare data for ML classification with PTIDs and lam_ prefix
clear all; close all;
addpath(genpath('/path/to/ADNI3/code/turbulence_fulldenoising'))

code_path = '/path/to/ADNI3/code/turbulence_fulldenoising/sch1000_N238rev';

%% Load ABeta status table (unchanged)
abeta_table = readtable('/path/to/ADNI3/participants/cross_sectional/ADNI3_N238rev_with_ABETA_Status_CL24.xlsx');

%% Load PTIDs (unchanged)
ptid_path = '/path/to/ADNI3/code/turbulence_fulldenoising/sch1000_N238rev';

load(fullfile(ptid_path, 'PTID_ADNI3_HC_MPRAGE_IRFSPGR_all.mat'))
PTID_HC = cellstr(PTID);

load(fullfile(ptid_path, 'PTID_ADNI3_MCI_MPRAGE_IRFSPGR_all.mat'))
PTID_MCI = cellstr(PTID);

load(fullfile(ptid_path, 'PTID_ADNI3_AD_MPRAGE_IRFSPGR_all.mat'))
PTID_AD = cellstr(PTID);

%% General Aβ mask helper (extended beyond MCI)
% NOTE: If your GROUP labels in the ABeta table use 'CN' instead of 'HC',
% change 'HC' below to 'CN'.
get_abeta_masks = @(group_ids, group_name) deal( ...
    ismember(group_ids, abeta_table.PTID(strcmp(abeta_table.GROUP, group_name) & abeta_table.ABeta_pvc == 1)), ...
    ismember(group_ids, abeta_table.PTID(strcmp(abeta_table.GROUP, group_name) & abeta_table.ABeta_pvc == 0)));

[isAbetaPos_HC,  isAbetaNeg_HC]  = get_abeta_masks(PTID_HC,  'HC');
[isAbetaPos_MCI, ~]              = get_abeta_masks(PTID_MCI, 'MCI');  % Only Aβ+ needed
[isAbetaPos_AD,  ~]              = get_abeta_masks(PTID_AD,  'AD');   % Only Aβ+ needed

%% Load results (unchanged)
data_path = '/path/to/ADNI3/timeseries/outputs_fulldenoising/';

results_HC = load(fullfile(data_path, 'Turbulence_ADNI3_HC_MPRAGE_IRFSPGR_N238rev_sch1000', ...
    'turbu_all_measurements_ADNI3_HC_MPRAGE_IRFSPGR_N238rev_sch1000.mat'));
results_MCI = load(fullfile(data_path, 'Turbulence_ADNI3_MCI_MPRAGE_IRFSPGR_N238rev_sch1000', ...
    'turbu_all_measurements_ADNI3_MCI_MPRAGE_IRFSPGR_N238rev_sch1000.mat'));
results_AD = load(fullfile(data_path, 'Turbulence_ADNI3_AD_MPRAGE_IRFSPGR_N238rev_sch1000', ...
    'turbu_all_measurements_ADNI3_AD_MPRAGE_IRFSPGR_N238rev_sch1000.mat'));

%% Assign data (unchanged)
% HC
Turbu_HC = results_HC.Turbulence_global_sub;
Info_flow_HC = results_HC.TransferLambda_sub;
Info_cascade_HC = results_HC.InformationCascade_sub;
Info_transfer_HC = results_HC.Transfer_sub;
Metastability_HC = results_HC.Meta;
gKoP_HC = results_HC.gKoP;

% MCI
Turbu_MCI = results_MCI.Turbulence_global_sub;
Info_flow_MCI = results_MCI.TransferLambda_sub;
Info_cascade_MCI = results_MCI.InformationCascade_sub;
Info_transfer_MCI = results_MCI.Transfer_sub;
Metastability_MCI = results_MCI.Meta;
gKoP_MCI = results_MCI.gKoP;

% AD
Turbu_AD = results_AD.Turbulence_global_sub;
Info_flow_AD = results_AD.TransferLambda_sub;
Info_cascade_AD = results_AD.InformationCascade_sub;
Info_transfer_AD = results_AD.Transfer_sub;
Metastability_AD = results_AD.Meta;
gKoP_AD = results_AD.gKoP;

%% Split by Aβ status
% HC → Aβ+ and Aβ−
Turbu_HC_ABpos         = Turbu_HC(:, isAbetaPos_HC);
Turbu_HC_ABneg         = Turbu_HC(:, isAbetaNeg_HC);
Info_flow_HC_ABpos     = Info_flow_HC(:, isAbetaPos_HC);
Info_flow_HC_ABneg     = Info_flow_HC(:, isAbetaNeg_HC);
Info_cascade_HC_ABpos  = Info_cascade_HC(:, isAbetaPos_HC);
Info_cascade_HC_ABneg  = Info_cascade_HC(:, isAbetaNeg_HC);
Info_transfer_HC_ABpos = Info_transfer_HC(:, isAbetaPos_HC);
Info_transfer_HC_ABneg = Info_transfer_HC(:, isAbetaNeg_HC);
Metastability_HC_ABpos = Metastability_HC(:, isAbetaPos_HC);
Metastability_HC_ABneg = Metastability_HC(:, isAbetaNeg_HC);
gKoP_HC_ABpos          = gKoP_HC(:, isAbetaPos_HC);
gKoP_HC_ABneg          = gKoP_HC(:, isAbetaNeg_HC);
PTID_HC_ABpos          = PTID_HC(isAbetaPos_HC);
PTID_HC_ABneg          = PTID_HC(isAbetaNeg_HC);

% MCI → only Aβ+
Turbu_MCI_ABpos         = Turbu_MCI(:, isAbetaPos_MCI);
Info_flow_MCI_ABpos     = Info_flow_MCI(:, isAbetaPos_MCI);
Info_cascade_MCI_ABpos  = Info_cascade_MCI(:, isAbetaPos_MCI);
Info_transfer_MCI_ABpos = Info_transfer_MCI(:, isAbetaPos_MCI);
Metastability_MCI_ABpos = Metastability_MCI(:, isAbetaPos_MCI);
gKoP_MCI_ABpos          = gKoP_MCI(:, isAbetaPos_MCI);
PTID_MCI_ABpos          = PTID_MCI(isAbetaPos_MCI);

% AD → only Aβ+
Turbu_AD_ABpos         = Turbu_AD(:, isAbetaPos_AD);
Info_flow_AD_ABpos     = Info_flow_AD(:, isAbetaPos_AD);
Info_cascade_AD_ABpos  = Info_cascade_AD(:, isAbetaPos_AD);
Info_transfer_AD_ABpos = Info_transfer_AD(:, isAbetaPos_AD);
Metastability_AD_ABpos = Metastability_AD(:, isAbetaPos_AD);
gKoP_AD_ABpos          = gKoP_AD(:, isAbetaPos_AD);
PTID_AD_ABpos          = PTID_AD(isAbetaPos_AD);

%% Lambda info (unchanged)
lambda_values = [0.27, 0.24, 0.21, 0.18, 0.15, 0.12, 0.09, 0.06, 0.03, 0.01];
lambda_labels = strcat('lam_', strrep(string(lambda_values), '.', '_'));
lamb_idx = 1:length(lambda_values);

%% Feature assembly (unchanged helpers)
make_features = @(T, IF, IT, idx) [T(idx,:)', IF(idx,:)', IT(idx,:)'];
add_1D = @(casc, meta, gkop) [casc(:), meta(:), gkop(:)];

% HC Aβ+
X_HC_ABpos = [make_features(Turbu_HC_ABpos, Info_flow_HC_ABpos, Info_transfer_HC_ABpos, lamb_idx), ...
              add_1D(Info_cascade_HC_ABpos, Metastability_HC_ABpos, gKoP_HC_ABpos)];

% HC Aβ−
X_HC_ABneg = [make_features(Turbu_HC_ABneg, Info_flow_HC_ABneg, Info_transfer_HC_ABneg, lamb_idx), ...
              add_1D(Info_cascade_HC_ABneg, Metastability_HC_ABneg, gKoP_HC_ABneg)];

% MCI Aβ+
X_MCI_ABpos = [make_features(Turbu_MCI_ABpos, Info_flow_MCI_ABpos, Info_transfer_MCI_ABpos, lamb_idx), ...
               add_1D(Info_cascade_MCI_ABpos, Metastability_MCI_ABpos, gKoP_MCI_ABpos)];

% AD Aβ+
X_AD_ABpos = [make_features(Turbu_AD_ABpos, Info_flow_AD_ABpos, Info_transfer_AD_ABpos, lamb_idx), ...
              add_1D(Info_cascade_AD_ABpos, Metastability_AD_ABpos, gKoP_AD_ABpos)];

%% PTIDs per requested groups (order as requested)
PTID_all = [PTID_HC_ABpos(:); PTID_HC_ABneg(:); PTID_MCI_ABpos(:); PTID_AD_ABpos(:)];
Group_labels = [ ...
    repmat("HC_ABpos", size(X_HC_ABpos,1), 1); ...
    repmat("HC_ABneg", size(X_HC_ABneg,1), 1); ...
    repmat("MCI_ABpos", size(X_MCI_ABpos,1), 1); ...
    repmat("AD_ABpos",  size(X_AD_ABpos,1), 1)];

%% Combine all (order matches PTID_all & Group_labels)
X_all = [X_HC_ABpos; X_HC_ABneg; X_MCI_ABpos; X_AD_ABpos];

%% Feature names (unchanged)
turbu_names = strcat("Turbu_", lambda_labels(lamb_idx));
infoflow_names = strcat("InfoFlow_", lambda_labels(lamb_idx));
infotrans_names = strcat("InfoTransfer_", lambda_labels(lamb_idx));
other_names = ["InfoCascade", "Metastability", "gKoP"];
var_names = [turbu_names, infoflow_names, infotrans_names, other_names];

%% Convert to table and save (unchanged except groups now reflect Aβ split)
T_features = array2table(X_all, 'VariableNames', matlab.lang.makeValidName(var_names));
T_features.PTID = string(PTID_all);
T_features.Group = categorical(Group_labels);

% Reorder columns: PTID, Group, then features
T_features = movevars(T_features, {'PTID','Group'}, 'Before', 1);

save(fullfile(code_path, 'ML_Input_ADNI3_4STAGINGBYABETA__N238rev_withPTID_sch1000.mat'), 'T_features');
writetable(T_features, fullfile(code_path, 'ML_Input_ADNI3_4STAGINGBYABETA_N238rev_withPTID_sch1000.csv'));
writetable(T_features, fullfile(code_path, 'ML_Input_ADNI3_4STAGINGBYABETA_N238rev_withPTID_sch1000.xlsx'), 'FileType', 'spreadsheet');

disp('✔ Table saved with PTIDs and lam_ naming convention for HC Aβ+, HC Aβ−, MCI Aβ+, AD Aβ+.');
