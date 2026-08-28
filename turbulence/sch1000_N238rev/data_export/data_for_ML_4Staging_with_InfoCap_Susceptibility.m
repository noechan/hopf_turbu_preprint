%% Export post-ComBat N145 dynamics plus Info_Cap and Susceptibility for ML
% Uses the same cohort and dynamics validation as data_for_ML_4Staging.m,
% then aligns the two participant-level Hopf measures by PTID. The output
% preserves the original combined ML layout:
%
%   PTID, Group, Info_Cap, Susceptibility, 32 dynamics features
%
% The script never joins by row position and does not alter either source.

clear; close all; clc;

sch1000_root = fileparts(fileparts(mfilename('fullpath')));
addpath(sch1000_root, '-begin');
P = setup_sch1000_paths();

%% Input and output names
harmonized_file = fullfile(P.harmonization_results, ...
    'Turbu_ComBat_ADNI3_allfeatures_N145.xlsx');
hopf_filename = 'Turbu_ComBat_ADNI3_2feat_IC_SUSC_N145.xlsx';

% Prefer an explicitly supplied source. Otherwise look beside the MATLAB
% exports and then in the local ML repository's data folder.
hopf_file = string(getenv('ADNI3_HOPF_ML_INPUT'));
if strlength(hopf_file) == 0
    code_root = fileparts(P.repo_root);
    hopf_candidates = [ ...
        string(fullfile(P.data_export, hopf_filename)), ...
        string(fullfile(code_root, 'ADNI3_machine_learning-main-det-int', ...
            'Data', 'turbu_hopf', hopf_filename))];
    candidate_index = find(isfile(hopf_candidates), 1, 'first');
    if isempty(candidate_index)
        error(['Could not locate the participant-level Hopf workbook.\n' ...
            'Set ADNI3_HOPF_ML_INPUT to its full path or place %s in:\n%s'], ...
            hopf_filename, P.data_export);
    end
    hopf_file = hopf_candidates(candidate_index);
end

base_name = ['ML_Input_ADNI3_4STAGINGBYABETA_ComBat_N145_' ...
    'with_infocap_suscep_sch1000'];

if ~isfile(harmonized_file)
    error(['Missing recomputed harmonized table: %s\n' ...
        'Run harmonization_allfeat/run_harmonization.py first.'], ...
        harmonized_file);
end
if ~isfile(hopf_file)
    error('Missing participant-level Hopf table: %s', hopf_file);
end
if ~isfolder(P.data_export)
    mkdir(P.data_export);
end

%% Expected original combined ML structure
group_labels = ["HC_ABneg", "HC_ABpos", "MCI_ABpos", "AD_ABpos"];
expected_counts = [51, 37, 31, 26];

lambda_values = [0.27, 0.24, 0.21, 0.18, 0.15, 0.12, 0.09, 0.06, 0.03, 0.01];
lambda_text = strrep(compose('%.2f', lambda_values), '.', '_');
turbu_names = "Turbu_lam_" + lambda_text;

% TransferLambda_sub(1,:) is an initialization row and is not an estimated
% information-flow feature.
infoflow_names = "InfoFlow_lam_" + lambda_text(2:end);
infotransfer_names = "InfoTransfer_lam_" + lambda_text;
dynamics_names = [turbu_names, infoflow_names, infotransfer_names, ...
    "InfoCascade", "Metastability", "gKoP"];
dynamics_variables = ["PTID", "Group", dynamics_names];
hopf_variables = ["PTID", "Group", "Info_Cap", "Susceptibility"];
combined_variables = ["PTID", "Group", "Info_Cap", ...
    "Susceptibility", dynamics_names];

%% Read both N145 inputs
T_dynamics = readtable(harmonized_file, 'VariableNamingRule', 'preserve');
T_hopf = readtable(hopf_file, 'VariableNamingRule', 'preserve');

actual_dynamics = string(T_dynamics.Properties.VariableNames);
missing_dynamics = setdiff(dynamics_variables, actual_dynamics, 'stable');
unexpected_dynamics = setdiff(actual_dynamics, dynamics_variables, 'stable');
if ~isempty(missing_dynamics) || ~isempty(unexpected_dynamics)
    error(['The harmonized dynamics table has an unexpected structure.\n' ...
        'Missing variables: %s\nUnexpected variables: %s'], ...
        strjoin(missing_dynamics, ', '), ...
        strjoin(unexpected_dynamics, ', '));
end

actual_hopf = string(T_hopf.Properties.VariableNames);
missing_hopf = setdiff(hopf_variables, actual_hopf, 'stable');
if ~isempty(missing_hopf)
    error('The Hopf table is missing variables: %s', ...
        strjoin(missing_hopf, ', '));
end

T_dynamics = T_dynamics(:, cellstr(dynamics_variables));
T_hopf = T_hopf(:, cellstr(hopf_variables));
T_dynamics.PTID = strip(string(T_dynamics.PTID));
T_hopf.PTID = strip(string(T_hopf.PTID));
dynamics_group = strip(string(T_dynamics.Group));
hopf_group = strip(string(T_hopf.Group));

%% Validate both canonical N145 cohorts
expected_n = sum(expected_counts);
if height(T_dynamics) ~= expected_n || height(T_hopf) ~= expected_n
    error(['Expected 145 rows in both inputs; found dynamics=%d and ' ...
        'Hopf=%d.'], height(T_dynamics), height(T_hopf));
end
if any(ismissing(T_dynamics.PTID)) || any(T_dynamics.PTID == "") || ...
        any(ismissing(T_hopf.PTID)) || any(T_hopf.PTID == "")
    error('One of the input tables contains a missing PTID.');
end
if numel(unique(T_dynamics.PTID)) ~= expected_n || ...
        numel(unique(T_hopf.PTID)) ~= expected_n
    error('One of the input tables contains duplicate PTIDs.');
end

unknown_dynamics_groups = setdiff(unique(dynamics_group), group_labels);
unknown_hopf_groups = setdiff(unique(hopf_group), group_labels);
if ~isempty(unknown_dynamics_groups) || ~isempty(unknown_hopf_groups)
    error(['Unexpected group labels. Dynamics: %s; Hopf: %s'], ...
        strjoin(unknown_dynamics_groups, ', '), ...
        strjoin(unknown_hopf_groups, ', '));
end

observed_counts = arrayfun(@(label) ...
    nnz(dynamics_group == label), group_labels);
hopf_counts = arrayfun(@(label) nnz(hopf_group == label), group_labels);
if ~isequal(observed_counts, expected_counts) || ...
        ~isequal(hopf_counts, expected_counts)
    error(['Unexpected group counts. Dynamics=[%s], Hopf=[%s], ' ...
        'expected=[%s].'], num2str(observed_counts), ...
        num2str(hopf_counts), num2str(expected_counts));
end

X_dynamics = T_dynamics{:, cellstr(dynamics_names)};
X_hopf = T_hopf{:, {'Info_Cap', 'Susceptibility'}};
if ~isnumeric(X_dynamics) || ~isequal(size(X_dynamics), [expected_n, 32])
    error('Expected a numeric 145-by-32 dynamics matrix.');
end
if ~isnumeric(X_hopf) || ~isequal(size(X_hopf), [expected_n, 2])
    error('Expected a numeric 145-by-2 Hopf matrix.');
end
if any(~isfinite(X_dynamics), 'all') || any(~isfinite(X_hopf), 'all')
    error('One of the input feature matrices contains NaN or Inf values.');
end

%% Align Hopf measures to the harmonized table by PTID
[matched, hopf_index] = ismember(T_dynamics.PTID, T_hopf.PTID);
if ~all(matched) || numel(unique(hopf_index)) ~= expected_n
    missing_ptids = T_dynamics.PTID(~matched);
    error('The Hopf table does not match all dynamics PTIDs: %s', ...
        strjoin(missing_ptids, ', '));
end
if any(dynamics_group ~= hopf_group(hopf_index))
    mismatch = T_dynamics.PTID(dynamics_group ~= hopf_group(hopf_index));
    error('Group labels disagree between inputs for PTIDs: %s', ...
        strjoin(mismatch, ', '));
end

info_cap = T_hopf.Info_Cap(hopf_index);
susceptibility = T_hopf.Susceptibility(hopf_index);
T_features = addvars(T_dynamics, info_cap, susceptibility, ...
    'After', 'Group', ...
    'NewVariableNames', {'Info_Cap', 'Susceptibility'});
T_features = T_features(:, cellstr(combined_variables));
T_features.Group = categorical(dynamics_group, group_labels, group_labels);

X_combined = T_features{:, 3:end};
if ~isequal(size(X_combined), [expected_n, 34]) || ...
        any(~isfinite(X_combined), 'all')
    error('The final combined feature matrix must be finite and 145-by-34.');
end

%% Export ML-ready copies
mat_file = fullfile(P.data_export, [base_name '.mat']);
csv_file = fullfile(P.data_export, [base_name '.csv']);
xlsx_file = fullfile(P.data_export, [base_name '.xlsx']);

save(mat_file, 'T_features');
writetable(T_features, csv_file);
writetable(T_features, xlsx_file, 'FileType', 'spreadsheet');

fprintf('\nCombined post-harmonization ML export completed.\n');
fprintf('Dynamics source: %s\n', harmonized_file);
fprintf('Hopf source: %s\n', hopf_file);
fprintf('Cohort: %d participants; features: %d.\n', ...
    height(T_features), width(T_features) - 2);
for g = 1:numel(group_labels)
    fprintf('  %s: %d\n', group_labels(g), observed_counts(g));
end
fprintf('XLSX to copy into the ML data folder:\n%s\n', xlsx_file);
