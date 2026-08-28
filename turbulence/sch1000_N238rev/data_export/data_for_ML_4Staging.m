%% Export post-ComBat N145 data for the original ML pipeline
% Reads the canonical all-feature table produced by
% harmonization_allfeat/run_harmonization.py. The exported table preserves
% the original ML layout: PTID, Group, followed by the 32 dynamics features.
%
% This script does not recompute or modify harmonized values. Copy the
% exported XLSX file into the ML repository's Data/turbu_hopf folder and
% update only the input filename in the original Python runner.

clear; close all; clc;

sch1000_root = fileparts(fileparts(mfilename('fullpath')));
addpath(sch1000_root, '-begin');
P = setup_sch1000_paths();

%% Input and output names
harmonized_file = fullfile(P.harmonization_results, ...
    'Turbu_ComBat_ADNI3_allfeatures_N145.xlsx');
base_name = ...
    'ML_Input_ADNI3_4STAGINGBYABETA_ComBat_N145_withPTID_sch1000';

if ~isfile(harmonized_file)
    error(['Missing recomputed harmonized table: %s\n' ...
        'Run harmonization_allfeat/run_harmonization.py first.'], ...
        harmonized_file);
end
if ~isfolder(P.data_export)
    mkdir(P.data_export);
end

%% Expected original ML structure
group_labels = ["HC_ABneg", "HC_ABpos", "MCI_ABpos", "AD_ABpos"];
expected_counts = [51, 37, 31, 26];

lambda_values = [0.27, 0.24, 0.21, 0.18, 0.15, 0.12, 0.09, 0.06, 0.03, 0.01];
lambda_text = strrep(compose('%.2f', lambda_values), '.', '_');
turbu_names = "Turbu_lam_" + lambda_text;

% TransferLambda_sub(1,:) is an initialization row and was not included in
% the harmonized feature table or the original N145 ML input.
infoflow_names = "InfoFlow_lam_" + lambda_text(2:end);
infotransfer_names = "InfoTransfer_lam_" + lambda_text;
feature_names = [turbu_names, infoflow_names, infotransfer_names, ...
    "InfoCascade", "Metastability", "gKoP"];
expected_variables = ["PTID", "Group", feature_names];

%% Read the post-harmonization table
T_features = readtable(harmonized_file, 'VariableNamingRule', 'preserve');
actual_variables = string(T_features.Properties.VariableNames);

missing_variables = setdiff(expected_variables, actual_variables, 'stable');
unexpected_variables = setdiff(actual_variables, expected_variables, 'stable');
if ~isempty(missing_variables) || ~isempty(unexpected_variables)
    error(['The harmonized table does not match the original ML structure.\n' ...
        'Missing variables: %s\nUnexpected variables: %s'], ...
        strjoin(missing_variables, ', '), ...
        strjoin(unexpected_variables, ', '));
end

% Enforce the exact original column order without changing any values.
T_features = T_features(:, cellstr(expected_variables));
T_features.PTID = strip(string(T_features.PTID));
group = strip(string(T_features.Group));

%% Validate the canonical N145 cohort
if height(T_features) ~= sum(expected_counts)
    error('Expected 145 rows after harmonization; found %d.', ...
        height(T_features));
end
if any(ismissing(T_features.PTID)) || any(T_features.PTID == "")
    error('The harmonized table contains missing PTIDs.');
end
if numel(unique(T_features.PTID)) ~= height(T_features)
    error('The harmonized table contains duplicate PTIDs.');
end
if any(ismissing(group)) || any(group == "")
    error('The harmonized table contains missing group labels.');
end

unknown_groups = setdiff(unique(group), group_labels);
if ~isempty(unknown_groups)
    error('Unexpected group labels: %s', strjoin(unknown_groups, ', '));
end
observed_counts = arrayfun(@(label) nnz(group == label), group_labels);
if ~isequal(observed_counts, expected_counts)
    error(['Unexpected post-harmonization group counts. Found [%s]; ' ...
        'expected [%s].'], ...
        num2str(observed_counts), num2str(expected_counts));
end

X = T_features{:, cellstr(feature_names)};
if ~isnumeric(X) || ~isequal(size(X), [145, numel(feature_names)])
    error('Expected a numeric 145-by-%d feature matrix.', ...
        numel(feature_names));
end
if any(~isfinite(X), 'all')
    error('The harmonized feature matrix contains NaN or Inf values.');
end

% Match the categorical Group representation used by the original exporter.
T_features.Group = categorical(group, group_labels, group_labels);

%% Export ML-ready copies
mat_file = fullfile(P.data_export, [base_name '.mat']);
csv_file = fullfile(P.data_export, [base_name '.csv']);
xlsx_file = fullfile(P.data_export, [base_name '.xlsx']);

save(mat_file, 'T_features');
writetable(T_features, csv_file);
writetable(T_features, xlsx_file, 'FileType', 'spreadsheet');

fprintf('\nPost-harmonization ML export completed.\n');
fprintf('Source: %s\n', harmonized_file);
fprintf('Cohort: %d participants; features: %d.\n', ...
    height(T_features), numel(feature_names));
for g = 1:numel(group_labels)
    fprintf('  %s: %d\n', group_labels(g), observed_counts(g));
end
fprintf('XLSX to copy into the ML data folder:\n%s\n', xlsx_file);
