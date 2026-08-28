%% Export raw feature tables for ComBat harmonization
% The four standalone amyloid-status result files are the only turbulence
% sources. Output tables retain PTID and Group so row order is explicit.

clear; close all; clc;

sch1000_root = fileparts(fileparts(mfilename('fullpath')));
addpath(sch1000_root, '-begin');
P = setup_sch1000_paths();

group_keys = ["HC_ABetaNeg", "HC_ABetaPos", "MCI_ABetaPos", "AD_ABetaPos"];
group_labels = ["HC_ABneg", "HC_ABpos", "MCI_ABpos", "AD_ABpos"];
expected_counts = [54, 39, 33, 26];
measurement_files = [ ...
    "turbu_all_measurements_ADNI3_HC_ABetaNeg_MPRAGE_IRFSPGR_N238rev_sch1000.mat", ...
    "turbu_all_measurements_ADNI3_HC_ABetaPos_MPRAGE_IRFSPGR_N238rev_sch1000.mat", ...
    "turbu_all_measurements_ADNI3_MCI_ABetaPos_MPRAGE_IRFSPGR_N238rev_sch1000.mat", ...
    "turbu_all_measurements_ADNI3_AD_ABetaPos_MPRAGE_IRFSPGR_N238rev_sch1000.mat"];
node_files = [ ...
    "turbu_by_node_ADNI3_HC_ABetaNeg_MPRAGE_IRFSPGR_N238rev_sch1000.mat", ...
    "turbu_by_node_ADNI3_HC_ABetaPos_MPRAGE_IRFSPGR_N238rev_sch1000.mat", ...
    "turbu_by_node_ADNI3_MCI_ABetaPos_MPRAGE_IRFSPGR_N238rev_sch1000.mat", ...
    "turbu_by_node_ADNI3_AD_ABetaPos_MPRAGE_IRFSPGR_N238rev_sch1000.mat"];

lambda_values = [0.27, 0.24, 0.21, 0.18, 0.15, 0.12, 0.09, 0.06, 0.03, 0.01];
% Export the three node scales used by the legacy Neuromaps analysis. The
% historical suffixes encoded physical lambda x 100: lam1 = 0.01,
% lam3 = 0.03 and lam6 = 0.06. Because the pipeline stores lambdas in
% descending order, these values occupy MATLAB array indices 10, 9 and 8.
% Active filenames always state the physical value explicitly.
node_lambda_indices = [10, 9, 8];
node_lambda_labels = ["0_01", "0_03", "0_06"];
lambda_text = strrep(compose('%.2f', lambda_values), '.', '_');
feature_names = [ ...
    "Turbu_lam_" + lambda_text, ...
    "InfoFlow_lam_" + lambda_text(2:end), ...
    "InfoTransfer_lam_" + lambda_text, ...
    "InfoCascade", "Metastability", "gKoP"];

PTID_all = strings(0, 1);
Group_all = strings(0, 1);
X_all = zeros(0, numel(feature_names));
node_all = cell(1, numel(node_lambda_indices));
for k = 1:numel(node_lambda_indices)
    node_all{k} = zeros(0, 1000);
end

for g = 1:numel(group_keys)
    output_dir = P.output.(char(group_keys(g)));
    measurement_file = fullfile(output_dir, measurement_files(g));
    node_file = fullfile(output_dir, node_files(g));
    if ~isfile(measurement_file) || ~isfile(node_file)
        error('Missing standalone result files for %s.', group_labels(g));
    end

    R = load(measurement_file);
    N = load(node_file);
    required_measurements = {'PTID', 'LAMBDA', 'Turbulence_global_sub', ...
        'TransferLambda_sub', 'Transfer_sub', 'InformationCascade_sub', ...
        'Meta', 'gKoP'};
    missing = required_measurements(~isfield(R, required_measurements));
    if ~isempty(missing)
        error('%s is missing: %s', measurement_file, strjoin(missing, ', '));
    end
    if ~isfield(N, 'PTID') || ~isfield(N, 'Turbulence_node_sub')
        error('%s must contain PTID and Turbulence_node_sub.', node_file);
    end

    ptid = string(R.PTID(:));
    node_ptid = string(N.PTID(:));
    n_subjects = numel(ptid);
    if n_subjects ~= expected_counts(g)
        error('%s contains %d subjects; expected %d.', ...
            group_labels(g), n_subjects, expected_counts(g));
    end
    if ~isequal(ptid, node_ptid)
        error('Measurement and node PTIDs differ for %s.', group_labels(g));
    end
    if numel(unique(ptid)) ~= n_subjects
        error('%s contains duplicate PTIDs.', group_labels(g));
    end
    if ~isequal(R.LAMBDA(:).', lambda_values)
        error('%s has an unexpected lambda vector.', measurement_file);
    end
    if ~isequal(size(N.Turbulence_node_sub), [10, 1000, n_subjects])
        error('%s has unexpected node-turbulence dimensions.', node_file);
    end

    X_group = [ ...
        R.Turbulence_global_sub.', ...
        R.TransferLambda_sub(2:end, :).', ...
        R.Transfer_sub.', ...
        R.InformationCascade_sub(:), ...
        R.Meta(:), ...
        R.gKoP(:)];
    if size(X_group, 2) ~= numel(feature_names) || ...
            any(~isfinite(X_group), 'all')
        error('%s produced an invalid feature matrix.', group_labels(g));
    end

    PTID_all = [PTID_all; ptid]; %#ok<AGROW>
    Group_all = [Group_all; repmat(group_labels(g), n_subjects, 1)]; %#ok<AGROW>
    X_all = [X_all; X_group]; %#ok<AGROW>
    for k = 1:numel(node_lambda_indices)
        values = squeeze(N.Turbulence_node_sub(node_lambda_indices(k), :, :)).';
        if any(~isfinite(values), 'all')
            error('%s has non-finite node values at lambda index %d.', ...
                group_labels(g), node_lambda_indices(k));
        end
        node_all{k} = [node_all{k}; values];
    end
end

if numel(PTID_all) ~= 152 || numel(unique(PTID_all)) ~= 152
    error('Expected 152 unique participants across the four groups.');
end

if ~isfolder(P.harmonization_raw_inputs)
    mkdir(P.harmonization_raw_inputs);
end

T_features = array2table(X_all, 'VariableNames', cellstr(feature_names));
T_features = addvars(T_features, PTID_all, categorical(Group_all), ...
    'Before', 1, 'NewVariableNames', {'PTID', 'Group'});
allfeatures_file = fullfile(P.harmonization_raw_inputs, ...
    'turbu_raw_ADNI3_ABeta_N152_sch1000.xlsx');
writetable(T_features, allfeatures_file, 'FileType', 'spreadsheet');

parcel_names = cellstr("Schaefer_" + string(1:1000));
for k = 1:numel(node_lambda_indices)
    T_nodes = array2table(node_all{k}, 'VariableNames', parcel_names);
    T_nodes = addvars(T_nodes, PTID_all, categorical(Group_all), ...
        'Before', 1, 'NewVariableNames', {'PTID', 'Group'});
    node_output = fullfile(P.harmonization_raw_inputs, sprintf( ...
        'turbu_by_node_raw_lambda_%s_sch1000_HC_MCI_AD_ABeta_precombat.xlsx', ...
        node_lambda_labels(k)));
    writetable(T_nodes, node_output, 'FileType', 'spreadsheet');
end

fprintf(['Saved the all-feature and node N152 tables for lambda = ' ...
    '0.01, 0.03 and 0.06 to:\n%s\n'], ...
    P.harmonization_raw_inputs);
