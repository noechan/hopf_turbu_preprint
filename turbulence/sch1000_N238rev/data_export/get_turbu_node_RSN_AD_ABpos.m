%% Extract raw turbulence values for AD Aβ+ only and export to Excel
% Rows = PTID
% Columns = Schaefer parcels (Schaefer_1 ... Schaefer_nNodes)

clear all; close all; clc;

sch1000_root = fileparts(fileparts(mfilename('fullpath')));
addpath(sch1000_root, '-begin');
P = setup_sch1000_paths();
output_path = P.timeseries_outputs;

%% 1. Load Aβ status table
abeta_table = readtable(P.abeta_status_file);

%% 2. Load PTIDs (only AD needed, but HC/MCI kept for consistency if extended later)
ptid_path = P.prepare_metadata;

% HC (not used here, but can be useful if extending the script)
load(fullfile(ptid_path, 'PTID_ADNI3_HC_MPRAGE_IRFSPGR_all.mat'))
PTID_HC = cellstr(PTID);

% MCI (not used here, but can be useful if extending the script)
load(fullfile(ptid_path, 'PTID_ADNI3_MCI_MPRAGE_IRFSPGR_all.mat'))
PTID_MCI = cellstr(PTID);

% AD (used in this script)
load(fullfile(ptid_path, 'PTID_ADNI3_AD_MPRAGE_IRFSPGR_all.mat'))
PTID_AD = cellstr(PTID);

%% 3. General Aβ mask helper (HC, MCI, AD)
% NOTE: If GROUP uses 'CN' instead of 'HC', replace 'HC' with 'CN' below.
get_abeta_masks = @(group_ids, group_name) deal( ...
    ismember(group_ids, abeta_table.PTID(strcmp(abeta_table.GROUP, group_name) & abeta_table.ABeta_pvc == 1)), ...
    ismember(group_ids, abeta_table.PTID(strcmp(abeta_table.GROUP, group_name) & abeta_table.ABeta_pvc == 0)));

[isAbetaPos_HC,  isAbetaNeg_HC]  = get_abeta_masks(PTID_HC,  'HC');  %#ok<NASGU>
[isAbetaPos_MCI, ~]              = get_abeta_masks(PTID_MCI, 'MCI'); %#ok<NASGU>
[isAbetaPos_AD,  ~]              = get_abeta_masks(PTID_AD,  'AD');  % AD Aβ+

%% 4. Load turbulence data for AD
turbu_file = fullfile(output_path, ...
    'Turbulence_ADNI3_AD_MPRAGE_IRFSPGR_N238rev_sch1000', ...
    'turbu_by_node_ADNI3_AD_MPRAGE_IRFSPGR_N238rev_sch1000.mat');

load(turbu_file, 'Turbulence_node_sub');   % expected: (nLambda x nNodes x nSubj_AD)

turbu = Turbulence_node_sub;
[nLambda, nNodes, nSubj_AD] = size(turbu);

fprintf('Loaded AD turbulence: %d lambda x %d nodes x %d AD subjects\n', ...
        nLambda, nNodes, nSubj_AD);

% Optional consistency check
if nSubj_AD ~= numel(PTID_AD)
    warning('Number of AD subjects in turbulence (%d) does not match PTID_AD length (%d).', ...
        nSubj_AD, numel(PTID_AD));
end

%% 5. Restrict to AD Aβ+ only (AD+)
mask_ADpos = isAbetaPos_AD(:);   % ensure column vector

% Safety check: mask length vs subjects
if numel(mask_ADpos) ~= nSubj_AD
    error('Length of isAbetaPos_AD (%d) does not match nSubj_AD (%d).', ...
        numel(mask_ADpos), nSubj_AD);
end

turbu_ADpos = turbu(:, :, mask_ADpos);    % (nLambda x nNodes x nSubj_ADpos)
PTID_ADpos  = PTID_AD(mask_ADpos);
nSubj_ADpos = numel(PTID_ADpos);

fprintf('AD Aβ+ subjects: %d\n', nSubj_ADpos);

%% 6. Prepare labels and output Excel file
% Hard-coded Schaefer parcel labels instead of node_#
Schaefer_labels = arrayfun(@(i) sprintf('Schaefer_%d', i), 1:nNodes, ...
                           'UniformOutput', false);

lambdaNames = arrayfun(@(i) sprintf('lambda_%d', i), 1:nLambda, ...
                       'UniformOutput', false);

outFile = fullfile(P.data_export, ...
    'turbulence_AD_AbetaPos_by_Schaefer_per_lambda_RAW.xlsx');
% If you prefer to write to the outputs folder instead, use:
% outFile = fullfile(output_path, ...
%    'turbulence_AD_AbetaPos_by_Schaefer_per_lambda_RAW.xlsx');

%% 7. Write one sheet per lambda (rows = PTID, columns = Schaefer parcels)
for l = 1:nLambda %1: lambda=0.27...10: Lambda=0.01
    % turbu(lambda, node, subj) → (nodes x subjects)
    % We want subjects x nodes → transpose
    lambdaData = squeeze(turbu_ADpos(l, :, :)).';   % size: nSubj_ADpos x nNodes

    T = array2table(lambdaData, ...
                    'VariableNames', Schaefer_labels, ...
                    'RowNames',      PTID_ADpos);   % PTIDs preserved exactly

    sheetName = sprintf('lambda_%d', l);
    writetable(T, outFile, 'Sheet', sheetName, 'WriteRowNames', true);
end

%% 8. Summary sheet: mean across AD Aβ+ per Schaefer parcel and lambda
% mean over subjects (3rd dim): result nLambda x nNodes
meanAcrossSubj = squeeze(mean(turbu_ADpos, 3));   % (nLambda x nNodes)

% We want rows = Schaefer parcels, columns = lambdas
meanAcrossSubj_T = meanAcrossSubj.';              % (nNodes x nLambda)

Tmean = array2table(meanAcrossSubj_T, ...
                    'VariableNames', lambdaNames, ...
                    'RowNames',      Schaefer_labels);

writetable(Tmean, outFile, 'Sheet', 'mean_per_Schaefer', 'WriteRowNames', true);

fprintf('Raw AD Aβ+ turbulence values (rows=PTID, cols=Schaefer) exported to:\n%s\n', outFile);
