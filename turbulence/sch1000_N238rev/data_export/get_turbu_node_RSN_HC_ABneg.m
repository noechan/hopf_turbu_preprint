% Rows = PTID
% Columns = Schaefer parcels (Schaefer_1 ... Schaefer_nNodes)

clear all; close all; clc;

sch1000_root = fileparts(fileparts(mfilename('fullpath')));
addpath(sch1000_root, '-begin');
P = setup_sch1000_paths();
output_path = P.timeseries_outputs;

%% 1. Load Aβ status table
abeta_table = readtable(P.abeta_status_file);

%% 2. Load PTIDs (HC, MCI, AD)
ptid_path = P.prepare_metadata;

% HC
load(fullfile(ptid_path, 'PTID_ADNI3_HC_MPRAGE_IRFSPGR_all.mat'))
PTID_HC = cellstr(PTID);

% MCI (not used here, but kept for symmetry)
load(fullfile(ptid_path, 'PTID_ADNI3_MCI_MPRAGE_IRFSPGR_all.mat'))
PTID_MCI = cellstr(PTID); %#ok<NASGU>

% AD (not used here, but kept for symmetry)
load(fullfile(ptid_path, 'PTID_ADNI3_AD_MPRAGE_IRFSPGR_all.mat'))
PTID_AD = cellstr(PTID); %#ok<NASGU>

%% 3. General Aβ mask helper (HC, MCI, AD)
% NOTE: If GROUP uses 'CN' instead of 'HC', replace 'HC' with 'CN' below.
get_abeta_masks = @(group_ids, group_name) deal( ...
    ismember(group_ids, abeta_table.PTID(strcmp(abeta_table.GROUP, group_name) & abeta_table.ABeta_pvc == 1)), ...
    ismember(group_ids, abeta_table.PTID(strcmp(abeta_table.GROUP, group_name) & abeta_table.ABeta_pvc == 0)));

[isAbetaPos_HC,  isAbetaNeg_HC]  = get_abeta_masks(PTID_HC,  'HC');
[isAbetaPos_MCI, ~]              = get_abeta_masks(PTID_MCI, 'MCI'); %#ok<NASGU>
[isAbetaPos_AD,  ~]              = get_abeta_masks(PTID_AD,  'AD');  %#ok<NASGU>

%% 4. Load turbulence data for HC
turbu_file = fullfile(output_path, ...
    'Turbulence_ADNI3_HC_MPRAGE_IRFSPGR_N238rev_sch1000', ...
    'turbu_by_node_ADNI3_HC_MPRAGE_IRFSPGR_N238rev_sch1000.mat');

load(turbu_file, 'Turbulence_node_sub');   % expected: (nLambda x nNodes x nSubj_HC)

turbu = Turbulence_node_sub;
[nLambda, nNodes, nSubj_HC] = size(turbu);

fprintf('Loaded HC turbulence: %d lambda x %d nodes x %d HC subjects\n', ...
        nLambda, nNodes, nSubj_HC);

% Optional consistency check
if nSubj_HC ~= numel(PTID_HC)
    warning('Number of HC subjects in turbulence (%d) does not match PTID_HC length (%d).', ...
        nSubj_HC, numel(PTID_HC));
end

%% 5. Restrict to HC Aβ− only (HC−)
mask_HCneg = isAbetaNeg_HC(:);   % ensure column vector

% Safety check: mask length vs subjects
if numel(mask_HCneg) ~= nSubj_HC
    error('Length of isAbetaNeg_HC (%d) does not match nSubj_HC (%d).', ...
        numel(mask_HCneg), nSubj_HC);
end

turbu_HCneg = turbu(:, :, mask_HCneg);    % (nLambda x nNodes x nSubj_HCneg)
PTID_HCneg  = PTID_HC(mask_HCneg);
nSubj_HCneg = numel(PTID_HCneg);

fprintf('HC Aβ- subjects: %d\n', nSubj_HCneg);

%% 6. Prepare labels and output Excel file
% Hard-coded Schaefer parcel labels instead of node_#
Schaefer_labels = arrayfun(@(i) sprintf('Schaefer_%d', i), 1:nNodes, ...
                           'UniformOutput', false);

lambdaNames = arrayfun(@(i) sprintf('lambda_%d', i), 1:nLambda, ...
                       'UniformOutput', false);

outFile = fullfile(P.data_export, ...
    'turbulence_HC_AbetaNeg_by_Schaefer_per_lambda_RAW.xlsx');
% If you prefer to write to the outputs folder instead, use:
% outFile = fullfile(output_path, ...
%    'turbulence_HC_AbetaNeg_by_Schaefer_per_lambda_RAW.xlsx');

%% 7. Write one sheet per lambda (rows = PTID, columns = Schaefer parcels)
for l = 1:nLambda   % 1: lambda=0.27 ... 10: lambda=0.01
    % turbu(lambda, node, subj) → (nodes x subjects)
    % We want subjects x nodes → transpose
    lambdaData = squeeze(turbu_HCneg(l, :, :)).';   % size: nSubj_HCneg x nNodes

    T = array2table(lambdaData, ...
                    'VariableNames', Schaefer_labels, ...
                    'RowNames',      PTID_HCneg);   % PTIDs preserved exactly

    sheetName = sprintf('lambda_%d', l);
    writetable(T, outFile, 'Sheet', sheetName, 'WriteRowNames', true);
end

%% 8. Summary sheet: mean across HC Aβ− per Schaefer parcel and lambda
% mean over subjects (3rd dim): result nLambda x nNodes
meanAcrossSubj = squeeze(mean(turbu_HCneg, 3));   % (nLambda x nNodes)

% We want rows = Schaefer parcels, columns = lambdas
meanAcrossSubj_T = meanAcrossSubj.';              % (nNodes x nLambda)

Tmean = array2table(meanAcrossSubj_T, ...
                    'VariableNames', lambdaNames, ...
                    'RowNames',      Schaefer_labels);

writetable(Tmean, outFile, 'Sheet', 'mean_per_Schaefer', 'WriteRowNames', true);

fprintf('Raw HC Aβ- turbulence values (rows=PTID, cols=Schaefer) exported to:\n%s\n', outFile);
