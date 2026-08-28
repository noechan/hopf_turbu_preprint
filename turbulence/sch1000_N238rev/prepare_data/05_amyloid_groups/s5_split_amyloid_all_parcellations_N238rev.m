%% Split ADNI3 time series by Aβ status across parcellations
% Also save PTIDs for each Aβ-defined group

clearvars; close all

sch1000_root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(sch1000_root, '-begin');
P = setup_sch1000_paths();

% ---------------------------------------------------------
% Parcellations
% ---------------------------------------------------------
parcellations = { ...
    'sch400', ...
    'sch1000', ...
    'sch100', ...
    'dbs80', ...
    'glasser360'};

% ---------------------------------------------------------
% Paths
% ---------------------------------------------------------
data_path = P.timeseries_inputs;

hc_path  = fullfile(data_path, 'ADNI3_HC_MPRAGE_IRFSPGR_N238rev');
mci_path = fullfile(data_path, 'ADNI3_MCI_MPRAGE_IRFSPGR_N238rev');
ad_path  = fullfile(data_path, 'ADNI3_AD_MPRAGE_IRFSPGR_N238rev');

ptid_path = P.prepare_metadata;

abeta_file = P.abeta_status_file;

% ---------------------------------------------------------
% Load Aβ status table
% ---------------------------------------------------------
abeta_table = readtable(abeta_file);

% Ensure PTID is cellstr/string-compatible
abeta_table.PTID = cellstr(string(abeta_table.PTID));
abeta_table.GROUP = cellstr(string(abeta_table.GROUP));

% ---------------------------------------------------------
% Load original PTIDs
% ---------------------------------------------------------
load(fullfile(ptid_path, 'PTID_ADNI3_HC_MPRAGE_IRFSPGR_all.mat'));
PTID_HC = cellstr(string(PTID));

load(fullfile(ptid_path, 'PTID_ADNI3_MCI_MPRAGE_IRFSPGR_all.mat'));
PTID_MCI = cellstr(string(PTID));

load(fullfile(ptid_path, 'PTID_ADNI3_AD_MPRAGE_IRFSPGR_all.mat'));
PTID_AD = cellstr(string(PTID));

% ---------------------------------------------------------
% Helper: masks by Aβ status aligned to group PTIDs
% ---------------------------------------------------------
get_abeta_masks = @(group_ids, gname) deal( ...
    ismember(group_ids, abeta_table.PTID(strcmpi(abeta_table.GROUP, gname) & abeta_table.ABeta_pvc == 1)), ...
    ismember(group_ids, abeta_table.PTID(strcmpi(abeta_table.GROUP, gname) & abeta_table.ABeta_pvc == 0)) ...
    );

[isAbetaPos_HC,  isAbetaNeg_HC] = get_abeta_masks(PTID_HC,  'HC');
[isAbetaPos_MCI, isAbetaNeg_MCI] = get_abeta_masks(PTID_MCI, 'MCI');
[isAbetaPos_AD,  isAbetaNeg_AD]  = get_abeta_masks(PTID_AD,  'AD');

% ---------------------------------------------------------
% Define PTID subsets
% ---------------------------------------------------------
PTID_HC_ABneg  = PTID_HC(isAbetaNeg_HC);
PTID_HC_ABpos  = PTID_HC(isAbetaPos_HC);

PTID_MCI_ABpos = PTID_MCI(isAbetaPos_MCI);
PTID_AD_ABpos  = PTID_AD(isAbetaPos_AD);

% Optional, saved for transparency
PTID_MCI_ABneg = PTID_MCI(isAbetaNeg_MCI);
PTID_AD_ABneg  = PTID_AD(isAbetaNeg_AD);

% ---------------------------------------------------------
% Print summary
% ---------------------------------------------------------
fprintf('HC Aβ−:   %d\n', numel(PTID_HC_ABneg));
fprintf('HC Aβ+:   %d\n', numel(PTID_HC_ABpos));
fprintf('MCI Aβ+:  %d\n', numel(PTID_MCI_ABpos));
fprintf('AD Aβ+:   %d\n', numel(PTID_AD_ABpos));

fprintf('MCI Aβ−:  %d\n', numel(PTID_MCI_ABneg));
fprintf('AD Aβ−:   %d\n', numel(PTID_AD_ABneg));

% ---------------------------------------------------------
% Save PTIDs once
% ---------------------------------------------------------
save(fullfile(ptid_path, 'PTID_ADNI3_HC_ABetaNeg_MPRAGE_IRFSPGR_all.mat'), ...
    'PTID_HC_ABneg');

save(fullfile(ptid_path, 'PTID_ADNI3_HC_ABetaPos_MPRAGE_IRFSPGR_all.mat'), ...
    'PTID_HC_ABpos');

save(fullfile(ptid_path, 'PTID_ADNI3_MCI_ABetaPos_MPRAGE_IRFSPGR_all.mat'), ...
    'PTID_MCI_ABpos');

save(fullfile(ptid_path, 'PTID_ADNI3_AD_ABetaPos_MPRAGE_IRFSPGR_all.mat'), ...
    'PTID_AD_ABpos');

save(fullfile(ptid_path, 'PTID_ADNI3_MCI_ABetaNeg_MPRAGE_IRFSPGR_all.mat'), ...
    'PTID_MCI_ABneg');

save(fullfile(ptid_path, 'PTID_ADNI3_AD_ABetaNeg_MPRAGE_IRFSPGR_all.mat'), ...
    'PTID_AD_ABneg');

fprintf('\nSaved PTID files by diagnostic group and Aβ status.\n');

% ---------------------------------------------------------
% Loop over parcellations
% ---------------------------------------------------------
for p = 1:numel(parcellations)

    parc = parcellations{p};

    fprintf('\nProcessing parcellation: %s\n', parc);

    %% -----------------------------------------------------
    % HC
    % ------------------------------------------------------
    hc_file = fullfile(hc_path, ...
        sprintf('tseries_ADNI3_HC_MPRAGE_IRFSPGR_%s_N238rev.mat', parc));

    if ~exist(hc_file, 'file')
        error('Missing HC file: %s', hc_file);
    end

    tmp = load(hc_file);

    if ~isfield(tmp, 'tseries')
        error('HC file does not contain variable "tseries": %s', hc_file);
    end

    tseries_HC = tmp.tseries;

    if numel(tseries_HC) ~= numel(PTID_HC)
        error('HC mismatch: %d tseries entries but %d PTIDs.', ...
            numel(tseries_HC), numel(PTID_HC));
    end

    tseries_HC_ABneg = tseries_HC(isAbetaNeg_HC, :);
    tseries_HC_ABpos = tseries_HC(isAbetaPos_HC, :);

    save(fullfile(hc_path, ...
        sprintf('tseries_ADNI3_HC_ABetaNeg_MPRAGE_IRFSPGR_%s_N238rev.mat', parc)), ...
        'tseries_HC_ABneg', 'PTID_HC_ABneg', '-v7.3');

    save(fullfile(hc_path, ...
        sprintf('tseries_ADNI3_HC_ABetaPos_MPRAGE_IRFSPGR_%s_N238rev.mat', parc)), ...
        'tseries_HC_ABpos', 'PTID_HC_ABpos', '-v7.3');

    fprintf('Saved HC Aβ− and HC Aβ+ for %s\n', parc);

    %% -----------------------------------------------------
    % MCI
    % ------------------------------------------------------
    mci_file = fullfile(mci_path, ...
        sprintf('tseries_ADNI3_MCI_MPRAGE_IRFSPGR_%s_N238rev.mat', parc));

    if ~exist(mci_file, 'file')
        error('Missing MCI file: %s', mci_file);
    end

    tmp = load(mci_file);

    if ~isfield(tmp, 'tseries')
        error('MCI file does not contain variable "tseries": %s', mci_file);
    end

    tseries_MCI = tmp.tseries;

    if numel(tseries_MCI) ~= numel(PTID_MCI)
        error('MCI mismatch: %d tseries entries but %d PTIDs.', ...
            numel(tseries_MCI), numel(PTID_MCI));
    end

    tseries_MCI_ABpos = tseries_MCI(isAbetaPos_MCI, :);

    save(fullfile(mci_path, ...
        sprintf('tseries_ADNI3_MCI_ABetaPos_MPRAGE_IRFSPGR_%s_N238rev.mat', parc)), ...
        'tseries_MCI_ABpos', 'PTID_MCI_ABpos', '-v7.3');

    fprintf('Saved MCI Aβ+ for %s\n', parc);

    %% -----------------------------------------------------
    % AD
    % ------------------------------------------------------
    ad_file = fullfile(ad_path, ...
        sprintf('tseries_ADNI3_AD_MPRAGE_IRFSPGR_%s_N238rev.mat', parc));

    if ~exist(ad_file, 'file')
        error('Missing AD file: %s', ad_file);
    end

    tmp = load(ad_file);

    if ~isfield(tmp, 'tseries')
        error('AD file does not contain variable "tseries": %s', ad_file);
    end

    tseries_AD = tmp.tseries;

    if numel(tseries_AD) ~= numel(PTID_AD)
        error('AD mismatch: %d tseries entries but %d PTIDs.', ...
            numel(tseries_AD), numel(PTID_AD));
    end

    tseries_AD_ABpos = tseries_AD(isAbetaPos_AD, :);

    save(fullfile(ad_path, ...
        sprintf('tseries_ADNI3_AD_ABetaPos_MPRAGE_IRFSPGR_%s_N238rev.mat', parc)), ...
        'tseries_AD_ABpos', 'PTID_AD_ABpos', '-v7.3');

    fprintf('Saved AD Aβ+ for %s\n', parc);

end

fprintf('\nDone.\n');
