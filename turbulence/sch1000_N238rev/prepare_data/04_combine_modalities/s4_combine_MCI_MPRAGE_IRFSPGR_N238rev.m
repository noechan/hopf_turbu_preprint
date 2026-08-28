%% Combine MCI MPRAGE and IRFSPGR time series across parcellations

clearvars
sch1000_root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(sch1000_root, '-begin');
P = setup_sch1000_paths();

% ---------------------------------------------------------
% Paths
% ---------------------------------------------------------
input_path = P.timeseries_inputs;

mprage_path  = fullfile(input_path, 'ADNI3_MCI_MPRAGE_batch_1_N238rev');
irfspgr_path = fullfile(input_path, 'ADNI3_MCI_IRFSPGR_batch_1_N238rev');

output_path = fullfile(input_path, 'ADNI3_MCI_MPRAGE_IRFSPGR_N238rev');

if ~exist(output_path, 'dir')
    mkdir(output_path);
end

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
% Combine per parcellation
% ---------------------------------------------------------
for p = 1:numel(parcellations)

    parc = parcellations{p};

    fprintf('\nCombining parcellation: %s\n', parc);

    % Input filenames
    mprage_file = fullfile( ...
        mprage_path, ...
        sprintf('tseries_ADNI3_MCI_MPRAGE_batch_1_%s_matching_QC.mat', parc));

    irfspgr_file = fullfile( ...
        irfspgr_path, ...
        sprintf('tseries_ADNI3_MCI_IRFSPGR_batch_1_%s_matching_QC.mat', parc));

    % Check files exist
    if ~exist(mprage_file, 'file')
        error('Missing MPRAGE file: %s', mprage_file);
    end

    if ~exist(irfspgr_file, 'file')
        error('Missing IRFSPGR file: %s', irfspgr_file);
    end

    % Load MPRAGE
    tmp_mprage = load(mprage_file);

    if ~isfield(tmp_mprage, 'tseries')
        error('MPRAGE file does not contain variable "tseries": %s', mprage_file);
    end

    tseries_mprage = tmp_mprage.tseries;

    % Load IRFSPGR
    tmp_irfspgr = load(irfspgr_file);

    if ~isfield(tmp_irfspgr, 'tseries')
        error('IRFSPGR file does not contain variable "tseries": %s', irfspgr_file);
    end

    tseries_irfspgr = tmp_irfspgr.tseries;

    % Combine
    tseries = [tseries_mprage; tseries_irfspgr];

    fprintf('MPRAGE subjects:  %d\n', numel(tseries_mprage));
    fprintf('IRFSPGR subjects: %d\n', numel(tseries_irfspgr));
    fprintf('Combined subjects: %d\n', numel(tseries));

    % Save combined file
    save_filename = sprintf( ...
        'tseries_ADNI3_MCI_MPRAGE_IRFSPGR_%s_N238rev.mat', ...
        parc);

    save(fullfile(output_path, save_filename), 'tseries', '-v7.3');

    fprintf('Saved: %s\n', fullfile(output_path, save_filename));

end
