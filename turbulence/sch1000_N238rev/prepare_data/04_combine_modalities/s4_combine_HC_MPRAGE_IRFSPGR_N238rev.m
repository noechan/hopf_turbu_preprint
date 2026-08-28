%% Combine HC MPRAGE and IRFSPGR time series across parcellations

clearvars
sch1000_root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(sch1000_root, '-begin');
P = setup_sch1000_paths();

% ---------------------------------------------------------
% Paths
% ---------------------------------------------------------
input_path = P.timeseries_inputs;

mprage_paths = { ...
    fullfile(input_path, 'ADNI3_HC_MPRAGE_batch_1_N238rev'), ...
    fullfile(input_path, 'ADNI3_HC_MPRAGE_batch_2_N238rev'), ...
    fullfile(input_path, 'ADNI3_HC_MPRAGE_batch_3_N238rev')};

irfspgr_path = fullfile(input_path, 'ADNI3_HC_IRFSPGR_batch_1_N238rev');

output_path = fullfile(input_path, 'ADNI3_HC_MPRAGE_IRFSPGR_N238rev');

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

    tseries_parts = cell(numel(mprage_paths) + 1, 1);

    % -----------------------------------------------------
    % Add HC MPRAGE batches 1-3
    % -----------------------------------------------------
    for b = 1:numel(mprage_paths)

        mprage_file = fullfile( ...
            mprage_paths{b}, ...
            sprintf('tseries_ADNI3_HC_MPRAGE_batch%d_%s_matching_QC.mat', b, parc));

        if ~exist(mprage_file, 'file')
            error('Missing MPRAGE batch %d file: %s', b, mprage_file);
        end

        tmp = load(mprage_file);

        if ~isfield(tmp, 'tseries')
            error('File does not contain variable "tseries": %s', mprage_file);
        end

        tseries_batch = tmp.tseries;

        fprintf('MPRAGE batch %d subjects: %d\n', b, numel(tseries_batch));

        tseries_parts{b} = tseries_batch;

    end

    % -----------------------------------------------------
    % Add HC IRFSPGR batch 1
    % -----------------------------------------------------
    irfspgr_file = fullfile( ...
        irfspgr_path, ...
        sprintf('tseries_ADNI3_HC_IRFSPGR_batch1_%s_matching_QC.mat', parc));

    if ~exist(irfspgr_file, 'file')
        error('Missing IRFSPGR file: %s', irfspgr_file);
    end

    tmp = load(irfspgr_file);

    if ~isfield(tmp, 'tseries')
        error('IRFSPGR file does not contain variable "tseries": %s', irfspgr_file);
    end

    tseries_irfspgr = tmp.tseries;

    fprintf('IRFSPGR batch 1 subjects: %d\n', numel(tseries_irfspgr));

    tseries_parts{end} = tseries_irfspgr;
    tseries = vertcat(tseries_parts{:});

    fprintf('Combined HC subjects: %d\n', numel(tseries));

    % -----------------------------------------------------
    % Save combined file
    % -----------------------------------------------------
    save_filename = sprintf( ...
        'tseries_ADNI3_HC_MPRAGE_IRFSPGR_%s_N238rev.mat', ...
        parc);

    save(fullfile(output_path, save_filename), 'tseries', '-v7.3');

    fprintf('Saved: %s\n', fullfile(output_path, save_filename));

end
