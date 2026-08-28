%% Combine AD MPRAGE and IRFSPGR time series across parcellations

clearvars
sch1000_root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(sch1000_root, '-begin');
P = setup_sch1000_paths();

% ---------------------------------------------------------
% Paths
% ---------------------------------------------------------
input_path = P.timeseries_inputs;

mprage_paths = { ...
    fullfile(input_path, 'ADNI3_AD_MPRAGE_batch_1_N238rev'), ...
    fullfile(input_path, 'ADNI3_AD_MPRAGE_batch_2_N238rev'), ...
    fullfile(input_path, 'ADNI3_AD_MPRAGE_batch_3_N238rev')};

irfspgr_paths = { ...
    fullfile(input_path, 'ADNI3_AD_IRFSPGR_batch_1_N238rev'), ...
    fullfile(input_path, 'ADNI3_AD_IRFSPGR_batch_2_N238rev')};

output_path = fullfile(input_path, 'ADNI3_AD_MPRAGE_IRFSPGR_N238rev');

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

    tseries_parts = cell(numel(mprage_paths) + numel(irfspgr_paths), 1);

    % -----------------------------------------------------
    % Add AD MPRAGE batches 1-3
    % -----------------------------------------------------
    for b = 1:numel(mprage_paths)

        mprage_file = fullfile( ...
            mprage_paths{b}, ...
            sprintf('tseries_ADNI3_AD_MPRAGE_batch_%d_%s_matching_QC.mat', b, parc));

        if ~exist(mprage_file, 'file')
            error('Missing MPRAGE batch %d file: %s', b, mprage_file);
        end

        tmp = load(mprage_file);

        if ~isfield(tmp, 'tseries')
            error('MPRAGE batch %d file does not contain variable "tseries": %s', b, mprage_file);
        end

        tseries_batch = tmp.tseries;

        fprintf('MPRAGE batch %d subjects: %d\n', b, numel(tseries_batch));

        tseries_parts{b} = tseries_batch;

    end

    % -----------------------------------------------------
    % Add AD IRFSPGR batches 1-2
    % -----------------------------------------------------
    for b = 1:numel(irfspgr_paths)

        irfspgr_file = fullfile( ...
            irfspgr_paths{b}, ...
            sprintf('tseries_ADNI3_AD_IRFSPGR_batch_%d_%s_matching_QC.mat', b, parc));

        if ~exist(irfspgr_file, 'file')
            error('Missing IRFSPGR batch %d file: %s', b, irfspgr_file);
        end

        tmp = load(irfspgr_file);

        if ~isfield(tmp, 'tseries')
            error('IRFSPGR batch %d file does not contain variable "tseries": %s', b, irfspgr_file);
        end

        tseries_batch = tmp.tseries;

        fprintf('IRFSPGR batch %d subjects: %d\n', b, numel(tseries_batch));

        tseries_parts{numel(mprage_paths) + b} = tseries_batch;

    end

    tseries = vertcat(tseries_parts{:});
    fprintf('Combined AD subjects: %d\n', numel(tseries));

    % -----------------------------------------------------
    % Save combined file
    % -----------------------------------------------------
    save_filename = sprintf( ...
        'tseries_ADNI3_AD_MPRAGE_IRFSPGR_%s_N238rev.mat', ...
        parc);

    save(fullfile(output_path, save_filename), 'tseries', '-v7.3');

    fprintf('Saved: %s\n', fullfile(output_path, save_filename));

end
