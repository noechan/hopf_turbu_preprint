%% Create timeseries for turbulence across multiple parcellations: AD_MPRAGE

clearvars
sch1000_root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(sch1000_root, '-begin');
P = setup_sch1000_paths();

% ---------------------------------------------------------
% Paths
% ---------------------------------------------------------
data_path = P.preprocessed_timeseries;
input_path = P.timeseries_inputs;

batch_path = { ...
    'ADNI3_AD_MPRAGE_batch_1_extended', ...
    'ADNI3_AD_MPRAGE_batch_2_extended', ...
    'ADNI3_AD_MPRAGE_batch_3_extended'};

tseries_path = { ...
    'ADNI3_AD_MPRAGE_batch_1_N238rev', ...
    'ADNI3_AD_MPRAGE_batch_2_N238rev', ...
    'ADNI3_AD_MPRAGE_batch_3_N238rev'};

% ---------------------------------------------------------
% Settings
% ---------------------------------------------------------
nTR = 197;
first_roi_idx = 168;

parcellations = struct();

parcellations(1).name      = 'sch400';
parcellations(1).nrois     = 400;
parcellations(1).start_idx = first_roi_idx;

parcellations(2).name      = 'sch1000';
parcellations(2).nrois     = 1000;
parcellations(2).start_idx = parcellations(1).start_idx + parcellations(1).nrois;

parcellations(3).name      = 'sch100';
parcellations(3).nrois     = 100;
parcellations(3).start_idx = parcellations(2).start_idx + parcellations(2).nrois;

parcellations(4).name      = 'dbs80';
parcellations(4).nrois     = 80;
parcellations(4).start_idx = parcellations(3).start_idx + parcellations(3).nrois;

parcellations(5).name      = 'glasser360';
parcellations(5).nrois     = 360;
parcellations(5).start_idx = parcellations(4).start_idx + parcellations(4).nrois;

% ---------------------------------------------------------
% Loop over batches
% ---------------------------------------------------------
for batch = 1:numel(batch_path)

    fprintf('\nProcessing AD_MPRAGE batch %d...\n', batch);

    data_path_batch  = fullfile(data_path, batch_path{batch});
    input_path_batch = fullfile(input_path, tseries_path{batch});

    if ~exist(input_path_batch, 'dir')
        mkdir(input_path_batch);
    end

    fnames = dir(data_path_batch);

    fnames_filt = fnames(~startsWith({fnames.name}, '._') & ...
                         ~startsWith({fnames.name}, '.')  & ...
                         ~startsWith({fnames.name}, '..'));

    connIDs_varname  = sprintf('connIDs_BIDS_MPRAGE_60_89_batch_%d_AD', batch);
    connIDs_filename = sprintf('%s.mat', connIDs_varname);

    load(fullfile(P.prepare_metadata, connIDs_filename), connIDs_varname);
    connIDs = eval(connIDs_varname);

    fnames_matching = fnames_filt(connIDs);

    fprintf('Number of matched subjects: %d\n', numel(fnames_matching));

    for p = 1:numel(parcellations)

        parc_name = parcellations(p).name;
        nrois     = parcellations(p).nrois;
        start_idx = parcellations(p).start_idx;
        end_idx   = start_idx + nrois - 1;

        fprintf('Extracting %s: data{%d:%d}\n', parc_name, start_idx, end_idx);

        tseries = cell(numel(fnames_matching), 1);

        for s = 1:numel(fnames_matching)

            subject_file = fullfile(data_path_batch, fnames_matching(s).name);
            tmp = load(subject_file);

            if ~isfield(tmp, 'data')
                error('File %s does not contain variable "data".', subject_file);
            end

            data = tmp.data;

            if numel(data) < end_idx
                error(['Subject %d file has only %d data entries, but %s requires index %d. ', ...
                       'Check whether you are using the *_extended folder and whether all ROIs were extracted.'], ...
                       s, numel(data), parc_name, end_idx);
            end

            data_parc = zeros(nrois, nTR);

            for roi = 1:nrois
                conn_idx = start_idx + roi - 1;
                roi_ts   = data{1, conn_idx};

                if numel(roi_ts) < nTR
                    error('Subject %d, ROI %d in %s has only %d time points; expected %d.', ...
                          s, roi, parc_name, numel(roi_ts), nTR);
                end

                data_parc(roi, :) = roi_ts(1:nTR)';
            end

            tseries{s} = data_parc;
        end

        save_filename = sprintf( ...
            'tseries_ADNI3_AD_MPRAGE_batch_%d_%s_matching_QC.mat', ...
            batch, parc_name);

        save(fullfile(input_path_batch, save_filename), 'tseries', '-v7.3');

        fprintf('Saved: %s\n', fullfile(input_path_batch, save_filename));

    end

end
