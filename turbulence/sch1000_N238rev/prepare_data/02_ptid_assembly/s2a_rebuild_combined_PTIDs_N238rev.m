%% Rebuild combined PTID files in the same order as the combined time series

clear; clc;
sch1000_root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(sch1000_root, '-begin');
P = setup_sch1000_paths();
if ~isfolder(P.prepare_metadata), mkdir(P.prepare_metadata); end

groups = {
    'HC', 109, { ...
        'PTID_BIDS_MPRAGE_60_89_batch_1_HC.mat', ...
        'PTID_BIDS_MPRAGE_60_89_batch_2_HC.mat', ...
        'PTID_BIDS_MPRAGE_60_89_batch_3_HC.mat', ...
        'PTID_BIDS_IRFSPGR_60_89_batch_1_HC.mat'};
    'MCI', 90, { ...
        'PTID_BIDS_MPRAGE_60_89_batch_1_MCI.mat', ...
        'PTID_BIDS_IRFSPGR_60_89_batch_1_MCI.mat'};
    'AD', 39, { ...
        'PTID_BIDS_MPRAGE_60_89_batch_1_AD.mat', ...
        'PTID_BIDS_MPRAGE_60_89_batch_2_AD.mat', ...
        'PTID_BIDS_MPRAGE_60_89_batch_3_AD.mat', ...
        'PTID_BIDS_IRFSPGR_60_89_batch_1_AD.mat', ...
        'PTID_BIDS_IRFSPGR_60_89_batch_2_AD.mat'};
    };

for group_index = 1:size(groups, 1)
    group_name = groups{group_index, 1};
    expected_count = groups{group_index, 2};
    batch_files = groups{group_index, 3};
    combined_ids = strings(0, 1);

    fprintf('\nRebuilding %s PTIDs...\n', group_name);
    for batch_index = 1:numel(batch_files)
        input_file = fullfile(P.prepare_metadata, batch_files{batch_index});
        if ~isfile(input_file)
            error('Missing batch PTID file: %s', input_file);
        end

        batch_data = load(input_file);
        variable_names = fieldnames(batch_data);
        if numel(variable_names) ~= 1
            error('Expected one variable in %s; found %d.', ...
                input_file, numel(variable_names));
        end

        batch_ids = string(batch_data.(variable_names{1}));
        batch_ids = batch_ids(:);
        if any(ismissing(batch_ids) | strlength(strtrim(batch_ids)) == 0)
            error('Missing or empty PTID in %s.', input_file);
        end

        fprintf('  %-55s %3d subjects\n', batch_files{batch_index}, ...
            numel(batch_ids));
        combined_ids = [combined_ids; batch_ids]; %#ok<AGROW>
    end

    if numel(combined_ids) ~= expected_count
        error('Unexpected %s count: %d; expected %d.', group_name, ...
            numel(combined_ids), expected_count);
    end
    if numel(unique(combined_ids)) ~= numel(combined_ids)
        error('Duplicate PTIDs found while combining the %s batches.', group_name);
    end

    % Retain the historical variable name and cell-array representation used
    % by all downstream scripts.
    PTID = cellstr(combined_ids);
    output_file = fullfile(P.prepare_metadata, sprintf( ...
        'PTID_ADNI3_%s_MPRAGE_IRFSPGR_all.mat', group_name));
    save(output_file, 'PTID');
    fprintf('Saved %d ordered PTIDs: %s\n', numel(PTID), output_file);
end

fprintf('\nCombined PTID reconstruction passed for all 238 participants.\n');
