function report = validate_sch1000_data(deep)
%VALIDATE_SCH1000_DATA Check cohort membership and MAT dimensions read-only.
%
% validate_sch1000_data(false) checks tables, PTIDs, input dimensions, and
% existing output dimensions. validate_sch1000_data(true) additionally
% reads every time-series cell and requires a 1000 x 197 numeric matrix.

if nargin < 1
    deep = false;
end

P = setup_sch1000_paths();
groups = {'HC', 'MCI', 'AD'};
selection_files = {
    'ADNI3_HC_ALL_MRI_sc_analysis_N238rev_idx4MRIselection.xlsx'
    'ADNI3_MCI_ALL_MRI_sc_analysis_N238rev_idx4MRIselection.xlsx'
    'ADNI3_AD_ALL_MRI_sc_analysis_N238rev_idx4MRIselection.xlsx'
    };
sheets = {
    'HC_removed_N238rev_idx'
    'MCI_removed_N238rev_idx'
    'AD_removed_N238rev_idx'
    };

abeta = readtable(P.abeta_status_file);
abeta.PTID = string(abeta.PTID);
abeta.GROUP = upper(string(abeta.GROUP));

retained = zeros(3, 1);
input_subjects = zeros(3, 1);
ptid_subjects = zeros(3, 1);
abeta_negative = zeros(3, 1);
abeta_positive = zeros(3, 1);
abeta_missing = zeros(3, 1);
existing_output_subjects = nan(3, 1);
node_output_subjects = nan(3, 1);

for i = 1:numel(groups)
    group = groups{i};
    selection = readtable(fullfile(P.participants_cross_sectional, ...
        selection_files{i}), 'Sheet', sheets{i});
    retained(i) = sum(selection.Total_removed ~= 1);

    input_file = fullfile(P.input.(group), sprintf( ...
        'tseries_ADNI3_%s_MPRAGE_IRFSPGR_sch1000_N238rev.mat', group));
    info = whos('-file', input_file, 'tseries');
    assert(~isempty(info), 'Missing tseries variable in %s', input_file);
    assert(info.size(2) == 1, 'Expected an N-by-1 tseries cell array in %s', input_file);
    input_subjects(i) = info.size(1);

    ptid_file = fullfile(P.prepare_metadata, sprintf( ...
        'PTID_ADNI3_%s_MPRAGE_IRFSPGR_all.mat', group));
    ptid_data = load(ptid_file);
    fields = fieldnames(ptid_data);
    assert(numel(fields) == 1, 'Expected one PTID variable in %s', ptid_file);
    ptids = string(ptid_data.(fields{1}));
    ptid_subjects(i) = numel(ptids);

    group_abeta = abeta(abeta.GROUP == group & ismember(abeta.PTID, ptids), :);
    abeta_negative(i) = sum(group_abeta.ABeta_pvc == 0);
    abeta_positive(i) = sum(group_abeta.ABeta_pvc == 1);
    abeta_missing(i) = sum(isnan(group_abeta.ABeta_pvc));

    assert(retained(i) == input_subjects(i), ...
        '%s retained count does not match its tseries count.', group);
    assert(ptid_subjects(i) == input_subjects(i), ...
        '%s PTID count does not match its tseries count.', group);
    assert(height(group_abeta) == input_subjects(i), ...
        '%s amyloid rows do not cover every selected PTID exactly once.', group);

    if deep
        data = matfile(input_file);
        for subject = 1:input_subjects(i)
            one_cell = data.tseries(subject, 1);
            ts = one_cell{1};
            assert(isnumeric(ts) && isequal(size(ts), [1000, 197]), ...
                '%s subject %d has size %s; expected [1000 197].', ...
                group, subject, mat2str(size(ts)));
        end
    end

    output_file = fullfile(P.output.(group), sprintf( ...
        'turbu_all_measurements_ADNI3_%s_MPRAGE_IRFSPGR_N238rev_sch1000.mat', ...
        group));
    if isfile(output_file)
        output_info = whos('-file', output_file, 'Turbulence_global_sub');
        if ~isempty(output_info)
            existing_output_subjects(i) = output_info.size(2);
            assert(existing_output_subjects(i) == input_subjects(i), ...
                '%s existing turbulence output has the wrong subject count.', group);
        end
    end

    node_output_file = fullfile(P.output.(group), sprintf( ...
        'turbu_by_node_ADNI3_%s_MPRAGE_IRFSPGR_N238rev_sch1000.mat', group));
    if isfile(node_output_file)
        node_info = whos('-file', node_output_file, 'Turbulence_node_sub');
        if ~isempty(node_info)
            assert(node_info.size(1) == 10 && node_info.size(2) == 1000, ...
                '%s node-wise turbulence must be 10 x 1000 x N.', group);
            node_output_subjects(i) = node_info.size(3);
            assert(node_output_subjects(i) == input_subjects(i), ...
                '%s node-wise output has the wrong subject count.', group);
        end
    end
end

report = table(string(groups(:)), retained, input_subjects, ptid_subjects, ...
    abeta_negative, abeta_positive, abeta_missing, ...
    existing_output_subjects, node_output_subjects, ...
    'VariableNames', {'Group', 'Retained', ...
    'InputSubjects', 'PTIDSubjects', 'ABetaNegative', 'ABetaPositive', ...
    'ABetaMissing', 'ExistingOutputSubjects', 'NodeOutputSubjects'});

disp(report)
fprintf('Cohort and dimension validation passed (deep=%d).\n', deep);
end
