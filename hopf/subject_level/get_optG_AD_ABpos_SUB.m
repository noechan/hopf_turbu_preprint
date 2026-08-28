clear all, clc

% Directory containing all files (no per-subject folders needed)
alldir = '/path/to/ADNI3/code/HPC_Hopf_SUB_DTI_1000_Staging/results/G_range/AD_ABpos/';

saveDir = fullfile(alldir, 'optG_AD_ABpos');

% Get list of all .mat files
allFiles = dir(fullfile(alldir, '*.mat'));
allNames = {allFiles.name};

% Extract subject identifiers (3 digits after WG_)
expr = '^WG_(\d{3})_';
tokens = regexp(allNames, expr, 'tokens');
subjectIDs = unique( cellfun(@(c) c{1}{1}, ...
                    tokens(~cellfun('isempty', tokens)), ...
                    'UniformOutput', false) );

fprintf('Found %d subjects.\n', numel(subjectIDs));

% -------------------------------------------------
% VECTOR TO COLLECT ALL SUBJECT-SPECIFIC optG
% -------------------------------------------------
optG_all_subjects = zeros(numel(subjectIDs), 1);

% -------------------------------------------------
% LOOP OVER SUBJECTS
% -------------------------------------------------
for s = 1:numel(subjectIDs)

    subj = subjectIDs{s};
    fprintf('Processing subject %s...\n', subj);

    % regex to select this subject's files
    subjExpr = ['^WG_' subj '_\d{3}_.*\.mat$'];
    subjFiles = allFiles(~cellfun('isempty', regexp(allNames, subjExpr)));

    nFiles = numel(subjFiles);
    meanErr = zeros(nFiles,1);
    Gvals   = zeros(nFiles,1);

    for i = 1:nFiles
        data = load(fullfile(alldir, subjFiles(i).name));
        meanErr(i) = mean(data.err_hete);
        Gvals(i)   = data.G;
    end

    % Sort by G and compute optimal G for this subject
    [G_sorted, idx] = sort(Gvals);
    err_sorted = meanErr(idx);

    [~, imin] = min(err_sorted);
    optG = G_sorted(imin);

    % Save subject-wise optimal G
    save(fullfile(saveDir, ['optG_subj_' subj '.mat']), ...
        'optG', 'G_sorted', 'err_sorted');

    % Store in group vector
    optG_all_subjects(s) = optG;

end

% -------------------------------------------------
% GROUP-LEVEL mean optimal G (AD_ABpos group)
% -------------------------------------------------
mean_optG_AD = mean(optG_all_subjects);
std_optG_AD  = std(optG_all_subjects);

fprintf('\nGroup-level mean optimal G (AD_ABpos): %.4f\n', mean_optG_AD);
fprintf('Standard deviation: %.4f\n', std_optG_AD);

% Save group summary
save(fullfile(saveDir, 'optG_AD_ABpos_group_stats.mat'), ...
    'optG_all_subjects', 'mean_optG_AD', 'std_optG_AD');
