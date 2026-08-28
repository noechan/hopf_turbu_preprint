function vcondition = load_condition_InfoCap_ADNI3_HC_ABneg_SUB(trials, cond, nSub, baseDir)
% Subject-level loader for HC_ABneg
% Loads files:
% results/Wtrials/HC_ABneg/Wtrials_%03d_%d_err_hete_ADNI3_HC_ABneg_GEC_sch1000_SUB%03d.mat
%
% Returns:
%   vcondition.infocap_all        [trials x nSub]
%   vcondition.suscep_all         [trials x nSub]
%   vcondition.infocap_subj_mean  [nSub x 1]
%   vcondition.suscep_subj_mean   [nSub x 1]
%   vcondition.infocap_group_mean scalar
%   vcondition.suscep_group_mean  scalar

if nargin < 4 || isempty(baseDir)
    baseDir = pwd;
end

groupLabel = 'HC_ABneg';
outDir = fullfile(baseDir, 'results', 'Wtrials', groupLabel);

if ~exist(outDir, 'dir')
    error('Folder not found: %s\nCheck baseDir and folder structure.', outDir);
end

vcondition.group = groupLabel;
vcondition.trials = trials;
vcondition.cond = cond;
vcondition.nSub = nSub;

vcondition.infocap_all = nan(trials, nSub);
vcondition.suscep_all  = nan(trials, nSub);

missing = 0;
expected = trials * nSub;

for n = 1:nSub
    for i = 1:trials
        fname = sprintf('Wtrials_%03d_%d_err_hete_ADNI3_%s_GEC_sch1000_SUB%03d.mat', ...
                        i, cond, groupLabel, n);
        fpath = fullfile(outDir, fname);

        if exist(fpath, 'file')
            S = load(fpath, 'infocapacity', 'susceptibility');
            if isfield(S, 'infocapacity')
                vcondition.infocap_all(i,n) = S.infocapacity;
            end
            if isfield(S, 'susceptibility')
                vcondition.suscep_all(i,n) = S.susceptibility;
            end
        else
            missing = missing + 1;
        end
    end
end

% Subject means across trials
vcondition.infocap_subj_mean = mean(vcondition.infocap_all, 1, 'omitnan')';
vcondition.suscep_subj_mean  = mean(vcondition.suscep_all,  1, 'omitnan')';

% Group means across subjects
vcondition.infocap_group_mean = mean(vcondition.infocap_subj_mean, 'omitnan');
vcondition.suscep_group_mean  = mean(vcondition.suscep_subj_mean,  'omitnan');

% Optional: SEM across subjects
nEffIC = sum(~isnan(vcondition.infocap_subj_mean));
nEffS  = sum(~isnan(vcondition.suscep_subj_mean));
vcondition.infocap_group_sem = std(vcondition.infocap_subj_mean, 'omitnan') / max(sqrt(nEffIC),1);
vcondition.suscep_group_sem  = std(vcondition.suscep_subj_mean,  'omitnan') / max(sqrt(nEffS),1);

% Report missing files
vcondition.missing_files = missing;
vcondition.expected_files = expected;

if missing > 0
    fprintf('[%s] Missing %d/%d files (%.2f%%). Means computed with available files.\n', ...
        groupLabel, missing, expected, 100*missing/expected);
end
end
