function G = aggregate_subject_level_group(groupName, baseDir, nSub, trials, condition, measure)
% aggregate_subject_level_group
%
% Deterministic filename-template implementation for grouped folders.
%
% Expected structure:
%   <baseDir>/<groupName>/Wtrials_<trial>_<condition>_<measure>_ADNI3_<groupName>_GEC_sch1000_SUB<subject>.mat

IC_subject_means = nan(nSub,1);
S_subject_means  = nan(nSub,1);

IC_subject_trials = cell(nSub,1);
S_subject_trials  = cell(nSub,1);

nTrials_found = zeros(nSub,1);

groupDir = fullfile(baseDir, groupName);

if ~exist(groupDir, 'dir')
    error('Group directory does not exist: %s', groupDir);
end

for subIdx = 1:nSub

    IC_trials = nan(trials,1);
    S_trials  = nan(trials,1);
    validTrial = false(trials,1);

    for trialIdx = 1:trials

        fname = sprintf( ...
            'Wtrials_%03d_%d_%s_ADNI3_%s_GEC_sch1000_SUB%03d.mat', ...
            trialIdx, condition, measure, groupName, subIdx);

        fpath = fullfile(groupDir, fname);

        if ~exist(fpath, 'file')
            warning('Missing file: %s', fpath);
            continue
        end

        D = load(fpath);

        if ~isfield(D, 'infocapacity') || ~isfield(D, 'susceptibility')
            warning('Missing required variables in file: %s', fpath);
            continue
        end

        IC_trials(trialIdx) = D.infocapacity;
        S_trials(trialIdx)  = D.susceptibility;
        validTrial(trialIdx) = true;
    end

    IC_trials = IC_trials(validTrial);
    S_trials  = S_trials(validTrial);

    nTrials_found(subIdx) = sum(validTrial);

    IC_subject_trials{subIdx} = IC_trials;
    S_subject_trials{subIdx}  = S_trials;

    if isempty(IC_trials)
        warning('Subject %03d in group %s has no valid trials.', subIdx, groupName);
        continue
    end

    IC_subject_means(subIdx) = mean(IC_trials, 'omitnan');
    S_subject_means(subIdx)  = mean(S_trials, 'omitnan');
end

validIC = ~isnan(IC_subject_means);
validS  = ~isnan(S_subject_means);

nValidIC = sum(validIC);
nValidS  = sum(validS);

IC_group_mean = mean(IC_subject_means(validIC), 'omitnan');
S_group_mean  = mean(S_subject_means(validS), 'omitnan');

if nValidIC > 1
    IC_group_sem = std(IC_subject_means(validIC), 'omitnan') / sqrt(nValidIC);
else
    IC_group_sem = NaN;
end

if nValidS > 1
    S_group_sem = std(S_subject_means(validS), 'omitnan') / sqrt(nValidS);
else
    S_group_sem = NaN;
end

G = struct;
G.groupName = groupName;
G.groupDir  = groupDir;

G.IC_subject_trials = IC_subject_trials;
G.S_subject_trials  = S_subject_trials;

G.IC_subject_means = IC_subject_means;
G.S_subject_means  = S_subject_means;

G.IC_group_mean = IC_group_mean;
G.IC_group_sem  = IC_group_sem;

G.S_group_mean = S_group_mean;
G.S_group_sem  = S_group_sem;

G.nTrials_found = nTrials_found;

G.nValidSubjects_IC = nValidIC;
G.nValidSubjects_S  = nValidS;
end