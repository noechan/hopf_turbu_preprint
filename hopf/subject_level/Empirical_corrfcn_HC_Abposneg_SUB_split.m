clear; clc;

% -------------------------------------------------------------------------
% Paths
% -------------------------------------------------------------------------
base_dir = '/path/to/ADNI3/code/HPC_Hopf_SUB_DTI_1000_Staging/';
cd(base_dir);
addpath(genpath(base_dir));

% -------------------------------------------------------------------------
% Load:
%   - Subject-level FC (from Script 1)
%   - Schaefer parcel centroids (for Euclidean distances)
% -------------------------------------------------------------------------
sub_level_file = fullfile(base_dir, ...
    'results_f_diff_fce_cond1_ADNI3_HC_sch1000_SUB.mat');

load(sub_level_file, 'FCemp_SUB', 'PTID_in_ts', 'idx_ABpos', 'idx_ABneg');
load('SchaeferCOG.mat');  % SchaeferCOG: [NPARCELLS x 3]

% -------------------------------------------------------------------------
% Basic parameters
% -------------------------------------------------------------------------
NPARCELLS = 1000;
NR        = 400;    % number of spatial bins
NCOND     = 1;      % here only one (rest) condition

NSUB = size(FCemp_SUB, 1);
fprintf('Number of HC subjects in FCemp_SUB: %d\n', NSUB);

% -------------------------------------------------------------------------
% Precompute distance matrix and radial bins
% -------------------------------------------------------------------------
rr = zeros(NPARCELLS, NPARCELLS, 'single');
for i = 1:NPARCELLS
    for j = 1:NPARCELLS
        rr(i,j) = norm(SchaeferCOG(i,:) - SchaeferCOG(j,:));
    end
end

range_dist = max(rr(:));
delta      = range_dist / NR;                       % bin width
xrange     = (delta/2) + delta * (0:NR-1);          % bin centers (1 x NR)

% -------------------------------------------------------------------------
% Allocate subject-level container for corrfcn
% CorrFcn_SUB: [NSUB x NPARCELLS x NR]
% -------------------------------------------------------------------------
CorrFcn_SUB = zeros(NSUB, NPARCELLS, NR);

% -------------------------------------------------------------------------
% Loop over conditions (kept for naming compatibility; NCOND = 1 here)
% -------------------------------------------------------------------------
for cond = 1:NCOND

    fprintf('Computing corrfcn for condition %d...\n', cond);

    % ------------------ SUBJECT-LEVEL LOOP ------------------------------
    parfor s = 1:NSUB
        FCs = squeeze(FCemp_SUB(s,:,:));  % [NPARCELLS x NPARCELLS]

        % Local container for this subject
        corrfcn_sub_i = zeros(NPARCELLS, NR);

        for ii = 1:NPARCELLS
            numind      = zeros(1, NR);  % number of pairs per bin
            corrfcn_1   = zeros(1, NR);  % accumulated correlations

            for jj = 1:NPARCELLS
                r = rr(ii,jj);

                % Bin index
                index = floor(r/delta) + 1;
                if index > NR
                    index = NR;
                end

                mcc = FCs(ii,jj);
                if ~isnan(mcc)
                    corrfcn_1(index) = corrfcn_1(index) + mcc;
                    numind(index)    = numind(index)    + 1;
                end
            end

            % Average over all pairs in each bin
            valid_bins = numind > 0;
            corrfcn_tmp        = nan(1,NR);
            corrfcn_tmp(valid_bins) = corrfcn_1(valid_bins)./numind(valid_bins);

            corrfcn_sub_i(ii,:) = corrfcn_tmp;
        end

        CorrFcn_SUB(s,:,:) = corrfcn_sub_i;
    end

    % ---------------------------------------------------------------------
    % Aβ status is already aligned to FCemp_SUB via idx_ABpos / idx_ABneg
    %   - idx_ABpos / idx_ABneg: indices in subject-space (1..NSUB)
    %   - PTID_in_ts: subject order used in FCemp_SUB and CorrFcn_SUB
    % ---------------------------------------------------------------------
    idx_ABpos = idx_ABpos(:);
    idx_ABneg = idx_ABneg(:);

    % Aβ-specific PTID lists
    PTID_ABpos = PTID_in_ts(idx_ABpos);
    PTID_ABneg = PTID_in_ts(idx_ABneg);

    % Subject-level containers per group
    CorrFcn_SUB_ABpos = CorrFcn_SUB(idx_ABpos,:,:);   % [Npos x NPARCELLS x NR]
    CorrFcn_SUB_ABneg = CorrFcn_SUB(idx_ABneg,:,:);   % [Nneg x NPARCELLS x NR]

    % ---------------------------------------------------------------------
    % GROUP-LEVEL MEANS (over subjects)
    % ---------------------------------------------------------------------
    corrfcn_mean_all   = squeeze(nanmean(CorrFcn_SUB,        1));  % [NPARCELLS x NR]
    corrfcn_mean_ABpos = squeeze(nanmean(CorrFcn_SUB_ABpos,  1));  % [NPARCELLS x NR]
    corrfcn_mean_ABneg = squeeze(nanmean(CorrFcn_SUB_ABneg,  1));  % [NPARCELLS x NR]

    % For backward compatibility with previous variable names
    corrfcn          = corrfcn_mean_all;
    corrfcn_HC_ABpos = corrfcn_mean_ABpos;
    corrfcn_HC_ABneg = corrfcn_mean_ABneg;

    % ---------------------------------------------------------------------
    % SAVE: subject-level + group-level, by Aβ status
    % ---------------------------------------------------------------------
    out_base = base_dir;

    % (A) Subject-level spatial correlations for ALL HC
    save(fullfile(out_base, ...
        sprintf('empirical_spacorr_rest_cond_%d_ADNI3_HC_sch1000_SUB.mat', cond)), ...
        'xrange', 'rr', ...
        'PTID_in_ts', 'idx_ABpos', 'idx_ABneg', ...
        'CorrFcn_SUB');

    % (B) Subject-level spatial correlations for HC Aβ+
    save(fullfile(out_base, ...
        sprintf('empirical_spacorr_rest_cond_%d_ADNI3_HC_ABpos_sch1000_SUB.mat', cond)), ...
        'xrange', 'rr', ...
        'PTID_ABpos', 'idx_ABpos', ...
        'CorrFcn_SUB_ABpos');

    % (C) Subject-level spatial correlations for HC Aβ−
    save(fullfile(out_base, ...
        sprintf('empirical_spacorr_rest_cond_%d_ADNI3_HC_ABneg_sch1000_SUB.mat', cond)), ...
        'xrange', 'rr', ...
        'PTID_ABneg', 'idx_ABneg', ...
        'CorrFcn_SUB_ABneg');

    % (D) Group-level averages (as before)
    save(fullfile(out_base, ...
        sprintf('empirical_spacorr_rest_cond_%d_ADNI3_HC_sch1000.mat', cond)), ...
        'xrange', 'rr', 'corrfcn');

    save(fullfile(out_base, ...
        sprintf('empirical_spacorr_rest_cond_%d_ADNI3_HC_ABpos_sch1000.mat', cond)), ...
        'xrange', 'rr', 'corrfcn_HC_ABpos');

    save(fullfile(out_base, ...
        sprintf('empirical_spacorr_rest_cond_%d_ADNI3_HC_ABneg_sch1000.mat', cond)), ...
        'xrange', 'rr', 'corrfcn_HC_ABneg');

    fprintf('\nSaved subject-level CorrFcn_SUB (all, Aβ+, Aβ−) and group-level corrfcn to:\n');
    fprintf('  empirical_spacorr_rest_cond_%d_ADNI3_HC*_sch1000[_SUB].mat\n', cond);
end
