%% SUBJECT-LEVEL spatial correlation function (corrfcn) for MCI (Schaefer-1000)

clear; clc;

base_dir = '/path/to/ADNI3/code/HPC_Hopf_SUB_DTI_1000_Staging/';
cd(base_dir);
addpath(genpath(base_dir));

% -------------------------------------------------------------------------
% Load subject-level FC and Aβ indices (from the MCI subject-level script)
% -------------------------------------------------------------------------
sub_level_file = fullfile(base_dir, ...
    'results_f_diff_fce_cond1_ADNI3_MCI_sch1000_SUB.mat');

% Expected variables in this file:
%   FCemp_SUB   : [NSUB x NPARCELLS x NPARCELLS]
%   PTID_in_ts  : {NSUB x 1} or cell array, subject order
%   idx_ABpos   : indices of Aβ+ subjects in 1..NSUB
%   idx_ABneg   : indices of Aβ− subjects in 1..NSUB
%   abeta_flag  : 1=Aβ+, 0=Aβ−, NaN=unknown
load(sub_level_file, 'FCemp_SUB', 'PTID_in_ts', 'idx_ABpos', 'idx_ABneg', 'abeta_flag');

% Structural coordinates for Schaefer-1000 parcels
load('SchaeferCOG.mat');    % SchaeferCOG: [NPARCELLS x 3]

% -------------------------------------------------------------------------
% Data / analysis parameters
% -------------------------------------------------------------------------
NPARCELLS = 1000;
NR        = 400;        % number of spatial bins
NCOND     = 1;          % here only one (rest) condition

NSUB = size(FCemp_SUB, 1);
fprintf('MCI subjects in FCemp_SUB: %d\n', NSUB);

% -------------------------------------------------------------------------
% Precompute Euclidean distances between parcels and radial bins
% -------------------------------------------------------------------------
rr = zeros(NPARCELLS, NPARCELLS);
for i = 1:NPARCELLS
    for j = 1:NPARCELLS
        rr(i,j) = norm(SchaeferCOG(i,:) - SchaeferCOG(j,:));
    end
end
range_dist = max(rr(:));
delta      = range_dist / NR;

xrange = zeros(1, NR);
for i = 1:NR
    xrange(i) = delta/2 + delta*(i-1);
end

% -------------------------------------------------------------------------
% Allocate subject-level container for corrfcn
% CorrFcn_SUB: [NSUB x NPARCELLS x NR]
% -------------------------------------------------------------------------
CorrFcn_SUB = zeros(NSUB, NPARCELLS, NR);

% -------------------------------------------------------------------------
% Loop over conditions (only cond=1, kept for naming compatibility)
% -------------------------------------------------------------------------
for cond = 1:NCOND

    fprintf('Computing corrfcn for MCI, condition %d...\n', cond);

    % ---------------------------------------------------------------------
    % SUBJECT-LEVEL LOOP (use FCemp_SUB directly)
    % ---------------------------------------------------------------------
    parfor s = 1:NSUB
        FCs = squeeze(FCemp_SUB(s,:,:));  % [NPARCELLS x NPARCELLS]

        % Local container for this subject
        corrfcn_sub_i = zeros(NPARCELLS, NR);

        for ii = 1:NPARCELLS
            numind      = zeros(1, NR);
            corrfcn_1   = zeros(1, NR);

            for jj = 1:NPARCELLS
                r = rr(ii,jj);
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

            % Average over bins; leave NaN where no pairs existed
            valid_bins = numind > 0;
            tmp        = nan(1,NR);
            tmp(valid_bins) = corrfcn_1(valid_bins)./numind(valid_bins);

            corrfcn_sub_i(ii,:) = tmp;
        end

        CorrFcn_SUB(s,:,:) = corrfcn_sub_i;
    end

    % ---------------------------------------------------------------------
    % Split subject-level corrfcn by Aβ status (using precomputed indices)
    % ---------------------------------------------------------------------
    idx_ABpos = idx_ABpos(:);
    idx_ABneg = idx_ABneg(:);

    PTID_ABpos = PTID_in_ts(idx_ABpos);
    PTID_ABneg = PTID_in_ts(idx_ABneg);

    CorrFcn_SUB_ABpos = CorrFcn_SUB(idx_ABpos,:,:);   % [Npos x NPARCELLS x NR]
    CorrFcn_SUB_ABneg = CorrFcn_SUB(idx_ABneg,:,:);   % [Nneg x NPARCELLS x NR]

    % Backward-compatible variable names (optional but convenient)
    corrfcnsub             = CorrFcn_SUB;
    CorrFcn_sub_MCI_ABpos  = CorrFcn_SUB_ABpos;
    CorrFcn_sub_MCI_ABneg  = CorrFcn_SUB_ABneg;

    % ---------------------------------------------------------------------
    % GROUP-LEVEL MEANS
    % ---------------------------------------------------------------------
    corrfcn           = squeeze(nanmean(CorrFcn_SUB,        1));  % [NPARCELLS x NR]
    corrfcn_MCI_ABpos = squeeze(nanmean(CorrFcn_SUB_ABpos,  1));
    corrfcn_MCI_ABneg = squeeze(nanmean(CorrFcn_SUB_ABneg,  1));

    % ---------------------------------------------------------------------
    % SAVE: subject-level + group-level
    % ---------------------------------------------------------------------
    cd(base_dir);

    % (A) Subject-level (ALL MCI)
    save(sprintf('empirical_spacorr_rest_cond_%d_ADNI3_MCI_sch1000_SUB.mat', cond), ...
        'xrange', 'rr', ...
        'PTID_in_ts', 'abeta_flag', 'idx_ABpos', 'idx_ABneg', ...
        'CorrFcn_SUB', ...                 % new name
        'corrfcnsub', ...                  % old name for compatibility
        'CorrFcn_sub_MCI_ABpos', 'CorrFcn_sub_MCI_ABneg', ...
        '-v7.3');

    % (B) Subject-level for MCI Aβ+
    save(sprintf('empirical_spacorr_rest_cond_%d_ADNI3_MCI_ABpos_sch1000_SUB.mat', cond), ...
        'xrange', 'rr', ...
        'PTID_ABpos', 'idx_ABpos', ...
        'CorrFcn_SUB_ABpos', 'CorrFcn_sub_MCI_ABpos', ...
        '-v7.3');

    % (C) Subject-level for MCI Aβ−
    save(sprintf('empirical_spacorr_rest_cond_%d_ADNI3_MCI_ABneg_sch1000_SUB.mat', cond), ...
        'xrange', 'rr', ...
        'PTID_ABneg', 'idx_ABneg', ...
        'CorrFcn_SUB_ABneg', 'CorrFcn_sub_MCI_ABneg', ...
        '-v7.3');

    % (D) Group-level (for compatibility with downstream scripts)
    save(sprintf('empirical_spacorr_rest_cond_%d_ADNI3_MCI_sch1000.mat', cond), ...
        'xrange', 'rr', 'corrfcn');

    save(sprintf('empirical_spacorr_rest_cond_%d_ADNI3_MCI_ABpos_sch1000.mat', cond), ...
        'xrange', 'rr', 'corrfcn_MCI_ABpos');

    save(sprintf('empirical_spacorr_rest_cond_%d_ADNI3_MCI_ABneg_sch1000.mat', cond), ...
        'xrange', 'rr', 'corrfcn_MCI_ABneg');

    fprintf('\nSaved MCI subject-level and group-level corrfcn files for cond %d.\n', cond);
end
