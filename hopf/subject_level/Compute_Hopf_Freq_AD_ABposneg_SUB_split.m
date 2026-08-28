%% SUBJECT-LEVEL Whole-brain modelling (AD, Schaefer-1000)
% Computes subject-level FC, power spectra, parcel-wise dominant frequency,
% and subject-level spatial correlation functions corrfcn_SUB.
% Also builds Aβ+ / Aβ− masks aligned to the subject order in `tseries`.

clear; clc;

% -------------------- Paths --------------------
cd('/path/to/ADNI3/code/HPC_Hopf_SUB_DTI_1000_Staging/');
addpath(genpath(pwd));

load('tseries_ADNI3_AD_MPRAGE_IRFSPGR_sch1000_N238rev.mat');   % -> tseries {sub,cond}

ptid_path = '/path/to/ADNI3/code/turbulence_fulldenoising/sch400_N238rev';
load(fullfile(ptid_path,'PTID_ADNI3_AD_MPRAGE_IRFSPGR_all.mat'));  % -> PTID (same order as tseries)
PTID_AD = cellstr(PTID);

abeta_table = readtable('/path/to/ADNI3/participants/cross_sectional/ADNI3_N238rev_with_ABETA_Status_CL24.xlsx');

% Structural coordinates for spatial correlation (Schaefer-1000)
load('SchaeferCOG.mat');   % SchaeferCOG: [NPARCELLS x 3]

% -------------------- Parameters --------------------
TR         = 3;         % s
NPARCELLS  = 1000;
NCOND      = 1;
TT_target  = 197;       % common length for spectra

% -------------------- Bandpass 0.008–0.08 Hz -----------------------------
fnq = 1/(2*TR); 
flp = 0.008; 
fhi = 0.08; 
Wn  = [flp/fnq fhi/fnq];
[bfilt,afilt] = butter(2, Wn);

% -------------------- Shared frequency grid --------------------------------
Ts      = TT_target * TR;
freq    = (0:TT_target/2-1)./Ts;     % nfreqs x 1
nfreqs  = numel(freq);
gauss_sigma = 0.01;

% -------------------- Spatial correlation settings -----------------------
NR       = 400;                          % number of distance bins
rr       = zeros(NPARCELLS,NPARCELLS);   % Euclidean distances
for i = 1:NPARCELLS
    for j = 1:NPARCELLS
        rr(i,j) = norm(SchaeferCOG(i,:) - SchaeferCOG(j,:));
    end
end
range_rr = max(rr(:));
delta    = range_rr/NR;                  % bin width

xrange   = zeros(1,NR);                  % distance-bin centres
for i = 1:NR
    xrange(i) = delta/2 + delta*(i-1);
end

% -------------------- Containers -------------------------------------------
xs   = tseries(:,NCOND);
sub_has_data = ~cellfun(@isempty, xs);
NSUB = nnz(sub_has_data);
fprintf('AD subjects with data: %d\n', NSUB);

FCemp_SUB    = zeros(NSUB, NPARCELLS, NPARCELLS);
PowSpect_SUB = zeros(nfreqs, NPARCELLS, NSUB);
f_diff_SUB   = zeros(NSUB, NPARCELLS);
corrfcn_SUB  = zeros(NSUB, NPARCELLS, NR);   % subject-level spatial corr

PTID_in_ts   = PTID_AD(sub_has_data);      % subject order used here
sub_idx      = find(sub_has_data);

% -------------------- Subject loop -----------------------------------------
parfor ii = 1:NSUB
    s  = sub_idx(ii);
    ts = xs{s};                             % N x Tsub
    if isempty(ts); continue; end
    [N, Tsub] = size(ts);
    if N ~= NPARCELLS
        error('Unexpected number of parcels for subject %d (found %d, expected %d).', ...
              s, N, NPARCELLS);
    end

    % Detrend + bandpass; enforce common length
    tss = zeros(NPARCELLS, TT_target);
    for p = 1:NPARCELLS
        x = ts(p,1:Tsub);
        x = detrend(x - mean(x,'omitnan'));
        x(~isfinite(x)) = 0;
        xf = filtfilt(bfilt, afilt, x);
        if Tsub >= TT_target
            tss(p,:) = xf(1:TT_target);
        else
            pad = zeros(1,TT_target); 
            pad(1:Tsub) = xf;
            tss(p,:) = pad;
        end
    end

    % Subject FC (filtered)
    FC_emp = corrcoef(tss(1:NPARCELLS,:)','Rows','pairwise');  % N x N
    FCemp_SUB(ii,:,:) = FC_emp;

    % Spectra and dominant frequency per parcel
    pw = abs(fft(tss, [], 2)).^2 / (TT_target/TR);   % N x TT_target
    pw = pw(:, 1:floor(TT_target/2));                % keep non-redundant half
    for p = 1:NPARCELLS
        sm = gaussfilt(freq, squeeze(pw(p,:))', gauss_sigma);
        PowSpect_SUB(:,p,ii) = sm;
        [~, idx_peak] = max(sm);
        f_diff_SUB(ii,p) = freq(idx_peak);
    end

    % Subject-level spatial correlation corrfcn_SUB
    local_corrfcn = zeros(NPARCELLS, NR);
    for i = 1:NPARCELLS
        numind    = zeros(1,NR);
        corrfcn_1 = zeros(1,NR);
        for j = 1:NPARCELLS
            r     = rr(i,j);
            index = floor(r/delta) + 1;
            if index == NR+1
                index = NR;
            end
            mcc = FC_emp(i,j);
            if ~isnan(mcc)
                corrfcn_1(index) = corrfcn_1(index) + mcc;
                numind(index)    = numind(index) + 1;
            end
        end
        corrfcn_1(numind>0) = corrfcn_1(numind>0)./numind(numind>0);
        local_corrfcn(i,:)  = corrfcn_1;
    end

    corrfcn_SUB(ii,:,:) = local_corrfcn;
end

% Replace any zeros defensively
mask_zero = (f_diff_SUB==0);
if any(mask_zero(:))
    f_diff_SUB(mask_zero) = mean(f_diff_SUB(~mask_zero),'all','omitnan');
end

% -------------------- Aβ masks aligned to PTID_in_ts -----------------------
isAD  = strcmp(abeta_table.GROUP,'AD');
tabAD = abeta_table(isAD, {'PTID','ABeta_pvc'});

[lia, loc] = ismember(PTID_in_ts, tabAD.PTID);
if ~all(lia)
    warning('Some PTIDs in tseries were not found in ABeta table. Marking them as NaN.');
end
abeta_flag      = nan(NSUB,1);                    % 1=Aβ+, 0=Aβ−, NaN=unknown
abeta_flag(lia) = tabAD.ABeta_pvc(loc(lia));      % 1=Aβ+, 0=Aβ−

idx_ABpos = find(abeta_flag==1);
idx_ABneg = find(abeta_flag==0);

fprintf('AD Aβ+ subjects with data: %d\n', numel(idx_ABpos));
fprintf('AD Aβ− subjects with data: %d\n', numel(idx_ABneg));

% -------------------- Group summaries from subject-level arrays ------------
FCemp_mean_all     = squeeze(mean(FCemp_SUB,1,'omitnan'));
FCemp_mean_ABpos   = squeeze(mean(FCemp_SUB(idx_ABpos,:,:),1,'omitnan'));
FCemp_mean_ABneg   = squeeze(mean(FCemp_SUB(idx_ABneg,:,:),1,'omitnan'));

PowSpect_mean_all   = mean(PowSpect_SUB, 3, 'omitnan');                  % nfreqs x N
PowSpect_mean_ABpos = mean(PowSpect_SUB(:,:,idx_ABpos), 3, 'omitnan');
PowSpect_mean_ABneg = mean(PowSpect_SUB(:,:,idx_ABneg), 3, 'omitnan');

[~, idx_all ] = max(PowSpect_mean_all,   [], 1);  f_diff_all      = freq(idx_all);
[~, idx_pos ] = max(PowSpect_mean_ABpos, [], 1);  f_diff_AD_ABpos = freq(idx_pos);
[~, idx_neg ] = max(PowSpect_mean_ABneg, [], 1);  f_diff_AD_ABneg = freq(idx_neg);

% Defensive zero handling
f_diff_all( f_diff_all==0 )         = mean(f_diff_all(f_diff_all~=0));
f_diff_AD_ABpos(f_diff_AD_ABpos==0) = mean(f_diff_AD_ABpos(f_diff_AD_ABpos~=0));
f_diff_AD_ABneg(f_diff_AD_ABneg==0) = mean(f_diff_AD_ABneg(f_diff_AD_ABneg~=0));

% Group-average spatial correlation functions
corrfcn_mean_all   = squeeze(mean(corrfcn_SUB,1,'omitnan'));              % [N x NR]
corrfcn_mean_ABpos = squeeze(mean(corrfcn_SUB(idx_ABpos,:,:),1,'omitnan'));
corrfcn_mean_ABneg = squeeze(mean(corrfcn_SUB(idx_ABneg,:,:),1,'omitnan'));

% -------------------- Aβ-specific SUBJECT-LEVEL subsets ------------------
PTID_ABpos          = PTID_in_ts(idx_ABpos);
PTID_ABneg          = PTID_in_ts(idx_ABneg);

f_diff_SUB_ABpos    = f_diff_SUB(idx_ABpos,:);
f_diff_SUB_ABneg    = f_diff_SUB(idx_ABneg,:);

corrfcn_SUB_ABpos   = corrfcn_SUB(idx_ABpos,:,:);   % [Npos x N x NR]
corrfcn_SUB_ABneg   = corrfcn_SUB(idx_ABneg,:,:);   % [Nneg x N x NR]

FCemp_SUB_ABpos     = FCemp_SUB(idx_ABpos,:,:);
FCemp_SUB_ABneg     = FCemp_SUB(idx_ABneg,:,:);

% -------------------- Save outputs ----------------------------------------
out_base = '/path/to/ADNI3/code/HPC_Hopf_SUB_DTI_1000_Staging';

% (A) Subject-level containers for ALL AD
save(fullfile(out_base,'results_f_diff_fce_cond1_ADNI3_AD_sch1000_SUB.mat'), ...
    'freq','xrange', ...
    'PTID_in_ts','idx_ABpos','idx_ABneg','abeta_flag', ...
    'FCemp_SUB','PowSpect_SUB','f_diff_SUB','corrfcn_SUB', ...
    '-v7.3');

% (B) Subject-level arrays for AD Aβ+ only
save(fullfile(out_base,'results_f_diff_fce_cond1_ADNI3_AD_ABpos_sch1000_SUB.mat'), ...
    'freq','xrange', ...
    'PTID_ABpos','idx_ABpos', ...
    'FCemp_SUB_ABpos','f_diff_SUB_ABpos','corrfcn_SUB_ABpos', ...
    '-v7.3');

% (C) Subject-level arrays for AD Aβ− only
save(fullfile(out_base,'results_f_diff_fce_cond1_ADNI3_AD_ABneg_sch1000_SUB.mat'), ...
    'freq','xrange', ...
    'PTID_ABneg','idx_ABneg', ...
    'FCemp_SUB_ABneg','f_diff_SUB_ABneg','corrfcn_SUB_ABneg', ...
    '-v7.3');

% (D) Backward-compatible group files (f_diff + mean FC)
save(fullfile(out_base,'results_f_diff_fce_cond1_ADNI3_AD_sch1000.mat'), ...
    'f_diff_all','FCemp_mean_all');

save(fullfile(out_base,'results_f_diff_fce_cond1_ADNI3_AD_ABpos_sch1000.mat'), ...
    'f_diff_AD_ABpos','FCemp_mean_ABpos');

save(fullfile(out_base,'results_f_diff_fce_cond1_ADNI3_AD_ABneg_sch1000.mat'), ...
    'f_diff_AD_ABneg','FCemp_mean_ABneg');

