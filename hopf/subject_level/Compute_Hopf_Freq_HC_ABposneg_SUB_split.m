%% SUBJECT-LEVEL Whole-brain modelling in the homogeneous case (Schaefer-1000)
% Computes subject-level FC, power spectra, peak frequency (f_diff) per parcel,
% AND subject-level spatial correlation functions corrfcn_SUB.
% Also builds Aβ+ / Aβ− masks aligned to the subject order in `tseries`,
% and saves separate subject-level files for HC Aβ+ and HC Aβ−.

clear; clc;

% -------------------- Paths --------------------
cd('/path/to/ADNI3/code/HPC_Hopf_SUB_DTI_1000_Staging/');
addpath(genpath(pwd));

load('tseries_ADNI3_HC_MPRAGE_IRFSPGR_sch1000_N238rev.mat');   % -> tseries
ptid_path = '/path/to/ADNI3/code/turbulence_fulldenoising/sch400_N238rev';
load(fullfile(ptid_path,'PTID_ADNI3_HC_MPRAGE_IRFSPGR_all.mat')); % -> PTID
PTID_HC = cellstr(PTID);

abeta_table = readtable('/path/to/ADNI3/participants/cross_sectional/ADNI3_N238rev_with_ABETA_Status_CL24.xlsx');

% Structural coordinates for spatial correlation
load('SchaeferCOG.mat');      % SchaeferCOG: [NPARCELLS x 3], coordinates

% -------------------- Acquisition / processing params --------------------
TR         = 3;              % s
NPARCELLS  = 1000;
NCOND      = 1;
TT_target  = 197;            % common length for spectral estimation

% -------------------- Bandpass 0.008–0.08 Hz -----------------------------
fnq = 1/(2*TR);
flp = 0.008; fhi = 0.08;
Wn  = [flp/fnq fhi/fnq];
k   = 2;
[bfilt,afilt] = butter(k,Wn);

% -------------------- Helpers --------------------
gauss_sigma = 0.01;          % for spectral smoothing


% -------------------- Prepare time-series containers ---------------------
xs   = tseries(:,NCOND);                      % cell array per subject
sub_has_data = ~cellfun(@isempty, xs);
NSUB = nnz(sub_has_data);
fprintf('HC subjects with data: %d\n', NSUB);

% Subject order used hereafter:
PTID_in_ts = PTID_HC(sub_has_data);          % same order as xs(sub_has_data)
sub_idx    = find(sub_has_data);             % indices into original tseries

% -------------------- Frequency grid (common to all subjects) ------------
Ts        = TT_target * TR;
freq      = (0:TT_target/2-1)./Ts;           % (nfreqs x 1)
nfreqs    = numel(freq);

% -------------------- Spatial correlation settings -----------------------
NR       = 400;                               % number of distance bins
rr       = zeros(NPARCELLS,NPARCELLS);        % Euclidean distances
for i = 1:NPARCELLS
    for j = 1:NPARCELLS
        rr(i,j) = norm(SchaeferCOG(i,:) - SchaeferCOG(j,:));
    end
end
range_rr = max(rr(:));
delta    = range_rr/NR;                       % bin width

xrange   = zeros(1,NR);
for i = 1:NR
    xrange(i) = delta/2 + delta*(i-1);        % bin centres (not strictly needed here, but useful later)
end

% -------------------- Preallocate subject-level containers ---------------
FCemp_SUB    = zeros(NSUB, NPARCELLS, NPARCELLS);
PowSpect_SUB = zeros(nfreqs, NPARCELLS, NSUB);
f_diff_SUB   = zeros(NSUB, NPARCELLS);        % peak frequency per parcel
corrfcn_SUB  = zeros(NSUB, NPARCELLS, NR);    % spatial corr. per subject, parcel, distance bin

% -------------------- Loop over subjects ------------------------
parfor ii = 1:NSUB

    s  = sub_idx(ii);           % index into original tseries array
    ts = xs{s};                 % NPARCELLS x Tsub
    if isempty(ts); continue; end

    [N, Tsub] = size(ts);
    if N ~= NPARCELLS
        error('Unexpected number of parcels for subject %d (found %d, expected %d).', s, N, NPARCELLS);
    end

    % Detrend + filter each parcel; enforce common length TT_target
    tss = zeros(NPARCELLS, TT_target);

    for p = 1:NPARCELLS
        x = ts(p,1:Tsub);
        x = detrend(x - mean(x,'omitnan'));
        x(~isfinite(x)) = 0;
        xf = filtfilt(bfilt, afilt, x);

        if Tsub >= TT_target
            tss(p,:) = xf(1:TT_target);
        else
            % zero-pad to TT_target
            pad        = zeros(1,TT_target);
            pad(1:Tsub) = xf;
            tss(p,:)    = pad;
        end
    end

    % ---------------- FC (filtered) --------------------------------------
    FC_emp = corrcoef(tss(1:NPARCELLS,:)','Rows','pairwise');   % [NPARCELLS x NPARCELLS]
    FCemp_SUB(ii,:,:) = FC_emp;

    % ---------------- Power spectra + peak frequency ---------------------
    % FFT of each parcel (common TT_target for everyone)
    pw = abs(fft(tss, [], 2)).^2 / (TT_target/TR);  % N x TT_target
    pw = pw(:, 1:floor(TT_target/2));               % keep non-redundant half

    local_PowSpect = zeros(nfreqs, NPARCELLS);
    local_f_diff   = zeros(1, NPARCELLS);

    for p = 1:NPARCELLS
        sm  = gaussfilt(freq, squeeze(pw(p,:))', gauss_sigma);   % smoothed spectrum
        local_PowSpect(:,p) = sm;
        [~, idx_max] = max(sm);
        local_f_diff(p) = freq(idx_max);
    end

    PowSpect_SUB(:,:,ii) = local_PowSpect;
    f_diff_SUB(ii,:)     = local_f_diff;

    % ---------------- Spatial correlation function corrfcn_SUB ----------
    % Follows the same logic used for empirical_spacorr_rest_cond_1_*.m
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

end % parfor

% -------------------- Guard against zeros at 0 Hz ------------------------
mask_zero = (f_diff_SUB==0);
if any(mask_zero(:))
    f_diff_SUB(mask_zero) = mean(f_diff_SUB(~mask_zero),'all','omitnan');
end

% -------------------- Build Aβ masks in SUBJECT ORDER --------------------
% abeta_table: variables PTID, GROUP (e.g., 'HC'), ABeta_pvc (1=Aβ+, 0=Aβ−)
isHC  = strcmp(abeta_table.GROUP,'HC');
tabHC = abeta_table(isHC, {'PTID','ABeta_pvc'});

% Align to the tseries order used above (PTID_in_ts)
[lia, loc] = ismember(PTID_in_ts, tabHC.PTID);
if ~all(lia)
    warning('Some PTIDs in tseries were not found in Aβ table. They will be marked as NaN and excluded from Aβ+ / Aβ− groups.');
end

abeta_flag        = nan(NSUB,1);           % 1=Aβ+, 0=Aβ−, NaN=unknown
abeta_flag(lia)   = tabHC.ABeta_pvc(loc(lia));

idx_ABpos = find(abeta_flag==1);
idx_ABneg = find(abeta_flag==0);

fprintf('Aβ+ HC subjects with data: %d\n', numel(idx_ABpos));
fprintf('Aβ− HC subjects with data: %d\n', numel(idx_ABneg));

% -------------------- Group means (optional, for compatibility) ---------
FCemp_mean_all    = squeeze(mean(FCemp_SUB,1,'omitnan'));           % [N x N]
FCemp_mean_ABpos  = squeeze(mean(FCemp_SUB(idx_ABpos,:,:),1,'omitnan'));
FCemp_mean_ABneg  = squeeze(mean(FCemp_SUB(idx_ABneg,:,:),1,'omitnan'));

PowSpect_mean_all   = mean(PowSpect_SUB, 3, 'omitnan');             % [nfreq x N]
PowSpect_mean_ABpos = mean(PowSpect_SUB(:,:,idx_ABpos), 3, 'omitnan');
PowSpect_mean_ABneg = mean(PowSpect_SUB(:,:,idx_ABneg), 3, 'omitnan');

% Group-average dominant frequency from mean spectra
[~, idx_all ] = max(PowSpect_mean_all,   [], 1);  f_diff_all   = freq(idx_all);
[~, idx_pos ] = max(PowSpect_mean_ABpos, [], 1);  f_diff_ABpos = freq(idx_pos);
[~, idx_neg ] = max(PowSpect_mean_ABneg, [], 1);  f_diff_ABneg = freq(idx_neg);

% Replace any zeros
f_diff_all  (f_diff_all  ==0) = mean(f_diff_all (f_diff_all ~=0));
f_diff_ABpos(f_diff_ABpos==0) = mean(f_diff_ABpos(f_diff_ABpos~=0));
f_diff_ABneg(f_diff_ABneg==0) = mean(f_diff_ABneg(f_diff_ABneg~=0));

% Group-average spatial correlation functions
corrfcn_mean_all   = squeeze(mean(corrfcn_SUB,1,'omitnan'));           % [N x NR]
corrfcn_mean_ABpos = squeeze(mean(corrfcn_SUB(idx_ABpos,:,:),1,'omitnan'));
corrfcn_mean_ABneg = squeeze(mean(corrfcn_SUB(idx_ABneg,:,:),1,'omitnan'));

% -------------------- Aβ-specific SUBJECT-LEVEL subsets -----------------
PTID_ABpos          = PTID_in_ts(idx_ABpos);
PTID_ABneg          = PTID_in_ts(idx_ABneg);

f_diff_SUB_ABpos    = f_diff_SUB(idx_ABpos,:);
f_diff_SUB_ABneg    = f_diff_SUB(idx_ABneg,:);

corrfcn_SUB_ABpos   = corrfcn_SUB(idx_ABpos,:,:);   % [Npos x N x NR]
corrfcn_SUB_ABneg   = corrfcn_SUB(idx_ABneg,:,:);   % [Nneg x N x NR]

FCemp_SUB_ABpos     = FCemp_SUB(idx_ABpos,:,:);
FCemp_SUB_ABneg     = FCemp_SUB(idx_ABneg,:,:);

% -------------------- SAVE: subject-level + group-level ------------------
out_base = '/path/to/ADNI3/code/HPC_Hopf_SUB_DTI_1000_Staging';

% (A) Subject-level arrays for ALL HC (including Aβ masks)
save(fullfile(out_base,'results_f_diff_fce_cond1_ADNI3_HC_sch1000_SUB.mat'), ...
    'freq','xrange', ...
    'PTID_in_ts','idx_ABpos','idx_ABneg','abeta_flag', ...
    'FCemp_SUB','PowSpect_SUB','f_diff_SUB','corrfcn_SUB');

% (B) Subject-level arrays for HC Aβ+ only
save(fullfile(out_base,'results_f_diff_fce_cond1_ADNI3_HC_ABpos_sch1000_SUB.mat'), ...
    'freq','xrange', ...
    'PTID_ABpos','idx_ABpos', ...
    'FCemp_SUB_ABpos','f_diff_SUB_ABpos','corrfcn_SUB_ABpos');

% (C) Subject-level arrays for HC Aβ− only
save(fullfile(out_base,'results_f_diff_fce_cond1_ADNI3_HC_ABneg_sch1000_SUB.mat'), ...
    'freq','xrange', ...
    'PTID_ABneg','idx_ABneg', ...
    'FCemp_SUB_ABneg','f_diff_SUB_ABneg','corrfcn_SUB_ABneg');

% (D) Group-level files (for compatibility with existing scripts)
save(fullfile(out_base,'results_f_diff_fce_cond1_ADNI3_HC_sch1000.mat'), ...
    'f_diff_all','FCemp_mean_all');

save(fullfile(out_base,'results_f_diff_fce_cond1_ADNI3_HC_ABpos_sch1000.mat'), ...
    'f_diff_ABpos','FCemp_mean_ABpos');

save(fullfile(out_base,'results_f_diff_fce_cond1_ADNI3_HC_ABneg_sch1000.mat'), ...
    'f_diff_ABneg','FCemp_mean_ABneg');

fprintf('\nSaved subject-level arrays for ALL HC, Aβ+, and Aβ−.\n');
fprintf('Group-level f_diff and FC/corrfcn saved with the same filenames as before.\n');
