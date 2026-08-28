function pert_infocapacity_susc_MCI_ABpos_GEC_SUB(nsub, s, condition)
% SUBJECT-LEVEL infocapacity/susceptibility in the Hopf model (MCI Aβ+, Schaefer-1000)
%
% nsub      : subject index within MCI Aβ+ subgroup (1..Nneg)
% s         : random seed index (integer, used for rng and filenames)
% condition : kept for compatibility (typically 1)
%
% Uses:
%   - f_diff_SUB_ABpos(nsub,:)            : subject-specific dominant frequency
%   - rr, xrange                          : distance bins (from *_SUB* spacorr file)
%   - Ceff from GEC_SUB_%03d_ADNI3_MCI_ABpos_GEC.mat : subject-specific EC
%   - optG_MCI_ABpos(nsub)                 : subject-specific optimal global coupling
%
% Saves:
%   Wtrials_s_condition_err_hete_ADNI3_MCI_ABpos_GEC_sch1000_SUB%03d.mat
%

%% Paths and basic setup
this_file  = mfilename('fullpath');
code_dir   = fileparts(this_file);
addpath(genpath(code_dir))

out_dir = fullfile(code_dir, 'results','Wtrials','MCI_ABpos');
if ~exist(out_dir,'dir')
    mkdir(out_dir);
end

measure = 'err_hete';  % kept for compatibility with original naming

%% Load subject-level f_diff (MCI Aβ+)
freq_file = fullfile(code_dir,'data','results_f_diff_fce_cond1_ADNI3_MCI_ABpos_sch1000_SUB.mat');
% Expected: f_diff_SUB_ABpos [Nneg x NPARCELLS]
load(freq_file, 'f_diff_SUB_ABpos');

[N_SUB_ABpos, NPARCELLS] = size(f_diff_SUB_ABpos);
if nsub < 1 || nsub > N_SUB_ABpos
    error('nsub=%d is out of range for MCI Aβ+ subjects (1..%d).', nsub, N_SUB_ABpos);
end

f_sub = squeeze(f_diff_SUB_ABpos(nsub,:));  % [1 x NPARCELLS] or [NPARCELLS x 1]

%% Load distance matrix rr and xrange (from subject-level spacorr file)
spacorr_file = fullfile(code_dir,'data','empirical_spacorr_rest_cond_1_ADNI3_MCI_ABpos_sch1000_SUB.mat');
% Expected: rr [NPARCELLS x NPARCELLS], xrange [1 x NR]
load(spacorr_file, 'rr', 'xrange');

NPARCELLS = size(rr,1);
NR        = numel(xrange); % kept for completeness

%% Load subject-specific effective connectivity
gec_file = fullfile(code_dir,'data','GEC_SUB_ADNI3_MCI_ABpos', ...
                    sprintf('GEC_SUB_%03d_ADNI3_MCI_ABpos_GEC.mat', nsub));
load(gec_file, 'Ceff');   % [NPARCELLS x NPARCELLS]

%% Global coupling G (subject-specific optimal G)
if condition == 1
    optG_file = fullfile(code_dir,'results','G_range','MCI_ABpos','optG_MCI_ABpos', ...
                         'optG_MCI_ABpos_group_stats.mat');
    % This file contains: optG_all_subjects, mean_optG_MCI, std_optG_MCI
    load(optG_file, 'optG_all_subjects');

    % rename to what this function expects
    optG_MCI_ABpos = optG_all_subjects(:);   % ensure column vector

    if numel(optG_MCI_ABpos) < nsub
        error('optG_MCI_ABpos has length %d, but nsub=%d.', numel(optG_MCI_ABpos), nsub);
    end

    G = optG_MCI_ABpos(nsub);
else
    error('Condition %d not implemented for subject-level info capacity.', condition);
end


%% Geometric kernel C1 (for enstrophy)
lambda  = 0.175;
lambda  = round(lambda,2);

LAMBDA   = [0.27 0.24 0.21 0.18 0.15 0.12 0.09 0.06 0.03 0.01];
NLAMBDA  = length(LAMBDA);
C1       = zeros(NLAMBDA,NPARCELLS,NPARCELLS);

[~, indsca] = min(abs(LAMBDA - lambda));

for ilam = 1:NLAMBDA
    lambda2 = LAMBDA(ilam);
    for i = 1:NPARCELLS
        for j = 1:NPARCELLS
            C1(ilam,i,j) = exp(-lambda2 * rr(i,j));
        end
    end
end

%% Data & filter parameters
TR  = 3;  % seconds
fnq = 1/(2*TR);
flp = 0.008;
fhi = 0.08;
Wn  = [flp/fnq fhi/fnq];
k   = 2;
[bfilt,afilt] = butter(k,Wn);

%% Hopf model parameters
Tmax  = 197;
omega = repmat(2*pi*f_sub(:),1,2); 
omega(:,1) = -omega(:,1);

dt   = 0.1*TR/2;
sig  = 0.01;
dsig = sqrt(dt)*sig;

NSUBSIM = 100;   % number of stochastic realizations

%% Connectivity scaling (subject-specific EC)
factor = max(Ceff(:));
C      = Ceff / factor * 0.2;

%% Containers for enstrophy and observables
lam_mean_spatime_enstrophy = zeros(NLAMBDA,NPARCELLS,Tmax);
ensspasub   = zeros(NSUBSIM,NPARCELLS);  % unperturbed
ensspasub1  = zeros(NSUBSIM,NPARCELLS);  % perturbed

%% MODEL SIMULATIONS (subject-level)
for sim = 1:NSUBSIM
    fprintf('  Simulation %d / %d\n', sim, NSUBSIM);

    wC   = G * C;
    sumC = repmat(sum(wC,2),1,2);

    %% HOPF SIMULATION (unperturbed)
    a   = -0.02*ones(NPARCELLS,2);
    xs  = zeros(Tmax,NPARCELLS);
    z   = 0.1*ones(NPARCELLS,2);  % state
    nn  = 0;

    % discard first 2000 time steps
    for t = 0:dt:2000
        suma = wC*z - sumC.*z;
        zz   = z(:,end:-1:1);
        z    = z + dt*(a.*z + zz.*omega - z.*(z.*z+zz.*zz) + suma) ...
                 + dsig*randn(NPARCELLS,2);
    end

    % actual modeling (unperturbed)
    for t = 0:dt:((Tmax-1)*TR)
        suma = wC*z - sumC.*z;
        zz   = z(:,end:-1:1);
        z    = z + dt*(a.*z + zz.*omega - z.*(z.*z+zz.*zz) + suma) ...
                 + dsig*randn(NPARCELLS,2);
        if abs(mod(t,TR)) < 0.01
            nn = nn+1;
            xs(nn,:) = z(:,1)';
        end
    end
    ts = xs';

    %% Kuramoto local parameter of the unperturbed model
    signal_filt = zeros(NPARCELLS, Tmax);
    Phases      = zeros(NPARCELLS, Tmax);

    for seed = 1:NPARCELLS
        ts(seed,:)          = detrend(ts(seed,:) - mean(ts(seed,:)));
        signal_filt(seed,:) = filtfilt(bfilt,afilt,ts(seed,:));
        Xanalytic           = hilbert(demean(signal_filt(seed,:)));
        Phases(seed,:)      = angle(Xanalytic);
    end

    for i = 1:NPARCELLS
        ilam = 1;
        for lam = LAMBDA
            w          = squeeze(C1(ilam,i,:));             % [NPARCELLS x 1]
            complex_ph = complex(cos(Phases),sin(Phases));  % [NPARCELLS x Tmax]
            enstrophy  = nansum(repmat(w,1,Tmax).*complex_ph) / sum(w);
            lam_mean_spatime_enstrophy(ilam,i,:) = abs(enstrophy);
            ilam = ilam+1;
        end
    end

    Rspatime          = squeeze(lam_mean_spatime_enstrophy(indsca,:,:));
    ensspasub(sim,:)  = (nanmean(Rspatime,2))';

    %% PERTURBATION: random perturbation of a for each node
    a = -0.02 + 0.02*repmat(rand(NPARCELLS,1),1,2).*ones(NPARCELLS,2);
    nn = 0;

    for t = 0:dt:((Tmax-1)*TR)
        suma = wC*z - sumC.*z;
        zz   = z(:,end:-1:1);
        z    = z + dt*(a.*z + zz.*omega - z.*(z.*z+zz.*zz) + suma) ...
                 + dsig*randn(NPARCELLS,2);
        if abs(mod(t,TR)) < 0.01
            nn = nn+1;
            xs(nn,:) = z(:,1)';
        end
    end
    ts = xs';

    %% Kuramoto local order parameter of the perturbed model
    signal_filt = zeros(NPARCELLS, Tmax);
    Phases      = zeros(NPARCELLS, Tmax);
    Rspatime1   = zeros(NPARCELLS,Tmax);

    for seed = 1:NPARCELLS
        ts(seed,:)          = detrend(ts(seed,:) - mean(ts(seed,:)));
        signal_filt(seed,:) = filtfilt(bfilt,afilt,ts(seed,:));
        Xanalytic           = hilbert(demean(signal_filt(seed,:)));
        Phases(seed,:)      = angle(Xanalytic);
    end

    for i = 1:NPARCELLS
        w          = squeeze(C1(indsca,i,:));             % [NPARCELLS x 1]
        complex_ph = complex(cos(Phases),sin(Phases));    % [NPARCELLS x Tmax]
        enstrophy  = nansum(repmat(w,1,Tmax).*complex_ph) / sum(w);
        Rspatime1(i,:) = abs(enstrophy);
    end

    ensspasub1(sim,:) = (nanmean(Rspatime1,2))';
end

%% Compute susceptibility and information capacity (subject-level)
diff_ens = ensspasub1 - ones(NSUBSIM,1)*nanmean(ensspasub); % NSUBSIM x NPARCELLS

infocapacity   = nanmean(nanstd(diff_ens));     % mean over parcels of trial-to-trial std
susceptibility = nanmean(nanmean(diff_ens));    % mean over parcels and trials

%% Save results
out_name = fullfile(out_dir,...
    sprintf('Wtrials_%03d_%d_%s_ADNI3_MCI_ABpos_GEC_sch1000_SUB%03d.mat', ...
            s, condition, measure, nsub));

save(out_name, 'infocapacity','susceptibility', ...
               'G','nsub','s','condition');

end
