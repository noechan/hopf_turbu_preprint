function hopf_DTI_Grange_MCI_ABpos_GEC_SUB(nsub, s, condition)
% SUBJECT-LEVEL Hopf model (MCI Aβ−, Schaefer-1000)
%
% nsub      : subject index within MCI Aβ− subgroup (1..Nneg)
% s         : index into G_range (global coupling, G_range(s))
% condition : kept for compatibility (typically 1)
%
% Uses:
%   - f_diff_SUB_ABpos(nsub,:)            : subject-specific dominant frequency
%   - CorrFcn_SUB_ABpos(nsub,:,:)         : subject-specific empirical corrfcn(r)
%   - Ceff from GEC_SUB_%03d_ADNI3_MCI_ABpos_GEC.mat : subject-specific EC
%
% Saves:
%   WG_nsub_s_condition_ADNI3_MCI_ABpos_GEC_sch1000.mat

%% Paths and basic setup

this_file  = mfilename('fullpath');
code_dir   = fileparts(this_file);
addpath(genpath(code_dir))
out_dir = fullfile(code_dir, 'results');

%% Load subject-level f_diff and empirical corrfcn (MCI Aβ−)
freq_file =fullfile(code_dir,'data','results_f_diff_fce_cond1_ADNI3_MCI_ABpos_sch1000_SUB.mat');
% Expected: f_diff_SUB_ABpos [Nneg x NPARCELLS]
load(freq_file, 'f_diff_SUB_ABpos');

spacorr_file =fullfile(code_dir,'data','empirical_spacorr_rest_cond_1_ADNI3_MCI_ABpos_sch1000_SUB.mat');
% Expected: CorrFcn_SUB_ABpos [Nneg x NPARCELLS x NR], rr [NPARCELLS x NPARCELLS], xrange [1 x NR]
load(spacorr_file, 'CorrFcn_SUB_ABpos', 'rr', 'xrange');

[NSUB_ABpos, NPARCELLS, NR] = size(CorrFcn_SUB_ABpos);
if nsub < 1 || nsub > NSUB_ABpos
    error('nsub=%d is out of range for MCI Aβ+ subjects (1..%d).', nsub, NSUB_ABpos);
end

% Empirical corrfcn for this subject: [NPARCELLS x NR]
empcorrfcn = squeeze(CorrFcn_SUB_ABpos(nsub,:,:));

%% Load subject-specific effective connectivity
gec_file = fullfile(code_dir,'data','GEC_SUB_ADNI3_MCI_ABpos',sprintf('GEC_SUB_%03d_ADNI3_MCI_ABpos_GEC.mat', nsub));
load(gec_file, 'Ceff');   % [NPARCELLS x NPARCELLS]

%% Parameters and G range
NRini   = 20;
NRfin   = 380;
NSUBSIM = 100;
lambda  = 0.175;
lambda  = round(lambda,2);

G_range = 0:0.01:0.4;
G = G_range(s);

fprintf('MCI Aβ+ subject nsub=%d | G index s=%d | G=%.2f\n', nsub, s, G);

%% Geometric kernel C1 (for enstrophy)
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
f_sub = squeeze(f_diff_SUB_ABpos(nsub,:));   % ensure vector [NPARCELLS x 1]
omega = repmat(2*pi*f_sub(:),1,2); 
omega(:,1) = -omega(:,1);

dt   = 0.1*TR/2;
sig  = 0.01;
dsig = sqrt(dt)*sig;

%% Containers
corrfcn                    = zeros(NPARCELLS,NR);
lam_mean_spatime_enstrophy = zeros(NLAMBDA,NPARCELLS,Tmax);
err_hete                   = zeros(1,NSUBSIM);
InfoFlow                   = zeros(NSUBSIM,NLAMBDA-1);
InfoCascade                = zeros(1,NSUBSIM);
Turbulence                 = zeros(1,NSUBSIM);
mutinfo                    = zeros(1,NLAMBDA);

%% Connectivity scaling
factor = max(Ceff(:));
C      = Ceff / factor * 0.2;

%% Distance bin width (consistent with rr/xrange)
range_dist = max(rr(:));
delta      = range_dist / NR;

%% MODEL SIMULATIONS
for sim = 1:NSUBSIM
    fprintf('  Simulation %d / %d\n', sim, NSUBSIM);

    wC   = G * C;
    sumC = repmat(sum(wC,2),1,2);

    %% Hopf Simulation
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

    % actual modelling
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

    %% Simulated FC and phases
    signal_filt = zeros(NPARCELLS, Tmax);
    Phases      = zeros(NPARCELLS, Tmax);
    for seed = 1:NPARCELLS
        sig_row           = detrend(ts(seed,:) - mean(ts(seed,:)));
        signal_filt(seed,:) = filtfilt(bfilt,afilt,sig_row);
        Xanalytic         = hilbert(demean(signal_filt(seed,:)));
        Phases(seed,:)    = angle(Xanalytic);
    end
    fcsimul = corrcoef(signal_filt');

    %% Simulated spatial corrfcn with explicit NaNs for empty bins
    for i = 1:NPARCELLS
        numind    = zeros(1,NR);
        corrfcn_1 = zeros(1,NR);
        for j = 1:NPARCELLS
            r     = rr(i,j);
            index = floor(r/delta)+1;
            if index > NR
                index = NR;
            end
            mcc = fcsimul(i,j);
            if ~isnan(mcc)
                corrfcn_1(index) = corrfcn_1(index) + mcc;
                numind(index)    = numind(index)    + 1;
            end
        end
        row        = nan(1,NR);
        valid_bins = numind > 0;
        row(valid_bins) = corrfcn_1(valid_bins)./numind(valid_bins);
        corrfcn(i,:) = row;

        % Enstrophy
        ilam = 1;
        for lam = LAMBDA
            w          = squeeze(C1(ilam,i,:));              % [NPARCELLS x 1]
            complex_ph = complex(cos(Phases),sin(Phases));   % [NPARCELLS x Tmax]
            enstrophy  = nansum(repmat(w,1,Tmax).*complex_ph) / sum(w);
            lam_mean_spatime_enstrophy(ilam,i,:) = abs(enstrophy);
            ilam = ilam+1;
        end
    end

    Rspatime = squeeze(lam_mean_spatime_enstrophy(indsca,:,:));
    Rsub     = nanstd(Rspatime(:));

    %% Mutual information across scales
    mutinfo(:) = 0;
    for ilam = 2:NLAMBDA
        A = squeeze(lam_mean_spatime_enstrophy(ilam,:,2:end))';
        B = squeeze(lam_mean_spatime_enstrophy(ilam-1,:,1:end-1))';
        [cc, pp] = corr(A,B);
        mutinfo(ilam) = nanmean(cc(pp(:)<0.05));
    end

    %% Error between simulated and empirical corrfcn
    err1 = nan(NPARCELLS,1);
    for i = 1:NPARCELLS
        sim_row = corrfcn(i,NRini:NRfin);
        emp_row = empcorrfcn(i,NRini:NRfin);

        mask = ~isnan(sim_row) & ~isnan(emp_row);
        if any(mask)
            diff_sq  = (sim_row(mask) - emp_row(mask)).^2;
            err1(i)  = mean(diff_sq);
        else
            err1(i)  = NaN;
        end
    end
    err_hete(sim) = sqrt(nanmean(err1));

    %% Observables
    Inflam2          = mutinfo(2:end);
    InfoFlow(sim,:)  = Inflam2;
    InfoCascade(sim) = nanmean(Inflam2);
    Turbulence(sim)  = Rsub;
end

%% Save results
out_name = fullfile(out_dir,'G_range','MCI_ABpos',sprintf('WG_%03d_%03d_%d_ADNI3_MCI_ABpos_GEC_sch1000.mat', ...
                   nsub, s, condition));
save(out_name, 'InfoCascade','InfoFlow','Turbulence','err_hete', ...
               'G','nsub','s','condition');
end
