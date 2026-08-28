%function Cnew=global_GEC_brave(cnd)

clear all
addpath(genpath('/path/to/ADNI3/code/GEC/ADNI3'));

load ('tseries_ADNI3_HC_MPRAGE_IRFSPGR_sch1000_N238rev.mat');
load('results_f_diff_fce_cond1_ADNI3_HC_sch1000.mat');
load('schaefer_MK.mat')

%% Load ABeta status table
abeta_table = readtable('/path/to/ADNI3/participants/cross_sectional/ADNI3_N238rev_with_ABETA_Status_CL24.xlsx'); % Adjust path
abeta_table.Properties.VariableNames

%% Load PTIDs
ptid_path='/path/to/ADNI3/code/turbulence_fulldenoising/sch400_N238rev';
load(fullfile(ptid_path, 'PTID_ADNI3_HC_MPRAGE_IRFSPGR_all.mat'))     % Loads variable: PTID
PTID_HC = PTID;

% Convert to cell array of strings if necessary
PTID_HC = cellstr(PTID_HC);

%% Match ABeta status for each group
% Helper function to generate masks by group and Abeta status
get_abeta_masks = @(group_ids, group_name) deal( ...
    ismember(group_ids, abeta_table.PTID(strcmp(abeta_table.GROUP, group_name) & abeta_table.ABeta_pvc == 1)), ... % Aβ+
    ismember(group_ids, abeta_table.PTID(strcmp(abeta_table.GROUP, group_name) & abeta_table.ABeta_pvc == 0)) ...  % Aβ−
    );

% Apply to each group
[isAbetaPos_HC, isAbetaNeg_HC] = get_abeta_masks(PTID_HC, 'HC');
tseries_HC_ABpos=tseries(isAbetaPos_HC,:);

cnd=1;
NSUB=size(tseries_HC_ABpos,1);
N=1000;
NPARCELLS=N;


Tau=1;%fit time-lagged statistics at 1 TR
sigma=0.01; %(additive) noise scale in the linear Hopf model

epsFC=0.00018; % learning rate for FC error
epsFCtau=0.00036; % learning rate for lagged-covariance error
maxC=0.1; %used to re-normalize the connectivity after each update, keeping its largest entry at 0.1

% Parameters of the data
TR=3;  % Repetition Time (seconds)

% Bandpass filter settings
fnq=1/(2*TR);                 % Nyquist frequency
flp = 0.008;                    % lowpass frequency of filter (Hz)
fhi = 0.08;                    % highpass
Wn=[flp/fnq fhi/fnq];         % butterworth bandpass non-dimensional frequency
k=2;                          % 2nd order butterworth filter
[bfilt,afilt]=butter(k,Wn);   % construct the filter

% Structural Connectivity
C = sc_schaefer; 
C = C/max(max(C))*0.2; %scale SC to have max weight 0.2 (initialization)

for nsub=1:NSUB
    ts=tseries_HC_ABpos{nsub,cnd}; 

    % Preprocessing signal
    clear signal_filt;
    for seed=1:N
        ts(seed,:)=detrend(ts(seed,:)-nanmean(ts(seed,:)));
        ts(isnan(ts))=0;ts(isinf(ts))=0;
        signal_filt(seed,:) =filtfilt(bfilt,afilt,ts(seed,:));
    end
    ts2=signal_filt(:,:);
    Tm=size(ts2,2);

    % Calculate FC and Covariance matrix

    FCemp=corrcoef(ts2');
    FC(nsub,:,:)=FCemp;
    COVemp=cov(ts2');
    % Compute time-lagged covariance at lag Tau
    tst=ts2';
    for i=1:N
        for j=1:N
            sigratio(i,j)=1/sqrt(COVemp(i,i))/sqrt(COVemp(j,j));
            [clag lags] = xcov(tst(:,i),tst(:,j),Tau); %cross-covariance
            indx=find(lags==Tau);
            COVtauemp(i,j)=clag(indx)/size(tst,1); % unbiased estimate
        end
    end
    COVtauemp=COVtauemp.*sigratio; % variance normalization
    COVtau(nsub,:,:)=COVtauemp;
end

FCemp=squeeze(nanmean(FC));
COVtauemp=squeeze(nanmean(COVtau));


Cnew= C;
olderror=100000;
for iter=1:5000
    % Linear Hopf FC
    [FCsim,COVsim,COVsimtotal,A]=hopf_int(Cnew,f_diff(:,:),sigma);
    % FCsim model predicted FC 
    % COVsim model predicted lagged covariance
    % COVsimtotal model covariance in the augmented state space used by the
    % integrator
    % A the linear system matrix of the Hopf approximation
    COVtausim=expm((Tau*TR)*A)*COVsimtotal; %lagged covariance prediction
    COVtausim=COVtausim(1:N,1:N);
    for i=1:N
        for j=1:N
            sigratiosim(i,j)=1/sqrt(COVsim(i,i))/sqrt(COVsim(j,j));
        end
    end
    COVtausim=COVtausim.*sigratiosim;
    errorFC(iter)=mean(mean((FCemp-FCsim).^2));
    errorCOVtau(iter)=mean(mean((COVtauemp-COVtausim).^2));
    % stopping rule
    if mod(iter,100)<0.1
        errornow=mean(mean((FCemp-FCsim).^2))+mean(mean((COVtauemp-COVtausim).^2));
        if  (olderror-errornow)/errornow<0.001
            break;
        end
        if  olderror<errornow
            break;
        end
        olderror=errornow;
    end
    %gradient descent-like learning
    for i=1:N  %% Optimization of structural information
        for j=1:N
            if (C(i,j)>0|| i==j+200 || j==i+200 )
                Cnew(i,j)=Cnew(i,j)+epsFC*(FCemp(i,j)-FCsim(i,j)) ...
                    +epsFCtau*(COVtauemp(i,j)-COVtausim(i,j));
                if Cnew(i,j)<0
                    Cnew(i,j)=0;
                end
            end
        end
    end
    Cnew = Cnew/max(max(Cnew))*maxC;
end
Ceffgroup=Cnew;

cd('/path/to/ADNI3/code/GEC/ADNI3/sch1000/inputs/')

save(sprintf('ADNI3_HC_ABpos_GECglo_cnd_%03d_sch1000_mod_eps_1_2_ratio_SC_max02.mat',cnd),'Ceffgroup');
