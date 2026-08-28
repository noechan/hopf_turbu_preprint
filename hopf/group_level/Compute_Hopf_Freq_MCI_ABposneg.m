%% Whole-brain modelling in the homogeneous case

%This model is based on Stuard-Landau oscillators, the local dynamics of single brain
%regions using the normal form of a Hopf bifurcation.
%The dynamics of the N brain regions were coupled through the connectivity matrix,
%which was given by the connectome of healthy subjects (C).
%coupling among areas is given by the SC-> by the g scaling parameter.
%For the homogeneous case, in which we seta=0 for all nodes. 
%This choice was based on previous studies which suggest that the best fit
%to the empirical data arises at the brink of the Hopf bifurcation 
%where a~0 (Deco et al. 2017). So, here the only free parameter is the g.

clear all
cd('/path/to/ADNI3/code/HPC_Hopf_DTI_1000_Staging/');
load('tseries_ADNI3_MCI_MPRAGE_IRFSPGR_sch1000_N238rev.mat')              % Denoised time-series           % Denoised time-series
addpath(genpath('/path/to/ADNI3/code/HPC_Hopf_DTI_1000_Staging/'))

% Parameters of the data
TR=3;  % Repetition Time (seconds)
NPARCELLS=1000;
Tmax=197;
NCOND=1;

%--------------------------------------------------------------------------
%COMPUTE POWER SPECTRA FOR
%NARROWLY FILTERED DATA WITH LOW BANDPASS (0.008 to 0.08 Hz)
%--------------------------------------------------------------------------

% Bandpass filter settings
fnq=1/(2*TR);                 % Nyquist frequency
flp = 0.008;                    % lowpass frequency of filter (Hz)
fhi = 0.08;                    % highpass
Wn=[flp/fnq fhi/fnq];         % butterworth bandpass non-dimensional frequency
k=2;                          % 2nd order butterworth filter
[bfilt,afilt]=butter(k,Wn);   % construct the filter
Isubdiag = find(tril(ones(NPARCELLS),-1));


for cond=1:NCOND
    xs=tseries(:,cond);
    NSUB=size(find(~cellfun(@isempty,xs)),1);
    fce1=zeros(NSUB,NPARCELLS,NPARCELLS);
    tss=zeros(NPARCELLS,Tmax);
    TT=Tmax;
    Ts = TT*TR;
    freq = (0:TT/2-1)/Ts;
    nfreqs=length(freq);
    PowSpect=zeros(nfreqs,NPARCELLS,NSUB);
    for sub=1:NSUB
        sub
        ts=xs{sub};
        [Ns, Tmax]=size(ts);
        TT=Tmax;
        Ts = TT*TR;
        freq = (0:TT/2-1)/Ts;
        nfreqs=length(freq);
        
        for seed=1:NPARCELLS
            x=detrend(ts(seed,:)-mean(ts(seed,:)));
            tss(seed,:)=filtfilt(bfilt,afilt,x);
            pw = abs(fft(tss(seed,:))); %Fourier transform
            PowSpect(:,seed,sub) = pw(1:floor(TT/2)).^2/(TT/TR); %calculates power spectrum
        end
        fce1(sub,:,:)=corrcoef(tss(1:NPARCELLS,:)','rows','pairwise'); 
    end
    
    fce=squeeze(mean(fce1)); %functional connectivity empirical of filtered time-series: mean for all subjects
    
    Power_Areas=squeeze(mean(PowSpect,3)); %power for each freq averaged for all subjects
    for seed=1:NPARCELLS
        Power_Areas(:,seed)=gaussfilt(freq,Power_Areas(:,seed)',0.01);
    end
    
    [maxpowdata,index]=max(Power_Areas);
    f_diff = freq(index);
    f_diff(find(f_diff==0))=mean(f_diff(find(f_diff~=0)));

    % Split by Abeta status
    %% Load ABeta status table
abeta_table = readtable('/path/to/ADNI3/participants/cross_sectional/ADNI3_N238rev_with_ABETA_Status_CL24.xlsx'); % Adjust path
abeta_table.Properties.VariableNames

%% Load PTIDs
ptid_path='/path/to/ADNI3/code/turbulence_fulldenoising/sch400_N238rev';
load(fullfile(ptid_path, 'PTID_ADNI3_MCI_MPRAGE_IRFSPGR_all.mat'))     % Loads variable: PTID
PTID_MCI = PTID;

% Convert to cell array of strings if necessary
PTID_MCI = cellstr(PTID_MCI);

%% Match ABeta status for each group
% Helper function to generate masks by group and Abeta status
get_abeta_masks = @(group_ids, group_name) deal( ...
    ismember(group_ids, abeta_table.PTID(strcmp(abeta_table.GROUP, group_name) & abeta_table.ABeta_pvc == 1)), ... % Aβ+
    ismember(group_ids, abeta_table.PTID(strcmp(abeta_table.GROUP, group_name) & abeta_table.ABeta_pvc == 0)) ...  % Aβ−
    );

% Apply to each group
[isAbetaPos_MCI, isAbetaNeg_MCI] = get_abeta_masks(PTID_MCI, 'MCI');
% fce by AB status
fce_MCI_ABpos=squeeze(mean(fce1(isAbetaPos_MCI,:,:)));
fce_MCI_ABneg=squeeze(mean(fce1(isAbetaNeg_MCI,:,:)));
% freq by AB status
 Power_Areas_MCI_ABpos=squeeze(mean(PowSpect(:,:,isAbetaPos_MCI),3)); %power for each freq averaged for all subjects
    for seed=1:NPARCELLS
        Power_Areas_MCI_ABpos(:,seed)=gaussfilt(freq,Power_Areas_MCI_ABpos(:,seed)',0.01);
    end
    
    [maxpowdata_MCI_ABpos,index_MCI_ABpos]=max(Power_Areas_MCI_ABpos);
    f_diff_MCI_ABpos = freq(index_MCI_ABpos);
    f_diff_MCI_ABpos(find(f_diff_MCI_ABpos==0))=mean(f_diff_MCI_ABpos(find(f_diff_MCI_ABpos~=0)));

     Power_Areas_MCI_ABneg=squeeze(mean(PowSpect(:,:,isAbetaNeg_MCI),3)); %power for each freq averaged for all subjects
    for seed=1:NPARCELLS
        Power_Areas_MCI_ABneg(:,seed)=gaussfilt(freq,Power_Areas_MCI_ABneg(:,seed)',0.01);
    end
    
    [maxpowdata_MCI_ABneg,index_MCI_ABneg]=max(Power_Areas_MCI_ABneg);
    f_diff_MCI_ABneg = freq(index_MCI_ABneg);
    f_diff_MCI_ABneg(find(f_diff_MCI_ABneg==0))=mean(f_diff_MCI_ABneg(find(f_diff_MCI_ABneg~=0)));

    
    cd('/path/to/ADNI3/code/HPC_Hopf_DTI_1000_Staging/')
    save (sprintf('results_f_diff_fce_cond%d_ADNI3_MCI_sch1000.mat', cond),...
    'f_diff', 'fce');
    save (sprintf('results_f_diff_fce_cond%d_ADNI3_MCI_ABpos_sch1000.mat', cond),...
    'f_diff_MCI_ABpos', 'fce_MCI_ABpos');
    save (sprintf('results_f_diff_fce_cond%d_ADNI3_MCI_ABneg_sch1000.mat', cond),...
    'f_diff_MCI_ABneg', 'fce_MCI_ABneg');
end

