clear all
cd('/path/to/ADNI3/code/HPC_Hopf_DTI_1000_Staging/');
load('tseries_ADNI3_HC_MPRAGE_IRFSPGR_sch1000_N238rev.mat')               % Denoised time-series           % Denoised time-series
addpath(genpath('/path/to/ADNI3/code/HPC_Hopf_DTI_1000_Staging/'))
load SchaeferCOG.mat;

% Parameters of the data
TR=3;  % Repetition Time (seconds)
NPARCELLS=1000;
NR=400; %inertial subrange
NCOND=1; % put here your total conditions
Tmax=197;

% Bandpass filter settings
fnq=1/(2*TR);                 % Nyquist frequency
flp = 0.008;                    % lowpass frequency of filter (Hz)
fhi = 0.08;                    % highpass
Wn=[flp/fnq fhi/fnq];         % butterworth bandpass non-dimensional frequency
k=2;                          % 2nd order butterworth filter
[bfilt,afilt]=butter(k,Wn);   % construct the filter

%----------------------------------------------------------------------------------------------------------------------
%COMPUTE THE EMPIRICAL FC AS THE SPATIAL CORRELATIONS OF TWO POINTS
%SEPARATED BY A EUCLIDEAN DISTANCE r WITHIN THE INERTIAL SUBRANGE
%----------------------------------------------------------------------------------------------------------------------

for i=1:NPARCELLS
    for j=1:NPARCELLS
        rr(i,j)=norm(SchaeferCOG(i,:)-SchaeferCOG(j,:));
    end
end
range=max(max(rr));
delta=range/NR;

for i=1:NR
    xrange(i)=delta/2+delta*(i-1);
end


for cond=1:NCOND
    xs=tseries(:,cond);
    NSUB=size(find(~cellfun(@isempty,xs)),1);
    ensspasub=zeros(NSUB,NPARCELLS);
    DTspatime=zeros(NPARCELLS,Tmax);
    Rsub=zeros(1,NSUB);
    DTsub=zeros(1,NSUB);
    corrfcnsub=zeros(NSUB,NPARCELLS,NR);
    for sub=1:NSUB
        sub
        ts=xs{sub};
        
        for seed=1:NPARCELLS
            ts(seed,:)=detrend(ts(seed,:)-mean(ts(seed,:)));
            signal_filt(seed,:) =filtfilt(bfilt,afilt,ts(seed,:));
        end
        
        fce=corrcoef(signal_filt');
        
        for i=1:NPARCELLS
            numind=zeros(1,NR);
            corrfcn_1=zeros(1,NR);
            for j=1:NPARCELLS
                r=rr(i,j);
                index=floor(r/delta)+1;
                if index==NR+1
                    index=NR;
                end
                mcc=fce(i,j);
                if ~isnan(mcc)
                    corrfcn_1(index)=corrfcn_1(index)+mcc;
                    numind(index)=numind(index)+1;
                end
            end
            corrfcnsub(sub,i,:)=corrfcn_1./numind;
        end
    end

    for sub=1:NSUB
        sub
        ts=xs{sub};
        
        for seed=1:NPARCELLS
            ts(seed,:)=detrend(ts(seed,:)-mean(ts(seed,:)));
            signal_filt(seed,:) =filtfilt(bfilt,afilt,ts(seed,:));
        end
        
        fce=corrcoef(signal_filt');
        
        for i=1:NPARCELLS
            numind=zeros(1,NR);
            corrfcn_1=zeros(1,NR);
            for j=1:NPARCELLS
                r=rr(i,j);
                index=floor(r/delta)+1;
                if index==NR+1
                    index=NR;
                end
                mcc=fce(i,j);
                if ~isnan(mcc)
                    corrfcn_1(index)=corrfcn_1(index)+mcc;
                    numind(index)=numind(index)+1;
                end
            end
            corrfcnsub(sub,i,:)=corrfcn_1./numind;
        end
    end
    
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
% --- Split corrfcn per subject by Aβ status---
CorrFcn_sub_HC_ABpos = corrfcnsub(isAbetaPos_HC,:,:);
CorrFcn_sub_HC_ABneg = corrfcnsub(isAbetaNeg_HC,:,:);

% Group
    
    corrfcn=squeeze(nanmean(corrfcnsub));
    corrfcn_HC_ABpos = squeeze(nanmean(CorrFcn_sub_HC_ABpos));
    corrfcn_HC_ABneg = squeeze(nanmean(CorrFcn_sub_HC_ABneg));

    
    cd('/path/to/ADNI3/code/HPC_Hopf_DTI_1000_Staging/')
    save (sprintf('empirical_spacorr_rest_cond_%d_ADNI3_HC_sch1000.mat', cond),'corrfcn', 'corrfcnsub');
    save (sprintf('empirical_spacorr_rest_cond_%d_ADNI3_HC_ABpos_sch1000.mat', cond),'corrfcn_HC_ABpos');
    save (sprintf('empirical_spacorr_rest_cond_%d_ADNI3_HC_ABneg_sch1000.mat', cond),'corrfcn_HC_ABneg');
end