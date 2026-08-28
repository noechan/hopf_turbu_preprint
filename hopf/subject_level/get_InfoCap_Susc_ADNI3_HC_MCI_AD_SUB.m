clear all
addpath(genpath('/path/to/ADNI3/code/HPC_Hopf_SUB_DTI_1000_Staging/'));
addpath('/path/to/ADNI3/code/HPC_Hopf_SUB_DTI_1000_Staging/'); %codepath
outdir='/path/to/ADNI3/code/HPC_Hopf_SUB_DTI_1000_Staging/plots/';

baseDir = '/path/to/ADNI3/code/HPC_Hopf_SUB_DTI_1000_Staging/';  % safer than pwd
trials = 100;
cond = 1;

groups=['HC_ABneg', 'MCI_ABpos', 'AD_ABpos'];
N_HCneg = 54;
N_HCpos = 39;
N_MCIpos = 33;
N_ADpos  = 26;

HCneg = load_condition_InfoCap_ADNI3_HC_ABneg_SUB(trials, cond, N_HCneg, baseDir);
HCpos = load_condition_InfoCap_ADNI3_HC_ABpos_SUB(trials, cond, N_HCpos, baseDir);
MCIpos = load_condition_InfoCap_ADNI3_MCI_ABpos_SUB(trials, cond, N_MCIpos, baseDir);
ADpos  = load_condition_InfoCap_ADNI3_AD_ABpos_SUB(trials, cond, N_ADpos, baseDir);


%% To export for further stats in R 
HC_ABneg_infocap_SUB=HCneg.infocap_subj_mean;
HC_ABpos_infocap_SUB=HCpos.infocap_subj_mean;
MCI_ABpos_infocap_SUB=MCIpos.infocap_subj_mean;
AD_ABpos_infocap_SUB=ADpos.infocap_subj_mean;

HC_ABneg_suscep_SUB=HCneg.suscep_subj_mean;
HC_ABpos_suscep_SUB=HCpos.suscep_subj_mean;
MCI_ABpos_suscep_SUB=MCIpos.suscep_subj_mean;
AD_ABpos_suscep_SUB=ADpos.suscep_subj_mean;

cd(outdir)
save ('perturbation_HC_ABneg_SUB.mat', 'HC_ABneg_infocap_SUB','HC_ABneg_suscep_SUB')
save ('perturbation_HC_ABpos_SUB.mat', 'HC_ABpos_infocap_SUB','HC_ABpos_suscep_SUB')
save ('perturbation_MCI_ABpos_SUB.mat', 'MCI_ABpos_infocap_SUB','MCI_ABpos_suscep_SUB')
save ('perturbation_AD_ABpos_SUB.mat', 'AD_ABpos_infocap_SUB','AD_ABpos_suscep_SUB')


%% Visualization

%Information Encoding Capacity

grpLabels = {'HC_ABneg','HC_ABpos','MCI_ABpos','AD_ABpos'};   % exactly 3 labels

x1 = HCneg.infocap_subj_mean; x1 = x1(~isnan(x1));
x2 = HCpos.infocap_subj_mean; x2 = x2(~isnan(x2));
x3 = MCIpos.infocap_subj_mean; x3 = x3(~isnan(x3));
x4 = ADpos.infocap_subj_mean;  x4 = x4(~isnan(x4));

data = [x1; x2; x3; x4];
groupIndex = [ ...
    ones(numel(x1),1);
    2*ones(numel(x2),1);
    3*ones(numel(x3),1);
    4*ones(numel(x4),1)];

figure
boxplot(data, groupIndex, 'Labels', grpLabels)
title('Information capacity')
ylabel('Information capacity (subject mean across 100 perturbations)')
box off; grid on


% pairwise permutation tests (SUBJECT LEVEL)
[pval,~,~] = permutationTest(x1,x2,1000); p1P(1)=pval;
[pval,~,~] = permutationTest(x1,x3,1000); p1P(2)=pval;
[pval,~,~] = permutationTest(x1,x4,1000); p1P(3)=pval;
[pval,~,~] = permutationTest(x2,x3,1000); p1P(4)=pval;
[pval,~,~] = permutationTest(x3,x4,1000); p1P(5)=pval;

% significance bars
H = sigstar({[1,2],[1,3],[1,4],[2,3],[3,4]}, p1P);
print('Info Cap ADNI3 N=238 Permutation 1000 both','-dpng')
saveas(gcf,'Info Cap ADNI3 N=238 Permutation 1000 both','fig')
close()

% FDR correction
[h_t,~,~,adj_p_t] = fdr_bh(p1P(1,1:5),0.05,'pdep','yes'); 
save('InfoCap_rejectedH0.mat','h_t','adj_p_t')


% Susceptibility

y1 = HCneg.suscep_subj_mean; y1 = y1(~isnan(y1));
y2 = HCpos.suscep_subj_mean; y2 = y2(~isnan(y2));
y3 = MCIpos.suscep_subj_mean; y3 = y3(~isnan(y3));
y4 = ADpos.suscep_subj_mean;  y4 = y4(~isnan(y4));

dataS = [y1; y2; y3; y4];
groupIndexS = [ ...
    ones(numel(y1),1);
    2*ones(numel(y2),1);
    3*ones(numel(y3),1);
    4*ones(numel(y4),1)];

figure
boxplot(dataS, groupIndexS, 'Labels', grpLabels)
title('Susceptibility')
ylabel('Susceptibility (subject mean across 100 perturbations)')
box off; grid on

% pairwise permutation tests (SUBJECT LEVEL)
[pval,~,~] = permutationTest(y1,y2,1000); p2P(1)=pval;
[pval,~,~] = permutationTest(y1,y3,1000); p2P(2)=pval;
[pval,~,~] = permutationTest(y1,y4,1000); p2P(3)=pval;
[pval,~,~] = permutationTest(y2,y3,1000); p2P(4)=pval;
[pval,~,~] = permutationTest(y3,y4,1000); p2P(5)=pval;

% significance bars
H = sigstar({[1,2],[1,3],[1,4],[2,3],[3,4]}, p2P);

print('Suscep ADNI3 N=238 Permutation 1000 both','-dpng')
saveas(gcf,'Suscep ADNI3 N=238 Permutation 1000 both','fig')
close()

% FDR correction
[h_t,~,~,adj_p_t] = fdr_bh(p2P(1,1:5),0.05,'pdep','yes'); 
save('Suscep_rejectedH0.mat','h_t','adj_p_t')

